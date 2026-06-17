"""Typer-based CLI for omics_target_prioritization.

This is the ``[project.scripts]`` entrypoint. Every pipeline stage is reachable
from the command line, so the repo is usable without writing Python:

    omics-target-prioritization simulate  --out results/locus.json
    omics-target-prioritization integrate --locus results/locus.json --out results/evidence.json
    omics-target-prioritization score     --evidence results/evidence.json --out results/scores.json
    omics-target-prioritization report     --scores results/scores.json --out results/dossier.html
    omics-target-prioritization run-all    --out results/
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import typer
from loguru import logger
from rich.console import Console
from rich.table import Table

from omics_target_prioritization import __version__
from omics_target_prioritization.integrate.v2g import integrate_locus
from omics_target_prioritization.models import (
    EvidenceItem,
    Gene,
    TargetScore,
)
from omics_target_prioritization.report.dossier import write_dossier
from omics_target_prioritization.score.prioritize import prioritize
from omics_target_prioritization.simulate import SimConfig, SimulatedLocus, simulate_locus

app = typer.Typer(
    name="omics-target-prioritization",
    help="omics-target-prioritization — GWAS x multi-omics QTL colocalization to target scores. See README.md.",
    no_args_is_help=True,
)
console = Console()


def _write_json(path: Path, payload: dict[str, object]) -> None:
    """Write ``payload`` as pretty JSON, creating parent directories."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


@app.command()
def version() -> None:
    """Print the installed version and exit."""
    console.print(f"omics-target-prioritization v{__version__}")


@app.command()
def simulate(
    out: Annotated[Path, typer.Option(help="Output JSON path for the locus.")] = Path(
        "results/locus.json"
    ),
    seed: Annotated[int, typer.Option(help="RNG seed for reproducibility.")] = 7,
    n_neighbour_genes: Annotated[int, typer.Option(help="Neighbouring genes.")] = 3,
) -> None:
    """Simulate a GWAS locus with planted multi-omics QTL evidence."""
    cfg = SimConfig(seed=seed, n_neighbour_genes=n_neighbour_genes)
    sim = simulate_locus(cfg)
    _write_json(out, json.loads(sim.model_dump_json()))
    logger.info("Wrote simulated locus to {} (causal gene={})", out, sim.causal_gene_id)
    console.print(f"[green]Simulated[/] locus -> {out}  (causal: {sim.causal_gene_id})")


@app.command()
def integrate(
    locus: Annotated[Path, typer.Option(help="Input locus JSON (from simulate).")] = Path(
        "results/locus.json"
    ),
    out: Annotated[Path, typer.Option(help="Output evidence JSON path.")] = Path(
        "results/evidence.json"
    ),
) -> None:
    """Run colocalization-based variant-to-gene integration over a locus."""
    sim = SimulatedLocus.model_validate_json(locus.read_text(encoding="utf-8"))
    evidence = integrate_locus(sim.gwas, sim.genes, sim.qtls)
    payload: dict[str, object] = {
        "locus_id": sim.gwas.locus_id,
        "genes": [g.model_dump() for g in sim.genes],
        "evidence": {
            gid: [item.model_dump(mode="json") for item in items] for gid, items in evidence.items()
        },
    }
    _write_json(out, payload)
    n_items = sum(len(v) for v in evidence.values())
    logger.info("Integrated {} genes, {} evidence items -> {}", len(sim.genes), n_items, out)
    console.print(f"[green]Integrated[/] {len(sim.genes)} genes -> {out}")


@app.command()
def score(
    evidence: Annotated[Path, typer.Option(help="Input evidence JSON (from integrate).")] = Path(
        "results/evidence.json"
    ),
    out: Annotated[Path, typer.Option(help="Output scores JSON path.")] = Path(
        "results/scores.json"
    ),
) -> None:
    """Aggregate evidence into ranked, calibrated target scores."""
    payload = json.loads(evidence.read_text(encoding="utf-8"))
    genes = [Gene.model_validate(g) for g in payload["genes"]]
    ev_map: dict[str, list[EvidenceItem]] = {
        gid: [EvidenceItem.model_validate(it) for it in items]
        for gid, items in payload["evidence"].items()
    }
    ranked = prioritize(genes, ev_map)
    out_payload: dict[str, object] = {
        "locus_id": payload.get("locus_id", ""),
        "scores": [s.model_dump(mode="json") for s in ranked],
    }
    _write_json(out, out_payload)
    _print_ranking(ranked)
    logger.info("Scored {} targets -> {}", len(ranked), out)


@app.command()
def report(
    scores: Annotated[Path, typer.Option(help="Input scores JSON (from score).")] = Path(
        "results/scores.json"
    ),
    out: Annotated[Path, typer.Option(help="Output HTML dossier for the top target.")] = Path(
        "results/dossier.html"
    ),
    rank: Annotated[int, typer.Option(help="1-based rank to render (default top).")] = 1,
) -> None:
    """Render an HTML dossier for a prioritized target."""
    payload = json.loads(scores.read_text(encoding="utf-8"))
    ranked = [TargetScore.model_validate(s) for s in payload["scores"]]
    if not ranked:
        raise typer.BadParameter("No scores to report.")
    target = ranked[min(rank - 1, len(ranked) - 1)]
    write_dossier(target, out, locus_id=str(payload.get("locus_id", "")))
    logger.info("Wrote dossier for {} -> {}", target.symbol, out)
    console.print(f"[green]Dossier[/] for {target.symbol} -> {out}")


@app.command(name="run-all")
def run_all(
    out: Annotated[Path, typer.Option(help="Output directory.")] = Path("results"),
    seed: Annotated[int, typer.Option(help="RNG seed.")] = 7,
) -> None:
    """Run the full pipeline: simulate -> integrate -> score -> report."""
    out.mkdir(parents=True, exist_ok=True)
    sim = simulate_locus(SimConfig(seed=seed))
    _write_json(out / "locus.json", json.loads(sim.model_dump_json()))

    evidence = integrate_locus(sim.gwas, sim.genes, sim.qtls)
    ranked = prioritize(sim.genes, evidence)
    _write_json(
        out / "scores.json",
        {
            "locus_id": sim.gwas.locus_id,
            "scores": [s.model_dump(mode="json") for s in ranked],
        },
    )

    dossier_path = write_dossier(ranked[0], out / "dossier.html", locus_id=sim.gwas.locus_id)
    _print_ranking(ranked)
    console.print(f"[green]run-all complete[/] · dossier -> {dossier_path}")
    logger.info("run-all complete; top target={}", ranked[0].symbol)


def _print_ranking(ranked: list[TargetScore]) -> None:
    """Pretty-print a ranked TargetScore list as a rich table."""
    table = Table(title="Prioritized targets")
    table.add_column("Rank", justify="right")
    table.add_column("Symbol")
    table.add_column("Score", justify="right")
    table.add_column("Confidence")
    table.add_column("Layers", justify="right")
    table.add_column("max PP.H4", justify="right")
    for i, s in enumerate(ranked, start=1):
        table.add_row(
            str(i),
            s.symbol,
            f"{s.total:.3f}",
            s.confidence,
            str(s.n_layers),
            f"{s.max_h4:.3f}",
        )
    console.print(table)


if __name__ == "__main__":
    app()
