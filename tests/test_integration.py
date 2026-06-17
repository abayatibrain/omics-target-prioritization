"""End-to-end integration tests for the full pipeline."""

from __future__ import annotations

from pathlib import Path

import pytest

from omics_target_prioritization.integrate.v2g import integrate_locus
from omics_target_prioritization.report.dossier import render_dossier, write_dossier
from omics_target_prioritization.score.prioritize import prioritize
from omics_target_prioritization.simulate import SimConfig, simulate_locus


@pytest.mark.integration
def test_full_pipeline_ranks_true_causal_first_and_writes_dossier(tmp_path: Path) -> None:
    # simulate -> integrate -> score -> report
    sim = simulate_locus(SimConfig(seed=11, n_snps=150, n_neighbour_genes=4))
    evidence = integrate_locus(sim.gwas, sim.genes, sim.qtls)
    ranked = prioritize(sim.genes, evidence)

    # Known answer: the true causal gene is ranked first.
    assert ranked[0].gene_id == sim.causal_gene_id
    assert ranked[0].confidence == "high"
    assert ranked[0].total > 0.0

    # Dossier renders and is written to disk.
    out = tmp_path / "dossier.html"
    write_dossier(ranked[0], out, locus_id=sim.gwas.locus_id)
    assert out.exists()
    html = out.read_text(encoding="utf-8")
    assert ranked[0].symbol in html
    assert "PP.H4" in html
    assert "Provenance" in html


@pytest.mark.integration
def test_dossier_contains_provenance_table() -> None:
    sim = simulate_locus(SimConfig(seed=3))
    evidence = integrate_locus(sim.gwas, sim.genes, sim.qtls)
    ranked = prioritize(sim.genes, evidence)
    html = render_dossier(ranked[0], locus_id=sim.gwas.locus_id)
    # Every evidence item's method should appear in the provenance table.
    assert "coloc.abf" in html
    assert "distance_decay" in html


@pytest.mark.integration
def test_pipeline_is_deterministic() -> None:
    a = simulate_locus(SimConfig(seed=42))
    b = simulate_locus(SimConfig(seed=42))
    ra = prioritize(a.genes, integrate_locus(a.gwas, a.genes, a.qtls))
    rb = prioritize(b.genes, integrate_locus(b.gwas, b.genes, b.qtls))
    assert [s.gene_id for s in ra] == [s.gene_id for s in rb]
    assert abs(ra[0].total - rb[0].total) < 1e-12
