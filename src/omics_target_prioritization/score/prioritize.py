"""Locus-level target prioritization: rank genes into ordered TargetScores.

This module ties the evidence layers together. Given the per-gene evidence
produced by :mod:`omics_target_prioritization.integrate.v2g`, it:

1. aggregates each gene's evidence with the Open-Targets harmonic sum
   (:mod:`omics_target_prioritization.evidence.aggregate`);
2. assigns a calibrated confidence label
   (:mod:`omics_target_prioritization.evidence.confidence`);
3. ranks the resulting :class:`~omics_target_prioritization.models.TargetScore`s
   by total score, descending.

The gene with colocalizing eQTL+pQTL evidence (high PP.H4 in multiple layers)
ranks above a gene that is merely close to the lead variant.
"""

from __future__ import annotations

from collections.abc import Mapping

from omics_target_prioritization.evidence.aggregate import aggregate_scores
from omics_target_prioritization.evidence.confidence import assign_confidence
from omics_target_prioritization.models import (
    EvidenceItem,
    Gene,
    MrDirection,
    TargetScore,
)


def _mr_summary(items: list[EvidenceItem]) -> tuple[MrDirection, float | None]:
    """Pull the strongest MR direction and p-value out of the evidence list."""
    mr_items = [it for it in items if it.source == "mr_directional"]
    if not mr_items:
        return "none", None
    # Strongest MR signal = lowest p-value recorded in provenance.
    best = min(mr_items, key=lambda it: float(it.provenance.parameters.get("pvalue", 1.0)))
    direction = str(best.provenance.parameters.get("direction", "none"))
    pvalue = best.provenance.parameters.get("pvalue")
    direction_typed: MrDirection
    if direction == "risk":
        direction_typed = "risk"
    elif direction == "protective":
        direction_typed = "protective"
    else:
        direction_typed = "none"
    return direction_typed, (float(pvalue) if pvalue is not None else None)


def score_target(gene: Gene, items: list[EvidenceItem]) -> TargetScore:
    """Aggregate and calibrate one gene's evidence into a :class:`TargetScore`.

    Parameters
    ----------
    gene
        The gene being scored (for ``gene_id`` / ``symbol``).
    items
        The gene's evidence items.

    Returns
    -------
    TargetScore
        The aggregated, calibrated, provenance-carrying score.
    """
    total, breakdown = aggregate_scores(items)
    confidence, n_layers, best_h4 = assign_confidence(items, total)
    mr_dir, mr_p = _mr_summary(items)
    return TargetScore(
        gene_id=gene.gene_id,
        symbol=gene.symbol,
        total=total,
        breakdown=breakdown,
        confidence=confidence,
        n_layers=n_layers,
        max_h4=best_h4,
        mr_direction=mr_dir,
        mr_pvalue=mr_p,
        evidence=items,
    )


def prioritize(
    genes: list[Gene],
    evidence: Mapping[str, list[EvidenceItem]],
) -> list[TargetScore]:
    """Rank all candidate genes at a locus, highest score first.

    Ties are broken deterministically by ``max_h4`` then ``gene_id`` so the
    ordering is stable across runs.

    Parameters
    ----------
    genes
        Candidate genes at the locus.
    evidence
        Mapping ``gene_id -> evidence items`` (as returned by
        :func:`omics_target_prioritization.integrate.v2g.integrate_locus`).

    Returns
    -------
    list[TargetScore]
        Target scores sorted by ``total`` descending.
    """
    scores: list[TargetScore] = []
    for gene in genes:
        items = evidence.get(gene.gene_id, [])
        scores.append(score_target(gene, items))

    scores.sort(key=lambda s: (s.total, s.max_h4, s.gene_id), reverse=True)
    return scores
