"""Tests for the Open-Targets harmonic-sum aggregation."""

from __future__ import annotations

from datetime import datetime

from omics_target_prioritization.evidence.aggregate import aggregate_scores, harmonic_sum
from omics_target_prioritization.models import UTC, EvidenceItem, Provenance


def _item(source: str, score: float, weight: float = 1.0) -> EvidenceItem:
    return EvidenceItem(
        gene_id="G",
        source=source,  # type: ignore[arg-type]
        score=score,
        weight=weight,
        provenance=Provenance(
            dataset="d",
            method="m",
            timestamp=datetime.now(UTC),
        ),
    )


def test_harmonic_sum_bounded_in_unit_interval() -> None:
    assert harmonic_sum([]) == 0.0
    assert 0.0 <= harmonic_sum([0.5]) <= 1.0
    # A single max source cannot saturate the score (Open-Targets behaviour:
    # the normalizer is the infinite Basel limit, not the per-call partial sum).
    assert 0.0 < harmonic_sum([1.0]) < 1.0
    # Many unit sources approach, but never exceed, 1.
    assert 0.0 < harmonic_sum([1.0] * 50) <= 1.0
    assert harmonic_sum([1.0] * 50) > harmonic_sum([1.0])
    assert 0.0 <= harmonic_sum([0.9, 0.8, 0.3, 0.1]) <= 1.0


def test_harmonic_sum_is_monotonic_in_added_evidence() -> None:
    base = harmonic_sum([0.6])
    more = harmonic_sum([0.6, 0.4])
    assert more >= base


def test_harmonic_sum_is_monotonic_in_score_magnitude() -> None:
    low = harmonic_sum([0.3, 0.2])
    high = harmonic_sum([0.9, 0.2])
    assert high > low


def test_top_signal_dominates() -> None:
    # One strong source should outscore many weak ones.
    strong = harmonic_sum([0.95])
    weak_many = harmonic_sum([0.2, 0.2, 0.2, 0.2, 0.2])
    assert strong > weak_many


def test_aggregate_takes_best_per_source() -> None:
    # Two eQTL tissues for the same source: only the strongest counts.
    items = [
        _item("eQTL_coloc", 0.95),
        _item("eQTL_coloc", 0.10),
        _item("distance_to_tss", 0.20, weight=0.3),
    ]
    total, breakdown = aggregate_scores(items)
    assert breakdown["eQTL_coloc"] == 0.95
    assert 0.0 <= total <= 1.0


def test_colocalizing_gene_beats_distance_only() -> None:
    coloc_gene = [
        _item("eQTL_coloc", 0.97),
        _item("pQTL_coloc", 0.93),
        _item("distance_to_tss", 0.9, weight=0.3),
    ]
    distance_gene = [_item("distance_to_tss", 1.0, weight=0.3)]
    total_coloc, _ = aggregate_scores(coloc_gene)
    total_dist, _ = aggregate_scores(distance_gene)
    assert total_coloc > total_dist
