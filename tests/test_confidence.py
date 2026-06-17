"""Tests for the calibrated confidence labelling."""

from __future__ import annotations

from omics_target_prioritization.evidence.confidence import (
    assign_confidence,
    count_corroborating_layers,
)
from omics_target_prioritization.models import EvidenceItem, Provenance


def _item(source: str, score: float) -> EvidenceItem:
    return EvidenceItem(
        gene_id="G",
        source=source,  # type: ignore[arg-type]
        score=score,
        provenance=Provenance(dataset="d", method="m"),
    )


def test_count_corroborating_layers() -> None:
    items = [
        _item("eQTL_coloc", 0.9),
        _item("pQTL_coloc", 0.8),
        _item("sQTL_coloc", 0.2),  # below threshold, not counted
        _item("distance_to_tss", 0.9),  # not a coloc layer
    ]
    assert count_corroborating_layers(items) == 2


def test_two_layers_high_score_is_high_confidence() -> None:
    items = [_item("eQTL_coloc", 0.95), _item("pQTL_coloc", 0.9)]
    label, n_layers, best_h4 = assign_confidence(items, total=0.8)
    assert label == "high"
    assert n_layers == 2
    assert best_h4 == 0.95


def test_single_strong_layer_is_medium() -> None:
    items = [_item("eQTL_coloc", 0.9)]
    label, n_layers, _ = assign_confidence(items, total=0.45)
    assert label == "medium"
    assert n_layers == 1


def test_distance_only_is_low_confidence() -> None:
    items = [_item("distance_to_tss", 1.0)]
    label, n_layers, best_h4 = assign_confidence(items, total=0.2)
    assert label == "low"
    assert n_layers == 0
    assert best_h4 == 0.0


def test_confidence_rises_with_more_layers() -> None:
    one = assign_confidence([_item("eQTL_coloc", 0.7)], total=0.4)[0]
    two = assign_confidence([_item("eQTL_coloc", 0.7), _item("pQTL_coloc", 0.7)], total=0.6)[0]
    rank = {"low": 0, "medium": 1, "high": 2}
    assert rank[two] >= rank[one]
    assert two == "high"
