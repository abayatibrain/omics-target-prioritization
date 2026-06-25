"""Tests for the drug-target MR / directionality evidence channel."""

from __future__ import annotations

import pytest

from omics_target_prioritization.evidence.mr_evidence import (
    MR_Z_MIDPOINT,
    build_mr_evidence,
    mr_direction,
    mr_strength,
)
from omics_target_prioritization.models import (
    EvidenceItem,
    Gene,
    Provenance,
)
from omics_target_prioritization.report.dossier import render_dossier
from omics_target_prioritization.score.prioritize import score_target


def _coloc_item(gene_id: str, layer: str, score: float) -> EvidenceItem:
    return EvidenceItem(
        gene_id=gene_id,
        source=f"{layer}_coloc",  # type: ignore[arg-type]
        score=score,
        weight=1.0,
        layer=layer,  # type: ignore[arg-type]
        tissue="whole_blood",
        provenance=Provenance(dataset="sim", method="coloc.abf"),
    )


def test_mr_strength_calibration() -> None:
    """Null -> low, |z|=2 -> 0.5, strong -> near 1, and monotone in |z|."""
    assert mr_strength(0.0, 1.0) < 0.2
    assert mr_strength(MR_Z_MIDPOINT, 1.0) == pytest.approx(0.5, abs=1e-9)
    assert mr_strength(6.0, 1.0) > 0.95
    # Monotone increasing in |z|.
    zs = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0]
    strengths = [mr_strength(z, 1.0) for z in zs]
    assert strengths == sorted(strengths)
    # Sign-symmetric.
    assert mr_strength(3.0, 1.0) == pytest.approx(mr_strength(-3.0, 1.0))


def test_mr_direction_classification() -> None:
    assert mr_direction(0.3, 0.05) == "risk"  # strong positive
    assert mr_direction(-0.3, 0.05) == "protective"  # strong negative
    assert mr_direction(0.01, 0.5) == "none"  # not significant


def test_build_mr_evidence_fields() -> None:
    item = build_mr_evidence("ENSG1", beta=0.4, se=0.05, exposure_layer="pQTL")
    assert item.source == "mr_directional"
    assert 0.0 <= item.score <= 1.0
    assert item.weight == pytest.approx(0.7)
    assert item.provenance.method == "drug_target_mr"
    assert item.provenance.parameters["direction"] == "risk"
    assert item.provenance.parameters["exposure_layer"] == "pQTL"


def test_mr_does_not_inflate_layer_count() -> None:
    """MR shares the QTL instrument, so it must not count as an extra omics layer."""
    gene = Gene(gene_id="ENSG1", symbol="GENEA", chromosome="1", tss=1000)
    coloc_only = [_coloc_item("ENSG1", "eQTL", 0.9)]
    with_mr = [*coloc_only, build_mr_evidence("ENSG1", beta=0.4, se=0.05)]
    s_coloc = score_target(gene, coloc_only)
    s_mr = score_target(gene, with_mr)
    # Same single corroborating layer either way.
    assert s_coloc.n_layers == s_mr.n_layers == 1
    # But MR sharpens the aggregate score and adds the direction readout.
    assert s_mr.total >= s_coloc.total
    assert s_mr.mr_direction == "risk"
    assert s_mr.mr_pvalue is not None and s_mr.mr_pvalue < 1e-10


def test_mr_direction_surfaces_on_target() -> None:
    gene = Gene(gene_id="ENSG2", symbol="GENEB", chromosome="2", tss=5000)
    items = [
        _coloc_item("ENSG2", "eQTL", 0.85),
        build_mr_evidence("ENSG2", beta=-0.5, se=0.06),
    ]
    target = score_target(gene, items)
    assert target.mr_direction == "protective"
    # Dossier renders the directional recommendation.
    html = render_dossier(target, locus_id="chr2:5kb")
    assert "protective" in html
    assert "agonise" in html


def test_no_mr_evidence_defaults() -> None:
    gene = Gene(gene_id="ENSG3", symbol="GENEC", chromosome="3", tss=9000)
    target = score_target(gene, [_coloc_item("ENSG3", "eQTL", 0.7)])
    assert target.mr_direction == "none"
    assert target.mr_pvalue is None
