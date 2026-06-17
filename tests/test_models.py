"""Validation tests for the pydantic evidence-chain models."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from omics_target_prioritization.models import (
    EvidenceItem,
    GwasLocus,
    Provenance,
    SummaryStat,
    TargetScore,
    Variant,
)


def test_summary_stat_rejects_nonpositive_se() -> None:
    v = Variant(rsid="rs1", chromosome="1", position=10)
    with pytest.raises(ValidationError):
        SummaryStat(variant=v, beta=0.1, se=0.0, maf=0.3)


def test_gwas_locus_requires_nonempty_summary_stats() -> None:
    v = Variant(rsid="rs1", chromosome="1", position=10)
    with pytest.raises(ValidationError):
        GwasLocus(locus_id="L", lead_variant=v, summary_stats=[])


def test_evidence_score_bounded() -> None:
    prov = Provenance(dataset="d", method="m")
    with pytest.raises(ValidationError):
        EvidenceItem(gene_id="G", source="eQTL_coloc", score=1.5, provenance=prov)


def test_provenance_timestamp_defaults_non_null() -> None:
    prov = Provenance(dataset="d", method="m")
    assert prov.timestamp is not None


def test_target_score_total_bounded() -> None:
    with pytest.raises(ValidationError):
        TargetScore(gene_id="G", symbol="G", total=2.0)
