"""Tests for the Giambartolomei-2014 coloc.abf implementation."""

from __future__ import annotations

import numpy as np

from omics_target_prioritization.coloc import (
    ColocPriors,
    coloc_abf,
    colocalize,
)
from omics_target_prioritization.models import (
    Gene,
    GwasLocus,
    QtlAssociation,
    SummaryStat,
    Variant,
)


def _stats(betas: list[float], se: float = 0.05) -> list[SummaryStat]:
    out = []
    for i, b in enumerate(betas):
        out.append(
            SummaryStat(
                variant=Variant(rsid=f"rs{i}", chromosome="1", position=1000 + i),
                beta=b,
                se=se,
                maf=0.3,
            )
        )
    return out


def test_posteriors_sum_to_one() -> None:
    rng = np.random.default_rng(0)
    log_abf1 = rng.normal(0, 1, size=50)
    log_abf2 = rng.normal(0, 1, size=50)
    res = coloc_abf(log_abf1, log_abf2)
    total = res.pp_h0 + res.pp_h1 + res.pp_h2 + res.pp_h3 + res.pp_h4
    assert abs(total - 1.0) < 1e-9
    assert res.n_snps == 50


def test_shared_causal_gives_high_h4() -> None:
    # Same SNP (index 25) carries the signal in both traits -> colocalization.
    n = 60
    g_betas = [0.02] * n
    q_betas = [0.02] * n
    g_betas[25] = 0.4  # z = 8
    q_betas[25] = 0.4
    gwas = GwasLocus(
        locus_id="L",
        lead_variant=Variant(rsid="rs25", chromosome="1", position=1025),
        summary_stats=_stats(g_betas),
    )
    gene = Gene(gene_id="G", symbol="G", chromosome="1", tss=1025)
    qtl = QtlAssociation(layer="eQTL", gene=gene, tissue="blood", summary_stats=_stats(q_betas))
    res = colocalize(gwas, qtl)
    assert res.pp_h4 > 0.9
    assert res.pp_h4 > res.pp_h3


def test_distinct_causal_gives_low_h4() -> None:
    # Different causal SNPs -> H3 (distinct) dominates, H4 is low.
    n = 60
    g_betas = [0.02] * n
    q_betas = [0.02] * n
    g_betas[10] = 0.5  # GWAS causal at SNP 10
    q_betas[45] = 0.5  # QTL causal at SNP 45
    gwas = GwasLocus(
        locus_id="L",
        lead_variant=Variant(rsid="rs10", chromosome="1", position=1010),
        summary_stats=_stats(g_betas),
    )
    gene = Gene(gene_id="G", symbol="G", chromosome="1", tss=1045)
    qtl = QtlAssociation(layer="eQTL", gene=gene, tissue="blood", summary_stats=_stats(q_betas))
    res = colocalize(gwas, qtl)
    assert res.pp_h4 < 0.1
    assert res.pp_h3 > res.pp_h4


def test_priors_are_configurable() -> None:
    priors = ColocPriors(p1=1e-4, p2=1e-4, p12=1e-6, sd_prior=0.2)
    log_abf1 = np.array([3.0, 0.1, 0.1])
    log_abf2 = np.array([3.0, 0.1, 0.1])
    res = coloc_abf(log_abf1, log_abf2, priors)
    assert 0.0 <= res.pp_h4 <= 1.0
