"""Deterministic simulator for a GWAS locus with multi-omics QTLs.

The simulator plants a *known answer* so the rest of the pipeline can be tested
against ground truth:

- a single **true causal gene** whose eQTL **and** pQTL signals share the GWAS
  causal variant (they will colocalize, high PP.H4);
- one or more **neighbouring genes** that either have a QTL with a *distinct*
  causal variant (no colocalization, low PP.H4) or have no QTL at all and earn
  only weak distance-to-TSS credit.

All randomness flows from a single ``numpy`` ``Generator`` seeded by
:class:`SimConfig.seed`, so a given config reproduces byte-for-byte.

The generative model is intentionally simple but faithful to what coloc.abf
consumes: a window of independent SNPs (no LD matrix needed for coloc.abf), one
designated causal SNP per signal, large ``z`` at the causal SNP and small noise
elsewhere. A *shared* causal SNP between two traits yields high PP.H4; a
*different* causal SNP yields high PP.H3 (distinct signals).
"""

from __future__ import annotations

import numpy as np
from pydantic import BaseModel

from omics_target_prioritization.models import (
    Gene,
    GwasLocus,
    QtlAssociation,
    QtlLayer,
    SummaryStat,
    Variant,
)


class SimConfig(BaseModel):
    """Configuration for the GWAS + multi-omics QTL simulation.

    Parameters
    ----------
    seed
        Master RNG seed; identical seeds give identical loci.
    n_snps
        Number of SNPs in the locus window.
    n_neighbour_genes
        Number of non-causal neighbouring genes to place around the causal one.
    chromosome
        Chromosome label for all simulated entities.
    locus_start
        1-based base-pair start of the locus window.
    snp_spacing
        Base pairs between consecutive SNPs.
    causal_z
        Approximate z-score injected at a causal SNP (controls signal strength).
    noise_se
        Standard error used for all SNPs (homoscedastic, for simplicity).
    causal_gene_tissues
        Tissues in which the causal gene has colocalizing eQTL/pQTL signals.
    """

    seed: int = 7
    n_snps: int = 200
    n_neighbour_genes: int = 3
    chromosome: str = "1"
    locus_start: int = 1_000_000
    snp_spacing: int = 500
    causal_z: float = 6.0
    noise_se: float = 0.05
    causal_gene_tissues: tuple[str, ...] = ("whole_blood",)


class SimulatedLocus(BaseModel):
    """Bundle returned by :func:`simulate_locus`.

    Parameters
    ----------
    gwas
        The simulated GWAS locus.
    genes
        All candidate genes (causal first, then neighbours).
    qtls
        All simulated QTL datasets across layers, genes and tissues.
    causal_gene_id
        Ground-truth causal gene id, for known-answer assertions.
    """

    gwas: GwasLocus
    genes: list[Gene]
    qtls: list[QtlAssociation]
    causal_gene_id: str


def _make_variants(cfg: SimConfig) -> list[Variant]:
    """Build the ordered SNP window shared by every trait in the locus."""
    variants: list[Variant] = []
    for i in range(cfg.n_snps):
        pos = cfg.locus_start + i * cfg.snp_spacing
        variants.append(
            Variant(
                rsid=f"rs{cfg.locus_start + i}",
                chromosome=cfg.chromosome,
                position=pos,
            )
        )
    return variants


def _signal_betas(
    rng: np.random.Generator,
    variants: list[Variant],
    causal_idx: int,
    cfg: SimConfig,
    *,
    has_signal: bool,
) -> list[SummaryStat]:
    """Generate per-SNP summary stats with an optional spike at ``causal_idx``.

    When ``has_signal`` is true a strong association is placed at
    ``causal_idx`` (``z ~ causal_z``) and weak null noise elsewhere. When false,
    every SNP is null noise (no QTL signal).

    Parameters
    ----------
    rng
        RNG instance.
    variants
        The shared SNP window.
    causal_idx
        Index of the causal SNP for this signal.
    cfg
        Simulation configuration.
    has_signal
        Whether to inject a causal spike.

    Returns
    -------
    list[SummaryStat]
        One :class:`SummaryStat` per SNP in ``variants``.
    """
    se = cfg.noise_se
    stats: list[SummaryStat] = []
    for i, variant in enumerate(variants):
        if has_signal and i == causal_idx:
            # Causal spike: beta = z * se, with mild jitter so runs differ.
            z = cfg.causal_z + rng.normal(0.0, 0.3)
            beta = z * se
        else:
            # Null SNP: small N(0, se) marginal noise -> |z| ~ O(1).
            beta = rng.normal(0.0, se)
        maf = float(rng.uniform(0.05, 0.5))
        stats.append(SummaryStat(variant=variant, beta=float(beta), se=se, maf=maf))
    return stats


def simulate_locus(config: SimConfig | None = None) -> SimulatedLocus:
    """Simulate one GWAS locus with planted multi-omics QTL evidence.

    The causal gene gets colocalizing eQTL **and** pQTL signals (same causal SNP
    as the GWAS). The first neighbour gets an eQTL with a *distinct* causal SNP
    (no colocalization). Remaining neighbours get no QTL at all (distance-only).

    Parameters
    ----------
    config
        Simulation configuration; defaults to :class:`SimConfig`.

    Returns
    -------
    SimulatedLocus
        The GWAS locus, candidate genes, QTL datasets, and the ground-truth
        causal gene id.
    """
    cfg = config or SimConfig()
    rng = np.random.default_rng(cfg.seed)
    variants = _make_variants(cfg)

    # The GWAS causal SNP sits near the centre of the window.
    gwas_causal_idx = cfg.n_snps // 2
    gwas_stats = _signal_betas(rng, variants, gwas_causal_idx, cfg, has_signal=True)
    lead = variants[gwas_causal_idx]
    gwas = GwasLocus(
        locus_id=f"locus_{cfg.chromosome}_{cfg.locus_start}",
        lead_variant=lead,
        summary_stats=gwas_stats,
        trait="simulated_disease",
    )

    # Causal gene's TSS sits right at the GWAS causal SNP position.
    causal_tss = variants[gwas_causal_idx].position
    causal_gene = Gene(
        gene_id="ENSG_CAUSAL",
        symbol="TRUEGENE",
        chromosome=cfg.chromosome,
        tss=causal_tss,
    )

    # Neighbour genes spread away from the causal SNP.
    genes: list[Gene] = [causal_gene]
    offset_step = max(1, cfg.n_snps // (cfg.n_neighbour_genes + 1))
    for j in range(cfg.n_neighbour_genes):
        snp_idx = min(cfg.n_snps - 1, gwas_causal_idx + (j + 1) * offset_step)
        genes.append(
            Gene(
                gene_id=f"ENSG_NBR{j + 1}",
                symbol=f"NEIGHBOR{j + 1}",
                chromosome=cfg.chromosome,
                tss=variants[snp_idx].position,
            )
        )

    qtls: list[QtlAssociation] = []

    # Causal gene: colocalizing eQTL + pQTL across the configured tissues.
    coloc_layers: tuple[QtlLayer, ...] = ("eQTL", "pQTL")
    for tissue in cfg.causal_gene_tissues:
        for layer in coloc_layers:
            stats = _signal_betas(rng, variants, gwas_causal_idx, cfg, has_signal=True)
            qtls.append(
                QtlAssociation(
                    layer=layer,
                    gene=causal_gene,
                    tissue=tissue,
                    summary_stats=stats,
                    dataset="GTEx_v8" if layer == "eQTL" else "UKB_PPP",
                )
            )

    # First neighbour: eQTL with a *distinct* causal SNP (does not colocalize).
    if cfg.n_neighbour_genes >= 1:
        distinct_idx = max(0, gwas_causal_idx - offset_step)
        stats = _signal_betas(rng, variants, distinct_idx, cfg, has_signal=True)
        qtls.append(
            QtlAssociation(
                layer="eQTL",
                gene=genes[1],
                tissue=cfg.causal_gene_tissues[0],
                summary_stats=stats,
                dataset="GTEx_v8",
            )
        )

    return SimulatedLocus(
        gwas=gwas,
        genes=genes,
        qtls=qtls,
        causal_gene_id=causal_gene.gene_id,
    )
