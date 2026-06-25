"""Pydantic data models for the GWAS-to-target evidence chain.

These models are the typed contract that flows through the whole pipeline:

    GwasLocus + [QtlAssociation ...]  ->  v2g  ->  [EvidenceItem ...]
                                                      |
                                            aggregate + confidence
                                                      v
                                                 TargetScore

Every model is immutable-by-convention (validated at construction) and every
piece of evidence carries explicit :class:`Provenance` so a reviewer can audit
*why* a target ranked where it did (see ADR-0003).
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone, tzinfo
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator

#: UTC tzinfo built from a zero offset. Equivalent to ``datetime.timezone.utc``
#: / the 3.11+ ``datetime.UTC`` alias, but constructed explicitly so the module
#: imports cleanly on both 3.10 and 3.11 (production / CI target 3.11 per
#: ``pyproject.toml``). Call sites read ``datetime.now(UTC)``.
UTC: tzinfo = timezone(timedelta(0))

#: The multi-omics QTL layers this pipeline understands.
QtlLayer = Literal["eQTL", "pQTL", "sQTL", "caQTL"]

#: The kinds of evidence a gene can accrue toward its target score.
EvidenceSource = Literal[
    "eQTL_coloc",
    "pQTL_coloc",
    "sQTL_coloc",
    "caQTL_coloc",
    "distance_to_tss",
    "functional_annotation",
    "mr_directional",
]

#: Direction of a drug-target MR effect: does *raising* the gene product raise
#: ("risk") or lower ("protective") disease risk? ``"none"`` when there is no
#: significant directional MR evidence.
MrDirection = Literal["risk", "protective", "none"]

#: Calibrated confidence labels (see ADR-0003 / ``evidence/confidence.py``).
ConfidenceLabel = Literal["high", "medium", "low"]


class Variant(BaseModel):
    """A single biallelic genetic variant.

    Parameters
    ----------
    rsid
        Variant identifier (e.g. ``"rs12345"`` or a ``chr:pos:ref:alt`` ID).
    chromosome
        Chromosome label, without the ``chr`` prefix (e.g. ``"1"``, ``"X"``).
    position
        1-based base-pair position on ``chromosome``.
    effect_allele
        The allele the GWAS / QTL effect size is reported with respect to.
    other_allele
        The reference (non-effect) allele.
    """

    model_config = ConfigDict(frozen=True)

    rsid: str
    chromosome: str
    position: int = Field(ge=1)
    effect_allele: str = "A"
    other_allele: str = "G"


class Gene(BaseModel):
    """A protein-coding gene candidate at a locus.

    Parameters
    ----------
    gene_id
        Stable identifier (e.g. an Ensembl gene ID).
    symbol
        HGNC gene symbol.
    chromosome
        Chromosome the gene is on.
    tss
        Transcription start site position (1-based). Used by the
        distance-to-TSS feature in :mod:`omics_target_prioritization.integrate.v2g`.
    """

    model_config = ConfigDict(frozen=True)

    gene_id: str
    symbol: str
    chromosome: str
    tss: int = Field(ge=1)


class SummaryStat(BaseModel):
    """Per-SNP association summary statistics for a single trait.

    These are the quantities :mod:`omics_target_prioritization.coloc` needs:
    the marginal effect estimate and its standard error. ``beta`` may be on a
    log-OR scale (binary GWAS) or a linear scale (molecular QTL); coloc.abf is
    invariant to that choice because it works from ``beta`` and ``se`` alone.

    Parameters
    ----------
    variant
        The variant these statistics describe.
    beta
        Marginal effect estimate.
    se
        Standard error of ``beta`` (strictly positive).
    maf
        Minor-allele frequency in ``(0, 0.5]``. Used for the variance prior in
        coloc when an effect-size prior on the standardized scale is requested.
    """

    model_config = ConfigDict(frozen=True)

    variant: Variant
    beta: float
    se: float = Field(gt=0.0)
    maf: float = Field(gt=0.0, le=0.5)


class GwasLocus(BaseModel):
    """A GWAS association locus: a lead variant and its per-SNP statistics.

    Parameters
    ----------
    locus_id
        Human-readable identifier for the locus.
    lead_variant
        The most significant variant at the locus.
    summary_stats
        Per-SNP GWAS summary statistics across the locus window. Order defines
        the SNP index shared with QTL datasets during colocalization.
    sample_size
        GWAS effective sample size (used only for documentation / provenance).
    trait
        Free-text trait / disease label.
    """

    model_config = ConfigDict(frozen=True)

    locus_id: str
    lead_variant: Variant
    summary_stats: list[SummaryStat]
    sample_size: int = Field(default=100_000, ge=1)
    trait: str = "disease"

    @field_validator("summary_stats")
    @classmethod
    def _non_empty(cls, value: list[SummaryStat]) -> list[SummaryStat]:
        if not value:
            raise ValueError("GwasLocus.summary_stats must be non-empty")
        return value


class QtlAssociation(BaseModel):
    """A molecular-QTL dataset for one gene in one tissue.

    Parameters
    ----------
    layer
        Which omics layer this QTL comes from (eQTL / pQTL / sQTL / caQTL).
    gene
        The gene whose molecular trait is being mapped.
    tissue
        Tissue / cell-type / context label (e.g. ``"whole_blood"``).
    summary_stats
        Per-SNP QTL summary statistics, indexed compatibly with the
        :class:`GwasLocus` they will be colocalized against.
    dataset
        Source dataset label, recorded into evidence provenance
        (e.g. ``"GTEx_v8"``, ``"UKB_PPP"``).
    sample_size
        QTL study sample size (documentation / provenance only).
    """

    model_config = ConfigDict(frozen=True)

    layer: QtlLayer
    gene: Gene
    tissue: str
    summary_stats: list[SummaryStat]
    dataset: str = "simulated"
    sample_size: int = Field(default=1_000, ge=1)

    @field_validator("summary_stats")
    @classmethod
    def _non_empty(cls, value: list[SummaryStat]) -> list[SummaryStat]:
        if not value:
            raise ValueError("QtlAssociation.summary_stats must be non-empty")
        return value


class Provenance(BaseModel):
    """Audit trail attached to every :class:`EvidenceItem` (ADR-0003).

    Parameters
    ----------
    dataset
        The data source the evidence was derived from.
    method
        The method / algorithm that produced the score
        (e.g. ``"coloc.abf"``, ``"distance_decay"``).
    timestamp
        UTC timestamp at which the evidence was generated. Defaults to *now*;
        never null, so audits always have a wall-clock anchor.
    parameters
        Free-form record of the parameters used (e.g. coloc priors).
    """

    model_config = ConfigDict(frozen=True)

    dataset: str
    method: str
    timestamp: datetime = Field(default_factory=lambda: datetime.now(UTC))
    parameters: dict[str, float | str] = Field(default_factory=dict)


class EvidenceItem(BaseModel):
    """One scored piece of evidence linking a variant signal to a gene.

    Parameters
    ----------
    gene_id
        The gene this evidence supports.
    source
        Which evidence channel produced it (coloc layer, distance, annotation).
    score
        Evidence strength in ``[0, 1]`` (e.g. a PP.H4 posterior, or a
        normalized distance-decay weight).
    weight
        Relative importance of this source in the harmonic-sum aggregation
        (ADR-0002). Strictly positive.
    layer
        QTL layer, when the evidence is a colocalization result.
    tissue
        Tissue / context, when applicable.
    provenance
        Mandatory audit trail (ADR-0003).
    """

    model_config = ConfigDict(frozen=True)

    gene_id: str
    source: EvidenceSource
    score: float = Field(ge=0.0, le=1.0)
    weight: float = Field(gt=0.0, default=1.0)
    layer: QtlLayer | None = None
    tissue: str | None = None
    provenance: Provenance


class TargetScore(BaseModel):
    """The aggregated, calibrated prioritization result for one gene.

    Parameters
    ----------
    gene_id
        The scored gene.
    symbol
        HGNC symbol, copied through for human-readable reports.
    total
        Overall target score in ``[0, 1]`` from the harmonic-sum aggregation.
    breakdown
        Per-source contribution map (source label -> raw evidence score).
    confidence
        Calibrated high/medium/low label (ADR-0003).
    n_layers
        Number of independent corroborating omics layers (distinct QTL layers
        with non-trivial colocalization). Drives ``confidence``.
    max_h4
        Maximum PP.H4 across all colocalization evidence for this gene.
    mr_direction
        Direction of the strongest drug-target MR effect for this gene: whether
        *raising* the gene product is predicted to raise (``"risk"``) or lower
        (``"protective"``) disease risk. ``"none"`` if no MR evidence is present.
        This is the "which way to drug it" readout; it does not change ``total``
        on its own, it annotates it.
    mr_pvalue
        P-value of the strongest MR effect, or ``None`` when absent.
    evidence
        The full list of :class:`EvidenceItem`s that produced ``total``.
    """

    model_config = ConfigDict(frozen=True)

    gene_id: str
    symbol: str
    total: float = Field(ge=0.0, le=1.0)
    breakdown: dict[str, float] = Field(default_factory=dict)
    confidence: ConfidenceLabel = "low"
    n_layers: int = Field(ge=0, default=0)
    max_h4: float = Field(ge=0.0, le=1.0, default=0.0)
    mr_direction: MrDirection = "none"
    mr_pvalue: float | None = None
    evidence: list[EvidenceItem] = Field(default_factory=list)
