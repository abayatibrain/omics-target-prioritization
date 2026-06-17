"""Colocalization-based variant-to-gene (V2G) scoring.

For each candidate gene at a locus we build a list of provenance-stamped
:class:`~omics_target_prioritization.models.EvidenceItem`s:

1. **Colocalization evidence** — run ``coloc.abf`` between the GWAS locus and
   every QTL dataset (layer x tissue) available for the gene. The posterior
   ``PP.H4`` becomes the evidence score. This is the dominant, mechanistic V2G
   signal (ADR-0001).
2. **Distance-to-TSS** — an exponential decay in the distance between the GWAS
   lead variant and the gene's transcription start site. Weak, always-available
   prior credit; never enough on its own to win a locus.
3. **Functional annotation** (optional) — a caller-supplied per-gene score in
   ``[0, 1]`` (e.g. a VEP / coding-consequence flag) folded in as a third channel.

Every item carries a :class:`~omics_target_prioritization.models.Provenance`
record (dataset, method, UTC timestamp, parameters) so the score is auditable.
"""

from __future__ import annotations

from collections.abc import Mapping

from omics_target_prioritization.coloc import ColocPriors, colocalize
from omics_target_prioritization.models import (
    EvidenceItem,
    Gene,
    GwasLocus,
    Provenance,
    QtlAssociation,
)

#: Default weights per evidence source for the harmonic-sum aggregation.
#: Mechanistic colocalization dominates; distance is a weak prior; annotation
#: is a moderate corroborator. Defended in ADR-0002.
DEFAULT_SOURCE_WEIGHTS: dict[str, float] = {
    "eQTL_coloc": 1.0,
    "pQTL_coloc": 1.0,
    "sQTL_coloc": 0.8,
    "caQTL_coloc": 0.6,
    "distance_to_tss": 0.3,
    "functional_annotation": 0.5,
}

#: Characteristic decay length (base pairs) for the distance-to-TSS feature.
#: ~50 kb half-credit scale, consistent with typical cis-regulatory windows.
DISTANCE_DECAY_BP: float = 50_000.0


def _layer_source(layer: str) -> str:
    """Map a QTL layer label to its evidence-source key (e.g. ``eQTL_coloc``)."""
    return f"{layer}_coloc"


def distance_to_tss_score(
    gwas: GwasLocus, gene: Gene, decay_bp: float = DISTANCE_DECAY_BP
) -> float:
    r"""Exponential distance-to-TSS decay score in ``[0, 1]``.

    .. math::

        s = \exp(-|pos_{lead} - tss| / \mathrm{decay\_bp})

    Parameters
    ----------
    gwas
        GWAS locus (uses the lead variant position).
    gene
        Candidate gene (uses its TSS).
    decay_bp
        Characteristic decay length in base pairs.

    Returns
    -------
    float
        ``1.0`` when the lead variant sits on the TSS, decaying toward ``0``.
    """
    distance = abs(gwas.lead_variant.position - gene.tss)
    import math

    return math.exp(-distance / decay_bp)


def score_gene(
    gwas: GwasLocus,
    gene: Gene,
    qtls: list[QtlAssociation],
    *,
    priors: ColocPriors | None = None,
    functional_annotation: float | None = None,
    weights: Mapping[str, float] | None = None,
) -> list[EvidenceItem]:
    """Produce all evidence items linking ``gwas`` to ``gene``.

    Parameters
    ----------
    gwas
        The GWAS locus.
    gene
        The candidate gene to score.
    qtls
        All QTL datasets at the locus; only those whose ``gene`` matches are
        used for this gene's colocalization evidence.
    priors
        Coloc priors; defaults to :class:`~omics_target_prioritization.coloc.ColocPriors`.
    functional_annotation
        Optional per-gene functional-annotation score in ``[0, 1]``. When
        provided, a ``functional_annotation`` evidence item is emitted.
    weights
        Per-source weights; defaults to :data:`DEFAULT_SOURCE_WEIGHTS`.

    Returns
    -------
    list[EvidenceItem]
        One coloc item per matching QTL dataset, one distance item, and
        optionally one functional-annotation item. Each carries full provenance.
    """
    priors = priors or ColocPriors()
    weight_map = dict(weights) if weights is not None else dict(DEFAULT_SOURCE_WEIGHTS)
    items: list[EvidenceItem] = []

    # 1. Colocalization evidence, one item per matching QTL dataset.
    for qtl in qtls:
        if qtl.gene.gene_id != gene.gene_id:
            continue
        result = colocalize(gwas, qtl, priors)
        source = _layer_source(qtl.layer)
        items.append(
            EvidenceItem(
                gene_id=gene.gene_id,
                source=source,  # type: ignore[arg-type]
                score=result.pp_h4,
                weight=weight_map.get(source, 1.0),
                layer=qtl.layer,
                tissue=qtl.tissue,
                provenance=Provenance(
                    dataset=qtl.dataset,
                    method="coloc.abf",
                    parameters={
                        "p1": priors.p1,
                        "p2": priors.p2,
                        "p12": priors.p12,
                        "sd_prior": priors.sd_prior,
                        "n_snps": float(result.n_snps),
                        "PP.H3": result.pp_h3,
                        "PP.H4": result.pp_h4,
                    },
                ),
            )
        )

    # 2. Distance-to-TSS evidence (always available, weak).
    dist_score = distance_to_tss_score(gwas, gene)
    items.append(
        EvidenceItem(
            gene_id=gene.gene_id,
            source="distance_to_tss",
            score=dist_score,
            weight=weight_map.get("distance_to_tss", 0.3),
            provenance=Provenance(
                dataset="locus_geometry",
                method="distance_decay",
                parameters={
                    "decay_bp": DISTANCE_DECAY_BP,
                    "lead_pos": float(gwas.lead_variant.position),
                    "tss": float(gene.tss),
                },
            ),
        )
    )

    # 3. Optional functional-annotation evidence.
    if functional_annotation is not None:
        items.append(
            EvidenceItem(
                gene_id=gene.gene_id,
                source="functional_annotation",
                score=float(functional_annotation),
                weight=weight_map.get("functional_annotation", 0.5),
                provenance=Provenance(
                    dataset="functional_annotation",
                    method="annotation_lookup",
                    parameters={"value": float(functional_annotation)},
                ),
            )
        )

    return items


def integrate_locus(
    gwas: GwasLocus,
    genes: list[Gene],
    qtls: list[QtlAssociation],
    *,
    priors: ColocPriors | None = None,
    functional_annotations: Mapping[str, float] | None = None,
    weights: Mapping[str, float] | None = None,
) -> dict[str, list[EvidenceItem]]:
    """Score every candidate gene at a locus.

    Parameters
    ----------
    gwas
        The GWAS locus.
    genes
        Candidate genes at the locus.
    qtls
        All QTL datasets at the locus.
    priors
        Coloc priors.
    functional_annotations
        Optional mapping ``gene_id -> annotation score in [0, 1]``.
    weights
        Per-source weights.

    Returns
    -------
    dict[str, list[EvidenceItem]]
        Mapping from ``gene_id`` to its list of evidence items.
    """
    annotations = functional_annotations or {}
    out: dict[str, list[EvidenceItem]] = {}
    for gene in genes:
        out[gene.gene_id] = score_gene(
            gwas,
            gene,
            qtls,
            priors=priors,
            functional_annotation=annotations.get(gene.gene_id),
            weights=weights,
        )
    return out
