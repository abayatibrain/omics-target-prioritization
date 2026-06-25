"""Drug-target Mendelian-randomization evidence channel.

Colocalization answers "does the disease act *through* this gene?" but it is
direction-agnostic: a high PP.H4 does not tell you whether you want an agonist
or an antagonist. Drug-target MR closes that gap. Using the gene's QTL as a
genetic instrument for its own expression / protein level, a two-sample MR
estimates the *causal* effect of perturbing the gene product on disease risk,
with a sign: does raising the product raise risk (so you would inhibit) or
lower it (so you would agonise)?

This module turns an MR estimate into an :class:`~omics_target_prioritization.models.EvidenceItem`
on the ``"mr_directional"`` channel. Two deliberate design choices:

* **Calibrated strength, not 1 - p.** The evidence score is a logistic function
  of the absolute MR Z-statistic, centred so that the conventional
  significance threshold (``|z| = 2``) maps to ``0.5`` rather than to ``0.95``.
  A barely-significant MR should not look like overwhelming evidence.
* **MR is not an independent omics layer.** Because the instrument is the same
  QTL that drives the colocalization, MR derived from it corroborates *the same*
  signal. It therefore contributes to the aggregate score (it sharpens and
  signs the mechanism) but is deliberately excluded from the multi-layer
  ``confidence`` count, which is reserved for genuinely distinct molecular
  layers.

References
----------
Schmidt, A. F. et al. (2020). Genetic drug target validation using Mendelian
randomisation. *Nature Communications* 11:3255.
"""

from __future__ import annotations

import math

from scipy import stats

from omics_target_prioritization.models import (
    EvidenceItem,
    MrDirection,
    Provenance,
    QtlLayer,
)

__all__ = ["MR_Z_MIDPOINT", "MR_Z_SLOPE", "build_mr_evidence", "mr_direction", "mr_strength"]

#: Absolute Z at which MR evidence strength crosses 0.5 (the conventional
#: two-sided significance threshold, ``p ≈ 0.05``).
MR_Z_MIDPOINT: float = 2.0

#: Logistic slope in the |z| -> strength map. Larger is sharper.
MR_Z_SLOPE: float = 1.0


def mr_strength(beta: float, se: float) -> float:
    """Map an MR estimate to a calibrated evidence strength in ``[0, 1]``.

    The strength is ``logistic(slope * (|z| - midpoint))`` with ``z = beta/se``,
    so a null effect (``z = 0``) is near zero, a borderline-significant effect
    (``|z| = 2``) is ``0.5``, and a strong effect (``|z| >= 5``) approaches 1.

    Parameters
    ----------
    beta
        MR causal-effect estimate (gene product -> disease).
    se
        Standard error of ``beta`` (must be positive).

    Returns
    -------
    float
        Calibrated evidence strength in ``[0, 1]``.
    """
    if se <= 0:
        raise ValueError("se must be positive")
    z = abs(beta / se)
    return 1.0 / (1.0 + math.exp(-MR_Z_SLOPE * (z - MR_Z_MIDPOINT)))


def mr_direction(beta: float, se: float, *, alpha: float = 0.05) -> MrDirection:
    """Classify the direction of an MR effect, or ``"none"`` if not significant.

    Parameters
    ----------
    beta
        MR causal estimate; positive means raising the gene product raises risk.
    se
        Standard error of ``beta``.
    alpha
        Two-sided significance threshold below which the direction is reported.

    Returns
    -------
    MrDirection
        ``"risk"`` (beta > 0), ``"protective"`` (beta < 0), or ``"none"``.
    """
    z = beta / se
    pvalue = 2.0 * stats.norm.sf(abs(z))
    if pvalue > alpha:
        return "none"
    return "risk" if beta > 0 else "protective"


def build_mr_evidence(
    gene_id: str,
    beta: float,
    se: float,
    *,
    exposure_layer: QtlLayer = "eQTL",
    dataset: str = "simulated",
    weight: float = 0.7,
) -> EvidenceItem:
    """Construct an ``mr_directional`` evidence item from an MR estimate.

    Parameters
    ----------
    gene_id
        The gene whose product was the MR exposure.
    beta, se
        MR causal estimate and its standard error.
    exposure_layer
        Which molecular layer supplied the instrument (recorded for audit).
    dataset
        Source dataset label for provenance.
    weight
        Aggregation weight. Defaults below 1.0 so MR corroborates and signs the
        mechanism without out-voting the primary colocalization evidence.

    Returns
    -------
    EvidenceItem
        A provenance-carrying evidence item on the ``"mr_directional"`` channel.
        The signed ``beta``, ``se``, ``pvalue`` and ``direction`` are stored in
        the provenance parameters so the dossier can state which way to drug it.
    """
    z = beta / se
    pvalue = float(2.0 * stats.norm.sf(abs(z)))
    direction = mr_direction(beta, se)
    score = mr_strength(beta, se)
    provenance = Provenance(
        dataset=dataset,
        method="drug_target_mr",
        parameters={
            "beta": float(beta),
            "se": float(se),
            "z": float(z),
            "pvalue": pvalue,
            "direction": direction,
            "exposure_layer": exposure_layer,
        },
    )
    return EvidenceItem(
        gene_id=gene_id,
        source="mr_directional",
        score=score,
        weight=weight,
        layer=exposure_layer,
        provenance=provenance,
    )
