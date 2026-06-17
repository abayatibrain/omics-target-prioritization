"""Calibrated confidence labelling for a target's evidence (ADR-0003).

A numeric score alone is hard to act on; reviewers want a calibrated
high/medium/low label. The label is a deterministic function of three
auditable quantities:

1. **Number of independent corroborating omics layers** (``n_layers``): how many
   *distinct* QTL layers (eQTL / pQTL / sQTL / caQTL) colocalize with the GWAS
   signal at a meaningful PP.H4. Independent molecular layers agreeing is the
   strongest form of corroboration.
2. **Maximum PP.H4** (``max_h4``): the strongest single colocalization posterior.
3. **Aggregate score** (``total``): the Open-Targets harmonic sum.

The calibration rule (documented in :data:`CONFIDENCE_RULE_DOC`) is intentionally
conservative: *high* confidence requires genuine multi-omics corroboration, not
a single strong layer or a strong distance prior.
"""

from __future__ import annotations

from collections.abc import Iterable

from omics_target_prioritization.models import ConfidenceLabel, EvidenceItem

#: PP.H4 above which a colocalization counts as a "corroborating layer".
H4_LAYER_THRESHOLD: float = 0.5

#: PP.H4 required for a single layer to reach "medium" on its own.
H4_MEDIUM_THRESHOLD: float = 0.8

#: Aggregate score required (with corroboration) to reach "high".
SCORE_HIGH_THRESHOLD: float = 0.5

#: Human-readable statement of the calibration rule, surfaced in the dossier.
CONFIDENCE_RULE_DOC: str = (
    "high   := (>=2 independent QTL layers with PP.H4 >= 0.5) AND total >= 0.5; "
    "medium := (>=1 QTL layer with PP.H4 >= 0.8) OR (>=2 layers with PP.H4 >= 0.5); "
    "low    := otherwise."
)

#: The colocalization evidence sources that count as "omics layers".
_COLOC_SOURCES = frozenset({"eQTL_coloc", "pQTL_coloc", "sQTL_coloc", "caQTL_coloc"})


def count_corroborating_layers(items: Iterable[EvidenceItem]) -> int:
    """Count distinct QTL layers colocalizing above :data:`H4_LAYER_THRESHOLD`.

    Parameters
    ----------
    items
        Evidence items for a single gene.

    Returns
    -------
    int
        Number of distinct colocalization layers with PP.H4 above threshold.
    """
    layers: set[str] = set()
    for item in items:
        if item.source in _COLOC_SOURCES and item.score >= H4_LAYER_THRESHOLD:
            layers.add(item.source)
    return len(layers)


def max_h4(items: Iterable[EvidenceItem]) -> float:
    """Return the maximum colocalization PP.H4 across evidence items.

    Parameters
    ----------
    items
        Evidence items for a single gene.

    Returns
    -------
    float
        Largest PP.H4 among colocalization items; ``0.0`` if there are none.
    """
    h4s = [item.score for item in items if item.source in _COLOC_SOURCES]
    return max(h4s, default=0.0)


def assign_confidence(
    items: Iterable[EvidenceItem],
    total: float,
) -> tuple[ConfidenceLabel, int, float]:
    """Assign a calibrated confidence label from evidence and aggregate score.

    See :data:`CONFIDENCE_RULE_DOC` for the exact rule.

    Parameters
    ----------
    items
        Evidence items for a single gene.
    total
        The gene's aggregate harmonic-sum score in ``[0, 1]``.

    Returns
    -------
    label : ConfidenceLabel
        ``"high"`` / ``"medium"`` / ``"low"``.
    n_layers : int
        Number of corroborating layers (also returned for the TargetScore).
    max_pp_h4 : float
        Maximum PP.H4 (also returned for the TargetScore).
    """
    item_list = list(items)
    n_layers = count_corroborating_layers(item_list)
    best_h4 = max_h4(item_list)

    if n_layers >= 2 and total >= SCORE_HIGH_THRESHOLD:
        label: ConfidenceLabel = "high"
    elif best_h4 >= H4_MEDIUM_THRESHOLD or n_layers >= 2:
        label = "medium"
    else:
        label = "low"

    return label, n_layers, best_h4
