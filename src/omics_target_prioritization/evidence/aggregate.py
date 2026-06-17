"""Open-Targets-style weighted harmonic-sum evidence aggregation.

Open Targets combines heterogeneous evidence into a single target score with a
**harmonic sum**: scores are sorted in descending order and combined as

.. math::

    \\mathrm{HS} = \\sum_{i=1}^{n} \\frac{s_{(i)}}{(i)^2},

where :math:`s_{(1)} \\ge s_{(2)} \\ge \\dots` are the sorted scores. The harmonic
sum rewards having a *strong top* signal while still giving partial credit for
additional corroborating evidence.

To make the score bounded **and monotonic under adding evidence**, we normalize
by the supremum of the harmonic sum over unboundedly many unit scores — the
Basel-series limit :math:`\\sum_{i=1}^{\\infty} 1/i^2 = \\pi^2/6 \\approx 1.6449`.
Normalizing by this *fixed* constant (rather than by the per-call partial sum)
is what Open Targets does, and it guarantees that appending a non-negative score
can only increase the result (ADR-0002).

We extend the vanilla harmonic sum with **per-source weights**: each evidence
score is multiplied by its source weight before sorting, which lets mechanistic
colocalization evidence outrank a weak distance prior even when both are present.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from omics_target_prioritization.models import EvidenceItem

#: Supremum of the harmonic sum over infinitely many unit scores (Basel series,
#: ``Σ 1/i² = π²/6``). Normalizing by this fixed constant keeps the score in
#: ``[0, 1]`` *and* monotonic when evidence is added (ADR-0002).
HARMONIC_SUM_MAX: float = math.pi**2 / 6.0


def harmonic_sum(scores: Sequence[float]) -> float:
    """Normalized harmonic sum of a set of scores, in ``[0, 1]``.

    Scores are sorted descending and combined as ``Σ s_(i) / i²`` (1-based
    ``i``), then divided by the fixed Basel-series limit
    :data:`HARMONIC_SUM_MAX` (``π²/6``) so the result is bounded in ``[0, 1]``
    and monotonic under adding evidence.

    Parameters
    ----------
    scores
        Evidence scores, each assumed to be in ``[0, 1]``.

    Returns
    -------
    float
        Normalized harmonic sum in ``[0, 1]``; ``0.0`` for an empty input.

    Notes
    -----
    The function is **monotonic**: adding a non-negative score, or increasing
    any existing score, never decreases the result. This is asserted in the
    tests and is the property that makes the aggregate interpretable as
    "more / stronger evidence -> higher score".
    """
    clean = [max(0.0, float(s)) for s in scores]
    if not clean:
        return 0.0
    ordered = sorted(clean, reverse=True)
    numerator = sum(s / (i**2) for i, s in enumerate(ordered, start=1))
    return min(1.0, numerator / HARMONIC_SUM_MAX)


def aggregate_scores(
    items: list[EvidenceItem],
) -> tuple[float, dict[str, float]]:
    """Aggregate evidence items into an overall score and a per-source breakdown.

    Each item contributes ``score * weight`` (clamped to ``[0, 1]``) to the
    harmonic sum. When several items share a source (e.g. eQTL coloc in multiple
    tissues), the *strongest* one is taken so a source cannot inflate the score
    merely by being measured in many tissues.

    Parameters
    ----------
    items
        Evidence items for a single gene.

    Returns
    -------
    total : float
        Overall target score in ``[0, 1]``.
    breakdown : dict[str, float]
        Mapping ``source -> best raw (unweighted) score`` for that source,
        for transparency in the dossier.
    """
    if not items:
        return 0.0, {}

    best_per_source: dict[str, float] = {}
    best_weighted_per_source: dict[str, float] = {}
    for item in items:
        weighted = min(1.0, max(0.0, item.score) * item.weight)
        if item.score > best_per_source.get(item.source, -1.0):
            best_per_source[item.source] = item.score
        if weighted > best_weighted_per_source.get(item.source, -1.0):
            best_weighted_per_source[item.source] = weighted

    total = harmonic_sum(list(best_weighted_per_source.values()))
    return total, best_per_source
