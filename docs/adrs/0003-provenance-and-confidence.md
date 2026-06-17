# ADR-0003: Explicit provenance + calibrated confidence design
Date: 2026-06-17
Status: Accepted (Armin 2026-06-17)

## Context

A prioritization score that cannot be audited is not actionable in a drug-
discovery setting. Reviewers need to know *why* a target ranked where it did and
*how much to trust it*. Two design decisions follow: how to record provenance,
and how to express confidence.

## Decision 1 — Mandatory provenance on every evidence item

Every `EvidenceItem` carries a non-optional `Provenance` record with:

- **dataset** — the source the evidence came from (e.g. `GTEx_v8`, `UKB_PPP`);
- **method** — the algorithm that produced the score (e.g. `coloc.abf`,
  `distance_decay`);
- **timestamp** — a UTC wall-clock anchor, defaulted to *now* so it is never
  null;
- **parameters** — the parameters used (e.g. coloc priors, PP.H3/H4).

Provenance is validated at construction; an evidence item cannot exist without
it. Tests assert non-null `dataset`, `method`, and `timestamp` on every item.

## Decision 2 — Calibrated, documented confidence rule

Rather than emit only a continuous score, every `TargetScore` carries a
**calibrated high/medium/low label** computed from three auditable quantities:
the number of independent corroborating QTL layers, the maximum `PP.H4`, and the
aggregate harmonic-sum score. The rule is:

```
high   := (>=2 independent QTL layers with PP.H4 >= 0.5) AND total >= 0.5
medium := (>=1 QTL layer with PP.H4 >= 0.8) OR (>=2 layers with PP.H4 >= 0.5)
low    := otherwise
```

The rule is conservative by design: *high* confidence requires genuine
multi-omics corroboration, not a single strong layer or a strong distance prior.
The thresholds are named constants in `confidence.py` and the rule text is
surfaced in the HTML dossier so a reviewer reads the calibration alongside the
label.

## Options considered (for confidence)

- **Continuous score only** — rejected; reviewers want a calibrated verdict, and
  a bare number invites over-interpretation of small differences.
- **Score-threshold-only label** — rejected; a strong distance prior could lift
  a non-colocalizing gene over a threshold without mechanistic support.
- **Layer-corroboration + PP.H4 + score (chosen)** — ties the label to the
  *kind* of evidence, not just its magnitude.

## Consequences

- The score is fully auditable; the dossier prints a provenance table and the
  calibration rule.
- Confidence rises with corroborating layers (asserted in tests).
- Recalibrating thresholds against a labelled benchmark is future work and would
  warrant a superseding ADR.
