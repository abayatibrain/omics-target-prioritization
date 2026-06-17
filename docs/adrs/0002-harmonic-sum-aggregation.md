# ADR-0002: Open-Targets harmonic-sum aggregation vs naive sum
Date: 2026-06-17
Status: Accepted (Armin 2026-06-17)

## Context

Each gene accrues heterogeneous evidence (coloc per layer, distance, annotation).
These must be combined into one score that ranks targets within a locus. The
combination rule shapes every downstream decision.

## Options considered

### Option A — Naive sum / mean of evidence scores
- **Pros:** trivial.
- **Cons:** a sum is unbounded and rewards *quantity* of weak evidence over
  *quality*; a gene measured in many tissues could out-sum a gene with one
  decisive colocalization. A mean dilutes a strong top signal with weak ones.
  Neither is bounded to a comparable `[0, 1]` scale across loci.

### Option B — Open-Targets weighted harmonic sum
Sort weighted scores descending, combine as `Σ s_(i) / i²`, normalize by the
theoretical maximum `Σ 1/i²` so the result is in `[0, 1]`.
- **Pros:** rewards a strong *top* signal while still giving diminishing credit
  to corroborating evidence; bounded and comparable across loci; monotonic in
  evidence. This is the established Open Targets association-scoring rule.
- **Cons:** the `1/i²` weighting is a modelling choice, not derived from first
  principles; per-source weights must be set deliberately.

## Decision

**Use the Open-Targets weighted harmonic sum (Option B).** Within a source
(e.g. eQTL across tissues) we take the *strongest* item so multi-tissue
measurement cannot inflate the score. Per-source weights privilege mechanistic
colocalization over the distance prior (coloc 1.0, distance 0.3, annotation 0.5).

## Consequences

- Scores are bounded `[0, 1]` and monotonic (asserted in tests), so they are
  interpretable as "more / stronger evidence -> higher score."
- The weighting and per-source-best rules are documented in `aggregate.py` and
  unit-tested.
- Changing weights is a modelling decision that should be justified, ideally by
  calibration against a labelled V2G benchmark (future work).

## References

- Ghoussaini M. *et al.* (2021) Open Targets Genetics. *Nucleic Acids Res*.
