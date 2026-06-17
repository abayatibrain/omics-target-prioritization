# ADR-0001: Colocalization-based variant-to-gene vs nearest-gene
Date: 2026-06-17
Status: Accepted (Armin 2026-06-17)

## Context

A GWAS locus implicates a region containing several genes. The pipeline must
decide *which* gene the causal variant acts through — the variant-to-gene (V2G)
problem. The two broad strategies are positional (nearest gene / distance) and
mechanistic (molecular-QTL colocalization).

## Options considered

### Option A — Nearest-gene / distance-to-TSS only
Assign the locus to the closest gene, or weight genes by distance.
- **Pros:** trivial, always computable, no QTL data required.
- **Cons:** wrong roughly half the time. Causal genes are frequently *not* the
  nearest gene; non-coding regulatory variants act at distance through chromatin
  looping. Distance carries no mechanistic evidence.

### Option B — Colocalization-based V2G
Test whether the GWAS signal and a gene's molecular QTL (eQTL/pQTL/sQTL/caQTL)
share a single causal variant; credit genes by `PP.H4`.
- **Pros:** mechanistic — directly tests "the disease variant changes *this*
  gene's expression/protein/splicing/accessibility." Multi-omics layers give
  independent corroboration. This is the approach used by Open Targets Genetics
  and standard statistical-genetics target-ID pipelines.
- **Cons:** needs QTL summary statistics; coloc.abf assumes a single causal
  variant per locus (mitigated in future work via SuSiE fine-mapping).

## Decision

**Primary V2G evidence is colocalization (`coloc.abf`, Option B).** Distance-to-TSS
is retained as a *weak prior* feature, never sufficient on its own to win a
locus. The harmonic-sum weights (ADR-0002) encode this: coloc weight 1.0,
distance weight 0.3.

## Consequences

- The pipeline requires QTL datasets per gene; the simulator and (future)
  real-data loaders must harmonize SNP order before coloc.
- A gene that is merely *close* to the lead variant cannot reach high confidence;
  it needs colocalizing molecular evidence.
- The single-causal-variant assumption of coloc.abf is a known limitation,
  recorded for the SuSiE-coloc extension in a future ADR.

## References

- Giambartolomei C. *et al.* (2014) *PLoS Genet* 10(5):e1004383.
- Ghoussaini M. *et al.* (2021) Open Targets Genetics. *Nucleic Acids Res*.
