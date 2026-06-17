# omics-target-prioritization

Integrate GWAS signals with multi-omics QTLs (eQTL / pQTL / sQTL / caQTL) via
**colocalization-based variant-to-gene mapping**, then aggregate the evidence
into a **calibrated, provenance-tracked target-prioritization score** and a
per-target dossier.

## The question

A GWAS locus implicates a *region*, not a *gene*. Drug discovery needs the gene.
This pipeline asks: given a disease GWAS locus, which gene is the causal variant
most plausibly acting through, and how confident should we be — backed by what
auditable evidence?

## How it works

1. **Colocalization** — `coloc.abf` (Giambartolomei 2014) tests whether the GWAS
   signal and each QTL share a single causal variant. `PP.H4` is the V2G signal.
2. **Variant-to-gene** — coloc per layer/tissue + distance-to-TSS decay +
   optional functional annotation, each a provenance-stamped evidence item.
3. **Aggregation** — Open-Targets harmonic sum into a score in `[0, 1]`.
4. **Confidence** — calibrated high/medium/low from corroborating layers + PP.H4.
5. **Dossier** — an HTML one-pager per target with evidence and provenance.

See the [biology primer](biology.md), the [architecture](architecture.md), the
[methods](methods.md), and the decision log.

## Quick start

```bash
pip install -e .
omics-target-prioritization run-all --out results/
```
