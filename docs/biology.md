# Biology primer

## Why "which gene?" is hard

Genome-wide association studies (GWAS) find common variants statistically
associated with a disease. Because nearby variants are correlated through
**linkage disequilibrium (LD)**, a single causal variant drags a whole haplotype
block of correlated variants to genome-wide significance. That block routinely
spans several genes. The lead (most significant) variant is usually
**non-coding** and **regulatory** — it changes *how much* of a nearby gene is
made, or *which isoform*, rather than altering a protein directly.

Picking the **nearest gene** to the lead variant is the historical default and
is wrong roughly half the time. Effector genes can sit hundreds of kilobases
away, with the causal variant looping to them through chromatin contacts.

## Molecular QTLs as the bridge

A **molecular quantitative trait locus (QTL)** is a variant associated with a
molecular phenotype:

- **eQTL** — gene **expression** level (GTEx, eQTL Catalogue).
- **pQTL** — **protein** abundance (UK Biobank Pharma Proteomics Project).
- **sQTL** — **splicing** ratios (intron excision).
- **caQTL** — **chromatin accessibility** (ATAC-seq peaks).

If the *same* causal variant that drives the disease signal also drives a gene's
eQTL, the variant plausibly causes disease *by changing that gene's expression*.
That gene is a credentialed candidate effector — and a candidate drug target.

## Colocalization, not just overlap

Two signals overlapping the same region is not enough: they might be **two
distinct causal variants** that happen to be in LD. **Colocalization** asks the
sharper question — are the disease and molecular signals best explained by **one
shared causal variant** (hypothesis H4) or **two distinct ones** (H3)? The
`coloc.abf` method returns posterior probabilities for both. We credit a gene
only when `PP.H4` is high.

## Corroboration across layers

A variant that lowers a gene's mRNA (eQTL) *and* its protein (pQTL) *and* shifts
its chromatin (caQTL), all colocalizing with disease, is far more credible than a
single-layer hit. The confidence calibration rewards independent molecular
layers agreeing, mirroring how human target-ID reviewers weigh evidence.
