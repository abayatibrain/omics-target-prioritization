# Methods

## Colocalization (coloc.abf)

We implement the approximate-Bayes-factor colocalization of
**Giambartolomei et al. (2014)**. For each SNP, the **Wakefield approximate
Bayes factor (ABF)** compares "this SNP is causal" to "null" using only the
marginal effect `beta` and its standard error `se`:

$$
\mathrm{ABF} = \sqrt{\frac{V}{V+W}}\,\exp\!\left(\frac{z^2}{2}\cdot\frac{W}{V+W}\right),
\qquad V = se^2,\; z = \beta/se,\; W = \mathrm{sd\_prior}^2 .
$$

Per-SNP ABFs for the two traits are combined, in log-space with the
log-sum-exp trick, into posterior probabilities of the five hypotheses
`H0..H4` under priors `p1, p2, p12`. **`PP.H4`** — one shared causal variant —
is the colocalization posterior used for variant-to-gene mapping. The default
priors (`p1 = p2 = 1e-4`, `p12 = 1e-5`) follow the coloc defaults.

## Variant-to-gene (V2G)

For each candidate gene we collect:

1. **Colocalization evidence** — `PP.H4` from `coloc.abf` against every QTL
   layer/tissue available for the gene. Modelled on **Open Targets Genetics**
   V2G, which uses molecular-QTL colocalization as a primary line of evidence.
   QTL sources mirror **GTEx** and the **eQTL Catalogue** (eQTL/sQTL) and the
   **UK Biobank Pharma Proteomics Project (UKB-PPP)** (pQTL).
2. **Distance-to-TSS** — an exponential decay `exp(-d / 50kb)` in the distance
   between the GWAS lead variant and the gene's TSS. A weak, always-available
   prior.
3. **Functional annotation** (optional) — a caller-supplied `[0,1]` score
   (e.g. coding-consequence / VEP severity) folded in as a third channel.

Each becomes a `Provenance`-stamped `EvidenceItem`.

## Aggregation (Open Targets harmonic sum)

Evidence is combined with the **Open Targets harmonic-sum** scoring. Per-source
scores (weighted, then the strongest per source) are sorted descending and
combined as

$$
\mathrm{HS} = \sum_{i=1}^{n} \frac{s_{(i)}}{i^2},
$$

then normalized by the theoretical maximum $\sum_i 1/i^2$ so the result lies in
`[0, 1]`. The harmonic sum rewards a strong top signal while giving diminishing
credit to additional corroborating evidence — and it is **monotonic**: more or
stronger evidence never lowers the score.

## Calibrated confidence

A deterministic, documented rule maps the number of corroborating QTL layers,
the maximum `PP.H4`, and the aggregate score to a **high / medium / low** label.
*High* confidence requires genuine multi-omics corroboration (≥2 layers with
`PP.H4 ≥ 0.5` and total ≥ 0.5), not a single strong layer or a strong distance
prior. See ADR-0003.

## References

- Giambartolomei C. *et al.* (2014) *PLoS Genet* 10(5):e1004383.
- Ghoussaini M. *et al.* (2021) Open Targets Genetics. *Nucleic Acids Res*.
- GTEx Consortium (2020) *Science* 369:1318–1330.
- Kerimov N. *et al.* (2021) eQTL Catalogue. *Nat Genet*.
- Sun B.B. *et al.* (2023) UK Biobank PPP plasma pQTLs. *Nature* 622:329–338.
