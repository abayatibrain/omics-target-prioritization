# Status — 2026-06-17
Repo: omics-target-prioritization
Phase: **v0.1.0 — full simulation-driven pipeline landed.**

## What ships in v0.1.0
- **Evidence-chain models** (`models.py`): `Variant`, `Gene`, `QtlAssociation`,
  `GwasLocus`, `Provenance`, `EvidenceItem`, `TargetScore` — pydantic v2,
  fully validated.
- **Simulator** (`simulate.py`): deterministic `SimConfig`-driven generation of
  a GWAS locus with a planted true causal gene that has colocalizing eQTL+pQTL
  signals, plus neighbouring genes with non-colocalizing or distance-only
  signals. This is the known-answer fixture the tests assert against.
- **Colocalization** (`coloc.py`): self-contained Giambartolomei-2014
  `coloc.abf` (per-SNP ABF, PP.H0..H4). Validated against the known-answer
  behaviour (shared causal -> high PP.H4; distinct causal -> low PP.H4).
- **V2G** (`integrate/v2g.py`): per-gene coloc against every QTL layer/tissue,
  distance-to-TSS decay, optional functional-annotation feature; emits
  provenance-stamped `EvidenceItem`s.
- **Aggregation** (`evidence/aggregate.py`): Open-Targets harmonic-sum,
  normalized to `[0, 1]`, with per-source weights.
- **Confidence** (`evidence/confidence.py`): calibrated high/medium/low rule.
- **Ranking** (`score/prioritize.py`) and **dossier** (`report/dossier.py`).
- **CLI** (`cli.py`): `simulate`, `integrate`, `score`, `report`, `run-all`.

## Tests (all green)
- Known-answer: colocalizing eQTL+pQTL gene outranks distance-only neighbour.
- Harmonic-sum is monotonic in evidence and bounded `[0, 1]`.
- Every `EvidenceItem` carries non-null provenance (source, method, timestamp).
- Confidence rises with more corroborating layers.
- Coloc PP.H4 high for shared-causal, low for distinct-causal QTLs.
- Integration: full `simulate -> integrate -> score -> report` ranks the true
  causal gene first and writes an HTML dossier.

## Next slices (post-v0.1.0)
- Real-data loaders for GTEx / eQTL Catalogue / UK Biobank PPP summary stats.
- LD-aware coloc (SuSiE fine-mapping) for loci with multiple causal signals.
- Calibration against a labelled gold-standard V2G benchmark.
