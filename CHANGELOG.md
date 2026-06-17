# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] — 2026-06-17

### Added
- **`models.py`** — pydantic v2 data models for the full evidence chain:
  `Variant`, `Gene`, `QtlAssociation` (eQTL/pQTL/sQTL/caQTL), `GwasLocus`,
  `Provenance`, `EvidenceItem`, and `TargetScore`.
- **`simulate.py`** — deterministic simulator (`SimConfig`) that plants a true
  causal gene with colocalizing eQTL+pQTL signals and surrounds it with
  neighbouring genes carrying non-colocalizing or distance-only signals.
- **`coloc.py`** — self-contained Giambartolomei-2014 `coloc.abf`
  implementation (per-SNP approximate Bayes factors, posteriors PP.H0..H4).
- **`integrate/v2g.py`** — colocalization-based variant-to-gene scoring with a
  distance-to-TSS decay feature and optional functional-annotation feature;
  emits fully provenance-tracked `EvidenceItem`s.
- **`evidence/aggregate.py`** — Open-Targets-style weighted harmonic-sum
  aggregation, normalized to `[0, 1]`.
- **`evidence/confidence.py`** — calibrated high/medium/low confidence label
  from corroborating omics layers, max PP.H4, and aggregate score.
- **`score/prioritize.py`** — locus-level ranking returning ordered
  `TargetScore`s.
- **`report/dossier.py`** — Jinja2 HTML one-pager per prioritized target.
- **`cli.py`** — Typer CLI: `simulate`, `integrate`, `score`, `report`,
  `run-all`.
- Unit + integration test suite; ADRs 0001–0003; mkdocs documentation.

[0.1.0]: https://github.com/abayatibrain/omics-target-prioritization/releases/tag/v0.1.0
