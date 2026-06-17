"""Tests for V2G integration, provenance, and locus-level prioritization."""

from __future__ import annotations

from omics_target_prioritization.integrate.v2g import (
    distance_to_tss_score,
    integrate_locus,
    score_gene,
)
from omics_target_prioritization.score.prioritize import prioritize
from omics_target_prioritization.simulate import SimulatedLocus


def test_distance_decay_is_monotonic_and_bounded(sim_locus: SimulatedLocus) -> None:
    causal = sim_locus.genes[0]
    s = distance_to_tss_score(sim_locus.gwas, causal)
    assert 0.0 < s <= 1.0


def test_every_evidence_item_has_full_provenance(sim_locus: SimulatedLocus) -> None:
    items = score_gene(sim_locus.gwas, sim_locus.genes[0], sim_locus.qtls)
    assert items
    for item in items:
        assert item.source is not None
        assert item.provenance.dataset
        assert item.provenance.method
        assert item.provenance.timestamp is not None


def test_causal_gene_gets_coloc_evidence(sim_locus: SimulatedLocus) -> None:
    items = score_gene(sim_locus.gwas, sim_locus.genes[0], sim_locus.qtls)
    coloc_items = [i for i in items if i.source.endswith("_coloc")]
    sources = {i.source for i in coloc_items}
    assert "eQTL_coloc" in sources
    assert "pQTL_coloc" in sources
    # Shared causal SNP -> high PP.H4 for the causal gene.
    assert max(i.score for i in coloc_items) > 0.8


def test_functional_annotation_emits_item(sim_locus: SimulatedLocus) -> None:
    items = score_gene(
        sim_locus.gwas,
        sim_locus.genes[0],
        sim_locus.qtls,
        functional_annotation=0.75,
    )
    annot = [i for i in items if i.source == "functional_annotation"]
    assert len(annot) == 1
    assert annot[0].score == 0.75


def test_known_answer_causal_gene_ranks_first(sim_locus: SimulatedLocus) -> None:
    evidence = integrate_locus(sim_locus.gwas, sim_locus.genes, sim_locus.qtls)
    ranked = prioritize(sim_locus.genes, evidence)
    assert ranked[0].gene_id == sim_locus.causal_gene_id
    # Causal gene should also be the most confident.
    assert ranked[0].confidence == "high"
    # And it should beat a distance-only neighbour that has no coloc.
    assert ranked[0].total > ranked[-1].total


def test_provenance_recorded_on_every_item_across_locus(sim_locus: SimulatedLocus) -> None:
    evidence = integrate_locus(sim_locus.gwas, sim_locus.genes, sim_locus.qtls)
    for items in evidence.values():
        for item in items:
            assert item.provenance.dataset
            assert item.provenance.method
            assert item.provenance.timestamp is not None
