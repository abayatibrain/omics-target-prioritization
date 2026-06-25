"""End-to-end statistical-genetics pipeline across all three repositories.

This script wires the three packages into one story, from raw genotypes to a
credentialed, directional drug target:

1. **gwas-meta-engine** simulates a biobank-style cohort *with LD blocks*, runs
   association testing, clumps the genome-wide hits into independent loci, and
   extracts the lead locus as summary statistics plus an in-sample LD matrix.
2. **post-gwas-causal** takes that locus, fine-maps it (SuSiE), and runs
   SuSiE-based colocalization against a cis-eQTL for a candidate gene to ask
   whether the disease acts *through* that gene.
3. **omics-target-prioritization** assembles the colocalization result and a
   drug-target Mendelian-randomization estimate (a Wald ratio on the lead cis
   instrument) into a calibrated target score and a one-page HTML dossier that
   states which direction to drug the target.

The eQTL is simulated on the *same* LD matrix the GWAS locus produced, sharing
the lead causal variant, so the pipeline has a known right answer (the gene
should colocalize). Everything else is real method code.

Run it (with all three packages installed):

    python examples/trio_end_to_end.py --out dossier.html

Requires the sibling packages ``gwas_meta_engine`` and ``post_gwas_causal`` to
be importable alongside ``omics_target_prioritization``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def run(seed: int, out: Path) -> None:
    # --- 1. GWAS: simulate with LD, associate, clump, extract a locus ----------
    from gwas_meta_engine.assoc.linear import linear_assoc
    from gwas_meta_engine.clump import clump, extract_locus
    from gwas_meta_engine.config import SimConfig
    from gwas_meta_engine.simulate import simulate

    cfg = SimConfig(
        n_samples=8000, n_snps=1000, n_causal=5, heritability=0.6, effect_sd=0.6,
        fst=0.02, strat_effect=0.1, ld_block_size=40, ld_rho=0.97, seed=seed,
    )
    sim = simulate(cfg)
    geno = sim.genotypes.astype(float)
    gwas = linear_assoc(geno, sim.phenotype_quant, sim.variant_info)
    loci = clump(gwas, geno, p_index=5e-8, r2_threshold=0.1, window_kb=200)
    if not loci:
        raise SystemExit("No genome-wide-significant loci; try another seed.")
    lead = loci[0]
    print(f"[1] GWAS -> {len(loci)} independent loci; lead {lead.lead_snp} (p={lead.lead_p:.1e})")

    locus = extract_locus(gwas, geno, lead.lead_idx, window_kb=120)
    ld = locus.ld
    z_gwas = locus.sumstats["Z"].to_numpy()
    lead_local = int(np.argmax(np.abs(z_gwas)))
    print(f"[1] extracted {len(z_gwas)}-SNP locus (lead at local index {lead_local})")

    # --- 2. Simulate a cis-eQTL sharing the GWAS lead causal variant -----------
    rng = np.random.default_rng(seed + 1)
    n_eqtl = 500
    b_eqtl = np.zeros(len(z_gwas))
    b_eqtl[lead_local] = 0.45
    chol = np.linalg.cholesky(ld + 1e-6 * np.eye(len(z_gwas)))
    z_eqtl = np.sqrt(n_eqtl) * (ld @ b_eqtl) + chol @ rng.standard_normal(len(z_gwas))

    # --- 3. Fine-mapping + SuSiE colocalization --------------------------------
    from post_gwas_causal.coloc.susie import coloc_susie
    from post_gwas_causal.finemap.susie import finemap_susie

    fm = finemap_susie(z_gwas, ld, cfg.n_samples)
    cs = coloc_susie(z_gwas, z_eqtl, ld, cfg.n_samples, n_eqtl, n_effects=5)
    print(
        f"[3] fine-mapping: {len(fm.credible_sets)} credible set(s); "
        f"coloc-SuSiE best PP.H4 = {cs.best_pp_h4:.3f}"
    )

    # --- 4. Drug-target MR (Wald ratio on the lead cis instrument) -------------
    beta_gwas = float(locus.sumstats["BETA"].to_numpy()[lead_local])
    se_gwas = float(locus.sumstats["SE"].to_numpy()[lead_local])
    beta_eqtl = float(z_eqtl[lead_local] / np.sqrt(n_eqtl))
    theta = beta_gwas / beta_eqtl
    se_theta = abs(se_gwas / beta_eqtl)
    print(f"[4] drug-target MR Wald ratio: theta = {theta:.3f} (se {se_theta:.3f})")

    # --- 5. Target prioritization + dossier ------------------------------------
    from omics_target_prioritization.evidence.mr_evidence import build_mr_evidence
    from omics_target_prioritization.models import EvidenceItem, Gene, Provenance
    from omics_target_prioritization.report.dossier import write_dossier
    from omics_target_prioritization.score.prioritize import score_target

    gene = Gene(
        gene_id="ENSG_DEMO",
        symbol="CANDGENE",
        chromosome="1",
        tss=int(locus.sumstats["BP"].to_numpy()[lead_local]),
    )
    coloc_ev = EvidenceItem(
        gene_id=gene.gene_id, source="eQTL_coloc", score=cs.best_pp_h4, weight=1.0,
        layer="eQTL", tissue="whole_blood",
        provenance=Provenance(
            dataset="sim_eQTL", method="coloc.susie", parameters={"PP.H4": cs.best_pp_h4}
        ),
    )
    mr_ev = build_mr_evidence(gene.gene_id, beta=theta, se=se_theta, exposure_layer="eQTL",
                              dataset="sim_eQTL")
    target = score_target(gene, [coloc_ev, mr_ev])
    print(
        f"[5] {gene.symbol}: score={target.total:.3f}  confidence={target.confidence}  "
        f"MR direction={target.mr_direction} (p={target.mr_pvalue:.1e})"
    )
    write_dossier(target, out, locus_id=lead.lead_snp)
    print(f"[5] wrote dossier -> {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=2024)
    parser.add_argument("--out", type=Path, default=Path("trio_dossier.html"))
    args = parser.parse_args()
    run(args.seed, args.out)


if __name__ == "__main__":
    main()
