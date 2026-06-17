"""Self-contained Giambartolomei-2014 ``coloc.abf`` colocalization.

This module implements the approximate-Bayes-factor (ABF) colocalization test
of Giambartolomei et al. (2014, *PLoS Genet* 10(5):e1004383). Given two sets of
per-SNP summary statistics over the *same* ordered set of variants — here, a
GWAS trait and a molecular QTL — it returns the five posterior probabilities:

    - ``PP.H0``: no causal variant for either trait,
    - ``PP.H1``: a causal variant for trait 1 only,
    - ``PP.H2``: a causal variant for trait 2 only,
    - ``PP.H3``: distinct causal variants (two signals in LD),
    - ``PP.H4``: a single *shared* causal variant (colocalization).

``PP.H4`` is the quantity that drives variant-to-gene scoring: a high ``PP.H4``
means the disease signal and the gene's molecular-trait signal are best
explained by the *same* causal variant, i.e. the variant plausibly acts through
that gene.

The implementation follows the original derivation exactly:

1. Per-SNP Wakefield approximate Bayes factor (ABF) from ``beta`` and ``se``,
   under a normal prior ``N(0, W)`` on the true effect size.
2. Per-hypothesis prior probabilities ``p1, p2, p12``.
3. Posterior combination over all SNP configurations, done in log-space with
   the log-sum-exp trick for numerical stability.

No external coloc library is used; ``numpy``/``scipy`` only.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

from omics_target_prioritization.models import GwasLocus, QtlAssociation


@dataclass(frozen=True)
class ColocPriors:
    """Prior probabilities for the five coloc hypotheses.

    Parameters
    ----------
    p1
        Prior probability that a given SNP is causal for trait 1 only.
    p2
        Prior probability that a given SNP is causal for trait 2 only.
    p12
        Prior probability that a given SNP is causal for *both* traits.
    sd_prior
        Standard deviation ``sqrt(W)`` of the normal prior on the true effect
        size used in the Wakefield ABF. The coloc default is ``0.15`` for
        quantitative traits and ``0.2`` for case/control on the log-OR scale;
        ``0.15`` is a reasonable shared default for summary-stat coloc.
    """

    p1: float = 1e-4
    p2: float = 1e-4
    p12: float = 1e-5
    sd_prior: float = 0.15


@dataclass(frozen=True)
class ColocResult:
    """Posterior probabilities from a coloc.abf run.

    Attributes
    ----------
    pp_h0, pp_h1, pp_h2, pp_h3, pp_h4
        Posterior probabilities of the five hypotheses; they sum to 1.
    n_snps
        Number of SNPs the test was run over.
    """

    pp_h0: float
    pp_h1: float
    pp_h2: float
    pp_h3: float
    pp_h4: float
    n_snps: int

    def as_dict(self) -> dict[str, float]:
        """Return the posteriors as a plain dict keyed ``H0``..``H4``."""
        return {
            "H0": self.pp_h0,
            "H1": self.pp_h1,
            "H2": self.pp_h2,
            "H3": self.pp_h3,
            "H4": self.pp_h4,
        }


def _log_abf(beta: np.ndarray, se: np.ndarray, sd_prior: float) -> np.ndarray:
    r"""Per-SNP log Wakefield approximate Bayes factor.

    For a marginal estimate ``beta`` with standard error ``se`` and a prior
    ``N(0, W)`` on the true effect (``W = sd_prior**2``), the Wakefield ABF for
    "this SNP is causal" vs "null" is

    .. math::

        \mathrm{ABF} = \sqrt{\frac{V}{V + W}}\,
                       \exp\!\left(\frac{z^2}{2}\,\frac{W}{V + W}\right),

    where ``V = se**2`` is the variance of the estimate and ``z = beta / se``.

    Parameters
    ----------
    beta
        Marginal effect estimates, shape ``(n_snps,)``.
    se
        Standard errors, shape ``(n_snps,)``.
    sd_prior
        Prior standard deviation on the true effect size.

    Returns
    -------
    numpy.ndarray
        Log ABF per SNP, shape ``(n_snps,)``.
    """
    v = se**2
    w = sd_prior**2
    z2 = (beta / se) ** 2
    # log ABF = 0.5*log(V/(V+W)) + 0.5*z^2*W/(V+W)
    return 0.5 * np.log(v / (v + w)) + 0.5 * z2 * (w / (v + w))


def coloc_abf(
    log_abf1: np.ndarray,
    log_abf2: np.ndarray,
    priors: ColocPriors | None = None,
) -> ColocResult:
    r"""Colocalization posterior probabilities from two log-ABF vectors.

    Implements the configuration-sum of Giambartolomei et al. (2014). For each
    hypothesis we sum, over all admissible SNP configurations, the product of
    per-SNP ABFs times the appropriate prior, working entirely in log-space.

    The per-hypothesis (un-normalized) log evidences are:

    - ``H0``: ``0`` (no causal SNP; ABF = 1).
    - ``H1``: ``log p1 + logsumexp(log_abf1)``.
    - ``H2``: ``log p2 + logsumexp(log_abf2)``.
    - ``H3``: ``log p1 + log p2 + log(S1 * S2 - sum_i abf1_i abf2_i)`` where the
      subtraction removes the *same-SNP* configurations (those belong to H4).
    - ``H4``: ``log p12 + logsumexp(log_abf1 + log_abf2)`` (shared SNP).

    Parameters
    ----------
    log_abf1, log_abf2
        Per-SNP log ABFs for trait 1 and trait 2, same length and SNP order.
    priors
        Hypothesis priors; defaults to :class:`ColocPriors`.

    Returns
    -------
    ColocResult
        Normalized posterior probabilities summing to 1.
    """
    if priors is None:
        priors = ColocPriors()
    if log_abf1.shape != log_abf2.shape:
        raise ValueError("log_abf1 and log_abf2 must have the same shape")
    n = int(log_abf1.shape[0])
    if n == 0:
        raise ValueError("coloc requires at least one SNP")

    log_p1 = np.log(priors.p1)
    log_p2 = np.log(priors.p2)
    log_p12 = np.log(priors.p12)

    # H1 / H2: exactly one causal SNP for one trait.
    l_h1 = log_p1 + logsumexp(log_abf1)
    l_h2 = log_p2 + logsumexp(log_abf2)

    # H4: one shared causal SNP -> sum of per-SNP products of ABFs.
    l_h4 = log_p12 + logsumexp(log_abf1 + log_abf2)

    # H3: distinct causal SNPs -> (sum_i abf1_i)(sum_j abf2_j) over i != j.
    # In log-space: log(S1 * S2 - sum_i abf1_i abf2_i). The bracketed term is
    # the full outer-product sum minus the diagonal (which is the H4 sum).
    s1 = logsumexp(log_abf1)
    s2 = logsumexp(log_abf2)
    log_full = s1 + s2  # log of S1 * S2
    log_diag = logsumexp(log_abf1 + log_abf2)  # log of diagonal sum
    # log(exp(log_full) - exp(log_diag)); log_full >= log_diag always.
    diff = log_full + np.log1p(-np.exp(np.minimum(log_diag - log_full, -1e-300)))
    l_h3 = log_p1 + log_p2 + diff

    log_evidences = np.array([0.0, l_h1, l_h2, l_h3, l_h4])
    log_norm = logsumexp(log_evidences)
    posteriors = np.exp(log_evidences - log_norm)

    return ColocResult(
        pp_h0=float(posteriors[0]),
        pp_h1=float(posteriors[1]),
        pp_h2=float(posteriors[2]),
        pp_h3=float(posteriors[3]),
        pp_h4=float(posteriors[4]),
        n_snps=n,
    )


def colocalize(
    gwas: GwasLocus,
    qtl: QtlAssociation,
    priors: ColocPriors | None = None,
) -> ColocResult:
    """Run coloc.abf between a GWAS locus and a QTL dataset.

    The two summary-statistic vectors must describe the *same* ordered set of
    variants. The simulator and real-data loaders are responsible for harmonizing
    SNP order before calling this function.

    Parameters
    ----------
    gwas
        The GWAS locus.
    qtl
        The QTL dataset (one gene, one tissue, one layer).
    priors
        Coloc priors; defaults to :class:`ColocPriors`.

    Returns
    -------
    ColocResult
        Posterior probabilities, with ``pp_h4`` the colocalization posterior.

    Raises
    ------
    ValueError
        If the GWAS and QTL summary stats are not the same length.
    """
    if priors is None:
        priors = ColocPriors()
    if len(gwas.summary_stats) != len(qtl.summary_stats):
        raise ValueError(
            "GWAS and QTL summary stats must cover the same SNPs "
            f"(got {len(gwas.summary_stats)} vs {len(qtl.summary_stats)})"
        )

    beta1 = np.array([s.beta for s in gwas.summary_stats], dtype=float)
    se1 = np.array([s.se for s in gwas.summary_stats], dtype=float)
    beta2 = np.array([s.beta for s in qtl.summary_stats], dtype=float)
    se2 = np.array([s.se for s in qtl.summary_stats], dtype=float)

    log_abf1 = _log_abf(beta1, se1, priors.sd_prior)
    log_abf2 = _log_abf(beta2, se2, priors.sd_prior)
    return coloc_abf(log_abf1, log_abf2, priors)
