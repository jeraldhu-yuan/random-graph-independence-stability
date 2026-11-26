# Proof boundaries and open questions

This is a research preprint. The finite verifier covers the numerical enclosure and enumerations. The asymptotic claims rely on the arguments in the manuscript and the cited theorems. An independent mathematical review of Theorem 1 is outstanding.

| Claim or boundary | What deserves scrutiny | How to settle it |
|---|---|---|
| Equation (5), $\sup_T \lvert M(T)-\mu(T)\rvert=o(\sqrt n)$ | The September preprint states a distribution-function CLT. Deducing a statement about a maximizer additionally needs Fourier decay uniform in the tree. | Check Proposition 5.1, Lemma 6.1 and equation (6.6) at activity one. Reproduce the uniform Fourier inversion and the comparison with the nearest integer. |
| Equations (9) and (10), $\mathbb E\mu(T_n)-cn=O((1+\log n)^4)$ | Using the ordinary unconditioned Poisson tree for the sampled root, or using fixed-radius convergence at a growing radius without an error bound, would give invalid centering. | Recount rooted balls with the boundary attachment restriction. Check both probability tails at $h=\lceil3\log_2 n\rceil$ and the uniform moment hypotheses. |
| Equations (11)–(13), centered variance approximation | An error small compared with $n$ can be large compared with $\sqrt n$. A random endpoint of a random edge is degree biased. | Derive the refresh inequality on the independent-coordinate sampler. Check the two distance ranges in (12), the degree weight, and the order of limits in $D,h,n$. |
| Section 5, $\operatorname{Var}F_{h,D}=n v_{h,D}+O(1)$ and its CLT | Intersecting balls can impose incompatible degree constraints. A high-degree vertex can affect many terms if the finite-size cutoff is dropped. | Check the identification argument in (14), repeated centers, the two-switch condition on the union graph, and the proof for $v_{h,D}=0$. |
| Section 6, $v\ge e^{-10}/216$ | Preserving only root occupation odds need not give an outside-independent shift of the mean. Conditional slot choices must have the claimed labelled multiplicities. | Verify the shift in **both** root states. Check the conditioning after rooting at label 1, the disjointness of equal-size slots, and (15). |
| Section 7, $0.30612549\le c\le0.30612646$ | The interpolation estimate uses unbiased splitting with unrounded ancestors. Mixing up the spine and ordinary messages, or confusing probability-mass units and grid-coordinate units, could invalidate the enclosure. | Recompute the recursion independently, derive the curvature bound along one path, and bound all update errors in total variation before rounding outward. |

## Quantitative questions still open in this manuscript

- Equation (2) gives a wide enclosure for $v$. A useful numerical value needs a covariance computation with an explicit truncation error, or another identity for the variance.
- Equation (5) leaves the size of the mode–mean difference below scale $\sqrt n$ undetermined. A uniform local expansion with sufficient accuracy to compare adjacent lattice probabilities could yield an order-one bound or a description of ties.
- Equation (1) gives convergence without a rate. A rate would require quantitative estimates through the two cutoffs, and a quantitative mode–mean comparison.
- Corollary 2 describes $M_n/\alpha(T_n)$ to first order. A fluctuation law for this ratio would require joint control of the mode and the independence number, including their covariance.

The finite-forest unimodality results and the fixed-forest analytic estimates are cited inputs. The accompanying verifier checks the finite claims specific to this manuscript. The external preprints' full proofs and computational certificates have not been independently audited as part of this work. [Fang et al.](https://arxiv.org/pdf/2609.20961) · [Li et al.](https://arxiv.org/pdf/2610.07943)

[Manuscript](manuscript.md) · [Verifier](verify.py) · [Independent review prompt](independent-review-prompt.md)
