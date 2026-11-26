# Independent review prompt

Review the research preprint **Gaussian fluctuations of the independence-sequence mode of a random labelled tree** in the public repository [jeraldhu-yuan/random-graph-independence-stability](https://github.com/jeraldhu-yuan/random-graph-independence-stability/tree/main/random-tree-mode). Record the exact Git commit reviewed. The tree study is in `random-tree-mode/`. The older papers at the repository root concern a different probability model.

Start with arXiv:2609.20961 and 2610.07943. Identify which law is conditional on a fixed tree and which varies over uniform labelled trees. Verify the status of unimodality, the part of the first-order mode law implied by prior work, and whether any primary source already proves the claimed tree-to-tree CLT. Check Heilman's published Question 8.2 and the exact versions of every theorem used. [Fang et al.](https://arxiv.org/pdf/2609.20961) · [Li et al.](https://arxiv.org/pdf/2610.07943) · [Heilman](https://doi.org/10.1002/jgt.23225)

Audit this assertion for the **largest mode** $M_n$ of the independent-set sequence of a uniform labelled tree: $(M_n-cn)/\sqrt n$ converges in distribution to a centered Gaussian with variance $v>0$, and $\operatorname{Var}(M_n)/n\to v$. Read `manuscript.md` and `GAPS.md` and run `python3 verify.py`. Look for counterexamples and incorrect uniformity or centering before accepting the proof.

1. Derive the fixed-tree mode–mean estimate from the stated analytic inputs. Check arbitrary tied modes, lattice Fourier tails, and the implication for variances.
2. Recompute the root law and centering. Distinguish the finite ordinary Galton–Watson tree from the one-spine law. Check the ball ratio with a growing radius and the moment tails.
3. Reprove the independent-coordinate cut inequality. Check degree weighting, high-degree trees, both distance ranges in the truncation bound and the centered $L^2$ approximation.
4. Verify the fixed-radius, fixed-size functional satisfies all conditions of the Isaev–Southwell–Zhukovskii theorem. Derive the variance expansion for overlapping and disjoint balls. Check vanishing component variance and the passage to the limit.
5. Enumerate the two rooted shapes independently. Check their conditional mean shifts, automorphism factors, number of eligible labelled blocks, and independence after conditioning.
6. Reimplement the density calculation, preferably with a different numerical representation. Check the truncation error, random interpolation error, explicit infinite-series tails, denominator corrections, and enclosure for $c/\omega$. Treat the seeded examples as simulations.

Return a verdict for each numbered item and the main theorem. For any failure cite the exact sentence or equation, give a concrete obstruction or counterexample, and identify the strongest surviving result. Separate results that follow from earlier papers from the contributions of this manuscript. State which computations you independently reproduced and which published inputs you took as given.
