# Are almost all graphs stable? No.

**Brown–Cameron Problem 2 has a negative answer.**

Brown and Cameron ([EJC 25(1), P1.46, 2018](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v25i1p46)) asked: *are almost all graphs stable?* That is, does the independence polynomial of $G(n,p)$ have all its roots in the closed left half-plane with probability tending to one? It does not.

## A certificate you can check by hand

Take the 86-vertex graph generated from the SHAKE256 seed `BC:86:622`. Counting its independent sets by size, starting with the empty set, gives

```
1, 86, 1810, 12478, 31706, 31559, 12792, 2225, 205, 21, 2
```

Its second Hurwitz determinant is

$$
21\cdot205 - 2\cdot2225 = -145 < 0,
$$

so its independence polynomial has a root in the open right half-plane. The roots are at about $2.199 \pm 9.008i$.

## The full picture

Let $u_n(p)$ be the probability that $G(n,p)$ has an independence root with positive real part.

| Edge probability | Behaviour of $u_n(p)$ as $n\to\infty$ | Status |
| --- | --- | --- |
| every $0<p<1$ | $\limsup u_n \ge 1-e^{-(1-p)^2}$, so the answer to Problem 2 is **no** | independently verified |
| $p = 1/2$ (uniform labelled graphs) | limit points fill $[0,\ 0.2714922345\ldots]$; average over orders $\to 0$ | independently verified |
| $p < 1-s_*^{1/4} \approx 0.2542775$ | $u_n \to 1$: almost all graphs are **unstable** ($s_*$ is the first spectral value of the partial theta function) | independently verified |
| all $p$ | five-phase classification with critical densities $p\approx0.273$ and $p\approx0.254$ | unrefereed, see [GAPS.md](GAPS.md) |

The mechanism is that the largest independent sets are rare: their count has a Poisson limit. Once the graph is rescaled, the top of the independence polynomial converges to a partial theta function, and one lucky large independent set pushes a pair of roots across the imaginary axis.

## Papers

1. **[Four-coefficient note](random-graph-stability.md)** ([PDF](random-graph-stability.pdf)): the negative answer to Problem 2, for every fixed $p$.
2. **[Follow-up](random-graph-stability-follow-up.md)** ([PDF](random-graph-stability-follow-up.pdf)): the exact set of limit points at $p=1/2$, and $u_n\to1$ for $p<0.254$.
3. **[Version 3](random-graph-stability-v3.md)** ([PDF](random-graph-stability-v3.pdf)): a limit law for the roots, the classification of every fixed edge probability, and unimodality of the independent-set sequence (Makowsky–Rakita Problem 1.1). Open points are listed in **[GAPS.md](GAPS.md)**.

## Verification status

| Document | What was checked independently |
| --- | --- |
| Four-coefficient note | Every step re-derived. The 86-vertex graph recounted by two separate methods. The right-half-plane roots certified by an exact Routh table. The bound at 260,000 vertices recomputed in exact rational arithmetic. |
| Follow-up | Sections 1–6 re-derived. The constants $\xi,\tau,L$ enclosed with separate rational code. The 65.6M-vertex bound and the 8,192-graph sample reproduced exactly. |
| Version 3 | The verifier passes. The predicted stable band of the limiting function was cross-checked with exact Routh counts on true finite-size coefficients up to degree 400. Sections 6–8 (boundary perturbation, adjacent orders, critical averaging rates) still need an expert referee. |

## Reproduce

Python 3.8+ and the standard library are enough. Run these from the repository root.

```sh
python3 verify_random_stability.py              # four-coefficient note, under 1 s
python3 -B verify_stability_followup.py         # follow-up, under 1 s
python3 -B verify_stability_followup.py --sample-86 8192   # optional sample, a few minutes
python3 -B verify_v3.py                         # version 3, about 1 min
```

## Authorship and disclosure

The mathematics and the verifiers were generated with OpenAI Codex, prompted by Jerald Hu Yuan. An independent AI referee (Claude) re-derived the proofs and reproduced the computations with separate code. None of this has been peer reviewed. Corrections and expert review are very welcome: please open an issue.

## Reference

J. I. Brown and B. Cameron, *On the Stability of Independence Polynomials*, Electron. J. Combin. **25**(1) (2018), #P1.46, Problem 2. [doi:10.37236/7280](https://doi.org/10.37236/7280)
