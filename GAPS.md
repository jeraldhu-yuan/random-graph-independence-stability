# Scope and referee questions for version 3

The checked base is 423f7949c981ca39a3721e6b3d98acfa94cedbe6. The reported independent check of that revision covers F2. The new phase diagram and its endpoint arguments still await independent mathematical review. The finite verifier checks the analytic inequalities it lists and the reproduced finite computations. It is not a formalization of the probability proof.

## Corrections to the proposed statement

- **Actual versus padded zeros, (2).** If the target measure on the whole plane were simply $Z(T+B_q)$, its atom at zero when $T=0$ would be spurious. At $\mu_k\to0$, with probability tending to one the degree of $I_G$ is $k-1$ and no actual reciprocal root is zero. Equation (2) subtracts the padded atom. Boundary-free count convergence or vague convergence is needed when a compact set cuts through a limiting zero.
- **The first-instability parameter, (3).** The equivalence between its vanishing and $q^4>s_*$ fails. At $q=37/50$ the stable interval is strictly inside $(0,1)$, although the odd part has only real squared-variable roots. At $Q$ a single stable constant survives. The even part at the first two odd zeros distinguishes these cases.

## Steps deserving a fresh derivation

- **Exterior control, (6).** The step from the full mean Laurent polynomial to a random polynomial with one exceptional coefficient is needed to exclude roots coming from infinity. A compact-disk argument alone would leave a gap. Check the direction of the $q$-spacing induction and the relative perturbation estimate for arbitrarily large radius and bounded $t$. The separate four-coefficient test handles unbounded $t$.
- **All higher odd zeros, (9)–(10).** A missed annulus could produce another active crossing. Check the winding number on both bounding circles, the two missing bilateral tails and the sign at $x/\sqrt{s}$ after shifting the central index. Independently evaluating the compact boxes with a rational or directed-rounding library would check the machine-assisted part.
- **Positive boundary mass, (12)–(13).** At an integer Poisson threshold a weak-limit inequality alone can lose a whole atom. The proof uses the sign of the $1/k$ displacement, including simultaneous atoms at opposite ends of a band. At $Q$ the odd zero is double but the entire-function zero is simple. The claimed inner and outer $k^{-1/2}$ bands depend on both these facts. A referee should check the expansion before treating a closed endpoint as attained.
- **Passage between consecutive orders, (14)–(15).** The count $X_k$ and its neighbor $X_{k-1}$ need not be independent. The proof uses unconditional anti-concentration of the latter around finitely many deterministic centers, coupled increments, and a smaller error from every lower count. The test here is whether $kV_{k-2}=o(1/\mu_{k-1})$ and the positive even-part slope suffice at the tangency as well as at both simple crossings. Removing this step would give extrema without proving the full intermediate interval.

## Finer questions left open

- **Individual resonant sequences, §7.** The table gives the entire set of subsequential limits. It does not assign a unique instability limit using only a resonant value of $\lambda$. A law on the finer approach scale, including the joint fluctuations of the adjacent counts, would resolve individual tuning sequences.
- **Decay above the spectral point, last row of Theorem B.** The bound $1-\bar u_N=O_q(N^{-1/4})$ proves convergence of the average to one. The sharp rate of the deficit in this regime is not established. It would require a sharper estimate for stable exceptions, including atypically large independent sets.
- **Varying densities and leading constants, §8.** The constants are for fixed $p$. The paper gives two-sided order estimates in the first four rows, not leading constants for the block sums. Simultaneous limits $q=q(n)\to q_c$ or $Q$ require uniform control of shrinking derivative margins and of the position within the last graph-order block.
- **Priority.** The literature section compares identifiable primary results with the new claims. It does not certify exhaustive priority. A reviewer should check for independent work on analytic limits of random clique or independence polynomials under reciprocal scaling.

## Prompt for an independent reviewer

> Act as a skeptical mathematical referee. Start from the version 3 commit in this repository and treat only the older follow-up at 423f794 as independently checked. Reproduce the rational certificates with a separate implementation. Try to disprove Theorem B at the two critical densities, checking the all-degree bound (6), higher odd zeros (9)–(10), displacement of endpoint atoms (12)–(13), and the coupled adjacent-order argument (14)–(15). Check the root measure in (2), verify the hypotheses of the cited clique-count results in the original sources, and check the final unimodality argument separately. Return a verdict for each theorem, the first unsupported inference with a precise location if one exists, and any counterexample or repair. Identify which parts you verified independently of the supplied Python program.
