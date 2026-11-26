# Gaussian fluctuations of the independence-sequence mode of a random labelled tree

Research preprint · 26 November 2026

## 1. The mode and the two sources of randomness

Choose a tree $T_n$ uniformly from the labelled trees on $[n]$. Write $i_k(T)$ for the number of independent sets of size $k$, and $M(T)$ for the largest index maximizing $i_k(T)$. Conditional on $T$, choose one of its independent sets uniformly and write $K$ for its size. Thus $\mu(T)=\mathbb E[K\mid T]=\sum_k k i_k(T)/\sum_k i_k(T)$. The mode is a function of the tree alone. Heilman's Question 8.2 asks whether it concentrates for this tree model. [Heilman, *Journal of Graph Theory* 109, Question 8.2](https://doi.org/10.1002/jgt.23225)

**Theorem 1 (fluctuations over the choice of tree).** There are constants $c$ and $v>0$ such that, for $M_n=M(T_n)$,

$$
\frac{M_n-cn}{\sqrt n}\ \xrightarrow{\ d\ }\ \mathcal N(0,v),
\qquad \frac{\operatorname{Var}(M_n)}{n}\longrightarrow v .
\tag{1}
$$

The same assertions hold with $\mu(T_n)$ in place of $M_n$. Every choice among tied modes has the same limits. The constant $c$ is the expected probability that the root is occupied under the activity-one hard-core measure on the one-spine Poisson(1) Galton–Watson tree. We obtain the rational enclosures

$$
0.30612549\le c\le0.30612646,\qquad
\frac1{12754584}<\frac{e^{-10}}{216}\le v\le384 .
\tag{2}
$$

In particular $M_n/n\to c$ in probability. A window $cn\pm a_n\sqrt n$ contains the mode with probability tending to one whenever $a_n\to\infty$. For fixed $a>0$ the limiting probability is $2\Phi(a/\sqrt v)-1$, where $\Phi$ is the standard normal distribution function. These statements follow from (1) by tightness and continuity. The variance enclosure certifies nondegeneracy. It is too wide to predict a practical confidence interval.

**Corollary 2 (location within the sequence).** Let $\alpha(T)$ be the maximum independent-set size. Then $M_n/\alpha(T_n)\to c/\omega$ in probability, where $\omega=e^{-\omega}$ and $0.539767\le c/\omega\le0.539770$. Indeed, the independence-number law gives $\alpha(T_n)/n\to\omega$ by taking offspring generating function $f(s)=e^{s-1}$ in Theorem 2 of Devroye, Goh and Zhao. The denominator counts ranks up to the largest nonzero coefficient. [Devroye, Goh and Zhao, Theorem 2](https://luc.devroye.org/Devroye%2BGoh%2BZhao-The_peel_number_and_the_leaf_height_of_a_Galton_Watson_tree-2021-2022.pdf)

**Corollary 3 (the coefficients near the peak).** Put $\sigma_T^2=\operatorname{Var}(K\mid T)$ and $\phi(x)=e^{-x^2/2}/\sqrt{2\pi}$. Uniformly over all $n$-vertex trees and all integers $k$, with $i_k=0$ outside its support,

$$
\frac{\sigma_T i_k(T)}{\sum_j i_j(T)}
=\phi\!\left(\frac{k-M(T)}{\sigma_T}\right)+o(1).
\tag{3}
$$

This is a consequence of the fixed-forest estimates in §2. Their strict log-concavity theorem also applies throughout $|k-cn|\le n/100$ for every sufficiently large tree. Indeed bipartition gives $\alpha(T)\ge n/2$, and (2) gives $c-1/100>1/5$ and $c+1/100<17/50$. With probability tending to one this interval contains all modes, so there can be at most two, consecutive in case of a tie. A strictly decreasing sequence of successive coefficient ratios can equal one at most once.

## 2. What was known about this question

Fang, Lu, Nevo, Yao and Zheng prove unimodality for every sufficiently large forest (Theorem 1.1) and strict log-concavity for $n/5\le k\le17\alpha/25$ (Theorem 1.2). Li, Vallier and Zhang give a proof for every finite forest (Theorem 1.1), crediting the earlier all-orders result of Zhang and Li. Thus unimodality has probability one for each $n$. Their Proposition 2.3 gives a deterministic location bracket. With $\nu$ the maximum matching size, the sequence increases weakly through $\lceil n/4\rceil$ and decreases weakly from $h_B\le\lceil n/2-\nu/3\rceil$, where $h_B$ is the first mode of the coefficient sequence of $(1+2x)^\nu(1+x)^{n-2\nu}$. Earlier, Basit and Galvin obtained high-probability monotonicity in an initial and a terminal segment for uniform labelled trees. [Fang et al., §§1, 8](https://arxiv.org/pdf/2609.20961) · [Li et al., §§1, 2.3](https://arxiv.org/pdf/2610.07943) · [Basit and Galvin, Theorems 1.4 and 1.7 in the preprint](https://arxiv.org/pdf/2006.12562)

Here is the mode-location implication of the September preprint. At activity one its Proposition 5.1 gives $\sigma_T^2\asymp n$ and a uniform central limit theorem for $(K-\mu(T))/\sigma_T$. Its Lemma 6.1 gives $|\mathbb E[e^{itK}\mid T]|\le e^{-b n\sin^2(t/2)}$ for an absolute $b>0$. Let $\chi_T$ be the characteristic function of the standardized variable. On $|u|\le\pi\sigma_T$ we have $|\chi_T(u)|\le e^{-b' u^2}$ for an absolute $b'>0$, using concavity of sine on $[0,\pi/2]$ and the variance bound. Extend the integrand by zero outside this interval. Along every sequence of forests, the CLT gives pointwise convergence to $e^{-u^2/2}$, and the bound gives convergence in $L^1$. Lattice Fourier inversion and the Gaussian Fourier formula now give

$$
\sup_{T,k}\left|\sigma_T\Pr(K=k\mid T)
-\phi\!\left((k-\mu(T))/\sigma_T\right)\right|\longrightarrow0.
\tag{4}
$$

Indeed the error is at most the $L^1$ distance divided by $2\pi$, independently of the Fourier phase. If uniformity failed one could pick a contradicting sequence of forests. Compare any mode to the nearest integer to $\mu(T)$. On $|k-\mu(T)|\ge\epsilon\sigma_T$ the Gaussian falls below its maximum by a fixed positive amount. Equation (4) therefore gives, for all choices of mode,

$$
\sup_T |M(T)-\mu(T)|/\sqrt n\ \longrightarrow\ 0 .
\tag{5}
$$

Uniform continuity of $\phi$ proves (3). The fluctuations in (4) arise from sampling an independent set after fixing its tree. Neither forest preprint gives a fluctuation law over the sampling of $T_n$ itself. [Fang et al., Proposition 5.1, Lemma 6.1 and equation (6.6)](https://arxiv.org/pdf/2609.20961) · [Li et al.](https://arxiv.org/pdf/2610.07943)

The first-order description of $c$ also follows from earlier ingredients. Janson's Theorem 7.1 gives the rooted local limit in the Poisson case. For each fixed ball shape its empirical frequency concentrates by the pattern variance estimate of Isaev, Southwell and Zhukovskii, Lemma 2.4 (distinguishing a center changes the fixed multiplicity). The hard-core root marginal follows from the tree recursion in equation (2.9) of Dembo, Montanari and Sun, together with loss of the boundary influence on the one-spine tree. For fixed radius one may first retain a finite collection of ball shapes of probability arbitrarily close to one. The omitted empirical mass then has arbitrarily small mean. Letting the radius grow yields the occupation-density law and (5) transfers it to the mode. We give a quantitative centering estimate in §3 because an error $o(n)$ alone is insufficient for (1). Ralaivaosaona, Šileikis and Wagner prove a CLT for $\log\sum_k i_k(T_n)$ (Theorem 10). Their argument with exchangeable fringe subtrees motivates §6. The derivative $\mu(T_n)$ and the maximizing index require additional estimates. [Janson, Theorem 7.1 and Example 9.2](https://arxiv.org/pdf/1112.0510) · [Isaev, Southwell and Zhukovskii](https://doi.org/10.1017/S0963548321000523) · [Dembo, Montanari and Sun, §§2.3–2.4](https://arxiv.org/pdf/1110.4821) · [Ralaivaosaona, Šileikis and Wagner, §6.1](https://arxiv.org/pdf/1810.00467)

## 3. The density and centering to smaller than the fluctuation scale

Root a finite tree at the vertex of interest and let $U$ be the probability that its root is vacant. The weights for excluding and including that vertex are respectively the product of the child partition functions and that product times $\prod_i U_i$. Set $x=-\log_2 U$. Its recursion and the occupation probability are

$$
x=g\!\left(\sum_i x_i\right),\qquad
g(s)=\log_2(1+2^{-s}),\qquad b(x)=1-2^{-x}.
\tag{6}
$$

All natural messages lie in $[0,1]$. For $s\ge0$, $|g'(s)|=1/(1+2^s)\le1/2$ and $|g''(s)|=(\log 2)2^s/(1+2^s)^2\le1/4$. Also $|b'|,|b''|\le1$ on $[0,1]$. Write $Z_j(v)$ for the number of vertices at distance $j$, $B_h(v)$ for the closed radius-$h$ ball, and $s_h(T,v)$ for the root occupation probability computed on the finite induced ball with free boundary. The full marginal is $s(T,v)$. Varying boundary messages one at a time and multiplying derivatives along their unique paths gives

$$
|s(T,v)-s_h(T,v)|\le 2^{-h} Z_h(v).
\tag{7}
$$

The bound applies to any two assignments in $[0,1]$ at depth $h$. On the limiting tree $\mathcal K$ there is one infinite path starting at the root. At each of its vertices an independent Poisson(1) number of ordinary critical Poisson bushes branches off. Such bushes are finite almost surely (the extinction probability is the least fixed point of $e^{z-1}$, which has only the root $1$ in $[0,1]$, as follows by iteration from zero). The mean size of each ordinary generation is one, so $\mathbb E Z_h(\mathcal K)=h+1$. Free balls form a Cauchy sequence of root marginals almost surely, since the bound in (7) tends to zero almost surely by summability of its expectations. Arbitrary admissible pinned boundary conditions one level farther out give messages in $[0,1]$ one level inward because a pinned occupied child forces vacancy, and otherwise the occupation odds of its parent are at most one. This also gives uniqueness of the local marginal and, by the same argument on finite connected sets (condition on their boundary messages), uniqueness of the Gibbs law. Define $c=\mathbb E s(\mathcal K,o)$. Section 7 encloses this expectation.

We use the following moment estimates uniformly in $n$. For each fixed positive integer $r$, with $V$ a uniform vertex of $T_n$, and also for the root of $\mathcal K$,

$$
\|Z_j(V)\|_r=O_r(j+1),\qquad \||B_h(V)|\|_r=O_r((h+1)^2).
\tag{8}
$$

The first is Lemma 3 of Ralaivaosaona, Šileikis and Wagner, with their moment parameter chosen large enough, possible for Poisson offspring. For a fixed root label, conditioning a Poisson(1) tree to have $n$ vertices gives the rooted Cayley shape. Indeed, summing the plane-tree probabilities $e^{-n}/\prod_v d_v!$ over child orders gives $e^{-n}/|\operatorname{Aut}(H)|$ for a rooted unordered shape $H$, proportional to its number of root-preserving labellings. Here $d_v$ counts children. The same law holds for a uniform root. Summing the first estimate by Minkowski gives the second. [Ralaivaosaona, Šileikis and Wagner, Lemma 3](https://arxiv.org/pdf/1810.00467)

For a radius-$h$ ball shape $H$ with $m$ vertices, $b_H$ boundary vertices at level $h$ and $a_H$ rooted automorphisms, let $P_n(H)$ and $P_\infty(H)$ be its probabilities at these two roots. If $m<n$,

$$
P_n(H)=\frac{(n-1)_{m-1}}{a_H}\,
 \frac{b_H(n-m+b_H)^{n-m-1}}{n^{n-2}},
\qquad P_\infty(H)=\frac{b_H}{a_H}e^{-(m-b_H)} .
\tag{9}
$$

Here $(z)_j=z(z-1)\cdots(z-j+1)$. For the first identity prescribe a labelled copy, allowing outside edges only at its boundary. Contract components of the prescribed forest, including singleton vertices outside the copy. If the components have $b_i>0$ possible attachment vertices, the number of completions is $(\prod_i b_i)(\sum_i b_i)^{r-2}$ for $r\ge2$ components. Each component of degree $d$ contributes $b_i^d$, and summation over Prüfer words of length $r-2$ gives the expression. A component with no attachment vertices gives zero. This is also Lemma 2.3 of Isaev, Southwell and Zhukovskii. In an ordinary Poisson tree the ball has probability $e^{-(m-b_H)}/a_H$, by the same ordering computation restricted to interior vertices. Size bias by the population at level $h$ gives the one-spine ball law: marking one of the $b_H$ boundary vertices size biases the offspring at each ancestor and chooses the continuing child uniformly. Size-biased Poisson offspring has law $1+\operatorname{Pois}(1)$. [Isaev, Southwell and Zhukovskii, Lemma 2.3](https://doi.org/10.1017/S0963548321000523)

For $m\le n^{1/8}$ and $b_H>0$, (9) yields $P_n(H)/P_\infty(H)=1+O(m^2/n)$ with an absolute constant. In the log ratio the falling factorial contributes $\sum_{i=1}^{m-1}\log(1-i/n)=O(m^2/n)$. The other contribution is $(n-m-1)\log(1-(m-b_H)/n)+(m-b_H)=O(m^2/n)$. These estimates follow from $|\log(1-x)+x|\le x^2/(1-x)$ for $0\le x<1$, by integrating $t/(1-t)$ from zero to $x$. The ratios are uniformly close to one, so exponentiation preserves the bound. For smaller balls with $b_H=0$ both probabilities vanish. Choose $h=\lceil 3\log_2 n\rceil$. By (8) and Markov's inequality the tail $m>n^{1/8}$ has probability $O((h+1)^{32}n^{-2})$ under either law. Integrating the small-ball error costs $O((h+1)^4/n)$ by (8). In (7) the expected truncation error is $O((h+1)2^{-h})$ under either law. Since $\mu(T)=\sum_v s(T,v)$, multiplication by $n$ proves

$$
\mathbb E\mu(T_n)=cn+O((1+\log n)^4).
\tag{10}
$$


## 4. Controlling the variance lost on truncation

For any forest functional $f$ used below, set $\Delta_e f(T)=f(T)-f(T-e)$, evaluating the same vertex scores within each component after deletion. A useful estimate for label-invariant $f$ is

$$
\operatorname{Var} f(T_n)\le 2\,\mathbb E\sum_{e\in T_n}|\Delta_e f(T_n)|^2 .
\tag{11}
$$

Use the independent coordinates of the Cayley sampler in §5, stages I–II, of Isaev, Southwell and Zhukovskii. For $1\le i<n$, join $i+1$ to $\min(i,U_i)$ with $U_i$ independent uniform on $[n]$, then permute the labels uniformly. The permutation does not affect $f$. On refreshing coordinate $i$, the two resulting trees agree after removing their respective $i$th edges. For independent coordinates the variance is at most half the sum of the expected squared refresh differences. To see this, take martingale differences on revealing successive coordinates, express the $i$th difference as the conditional expectation of $f-\mathbb E[f\mid(U_j)_{j\ne i}]$ given the first $i$ coordinates, then apply Jensen and sum. Each squared refresh difference is bounded by twice the sum of the two squared cut differences. The refreshed tree has the original coordinate law for each $i$, so the two sums have equal expectation. Summing over $i$ visits every edge once and gives (11). This argument also covers refreshes that leave the edge unchanged. [Isaev, Southwell and Zhukovskii, §5, stages I and II](https://doi.org/10.1017/S0963548321000523)

For a nonnegative vertex quantity $W$, write $\|W\|_*^2=\mathbb E[\deg(V)W(V)^2]$ for $V$ uniform in $T_n$. If $|\Delta_{xy}f|\le W(x)+W(y)$, then (11) gives $\operatorname{Var}(f)/n\le4\|W\|_*^2$. Indeed the squared sum is at most $2(W(x)^2+W(y)^2)$ and each vertex is counted once per incident edge. By Cauchy–Schwarz, $\|W\|_*\le\|\deg(V)\|_2^{1/2}\|W(V)\|_4$. Degree has bounded moments by (8) at $j=1$.

Put $F_h(T)=\sum_v s_h(T,v)$. Cut an edge $xy$ and consider a vertex on the $x$ side at distance $j$ from $x$. Its full score changes by at most $2^{-j}$ because only one child message can change at each step of the path from $x$, starting with a difference at most one. Its radius-$h$ score changes by at most the same amount, and is unchanged for $j\ge h$. These comparisons hold even for a high-degree endpoint. For $j>\lfloor h/2\rfloor$ use the two change bounds. For $j\le\lfloor h/2\rfloor$ use (7) before and after the cut. Each radius-$h$ ball involved lies in $B_{2h}(x)$ of the original tree and deletion only removes vertices from balls. Thus an endpoint bound for the cut of $\mu-F_h$, for $h\ge1$, is

$$
W_h(x)=2^{1-h}|B_{2h}(x)|^2+
2\sum_{j>\lfloor h/2\rfloor}2^{-j}Z_j(x).
\tag{12}
$$

By (8), Hölder and summation of the geometric tail with a linear factor, $\|W_h\|_*=O(2^{-h}(h+1)^4+2^{-h/2}(h+1))$ uniformly in $n$. Consequently $\sup_n \operatorname{Var}(\mu-F_h)/n\to0$.

For integers $h,D\ge1$, put $F_{h,D}(T)=\sum_v s_h(T,v)\mathbf1_{\{|B_h(v)|\le D\}}$. In $F_h-F_{h,D}$ only vertices within distance $h-1$ of the cut on their side can change. If any such summand changes on the $x$ side, its ball in the original tree has size greater than $D$, and is contained in $B_{2h}(x)$. An endpoint bound is $|B_{2h}(x)|\mathbf1_{\{|B_{2h}(x)|>D\}}$, since each summand changes by at most one. At fixed $h$ its weighted norm tends to zero uniformly in $n$ since on the tail $|B|^4\le |B|^8/D^4$ and (8) applies. In particular,

$$
\lim_{h\to\infty}\limsup_{D\to\infty}\sup_{n\ge2}
\frac{\operatorname{Var}(\mu(T_n)-F_{h,D}(T_n))}{n}=0 .
\tag{13}
$$

We record an explicit upper variance bound. An endpoint bound for cutting $\mu$ is $L(x)=\sum_{j\ge0}2^{-j} Z_j(x)$. For a fixed root, two ordered length-$j$ paths from it have a common initial segment of length $\ell\in[0,j]$ and otherwise are disjoint. Their union has $m=2j-\ell+1$ vertices. For any one abstract marked union, at most $(n-1)_{m-1}$ labellings are possible. The probability that the tree contains its prescribed edges is $m/n^{m-1}$ by the completion formula with all $m$ vertices available (also valid for $m=n$ by Cayley's formula). Thus each $\ell$ contributes at most $m$ to $\mathbb E Z_j^2$. To include $\deg(V)$, mark a neighbor of the root. At most two choices are already in the union and adjacent to the root. Any other vertex of the union would create a cycle. A new neighbor gives an abstract union on $m+1$ vertices and contributes at most $m+1$. Counts requiring more than $n$ vertices vanish. For each $\ell$ this total is at most $3m+1\le6(j+1)$, including $j=0$. Hence $\mathbb E[\deg(V)Z_j^2]\le6(j+1)^2$. Minkowski's inequality gives $\|L\|_*\le\sqrt6\sum_{j\ge0}(j+1)2^{-j}=4\sqrt6$, and (11) gives $\operatorname{Var}(\mu(T_n))\le384n$.

## 5. A Gaussian law for the random center

Fix $h,D$. We verify the hypotheses of Theorem 1.1 of Isaev, Southwell and Zhukovskii for $F=F_{h,D}$. Its switch deletes $ij$ and inserts $ik$, keeping the result a tree. First consider a single deletion. A changed, nonzero contribution on the $x$ side of $xy$ must have its center within distance $h-1$ of $x$ and its radius-$h$ ball after deletion has size at most $D$ (the ball after deletion is the smaller one). Every vertex on the path from $x$ to that center has degree at most $D$ within this component, because all of its neighbors belong to the center's ball. There are at most $1+D+\cdots+D^{h-1}$ candidate paths, obtained by extending only from vertices of degree at most $D$. Each contribution changes by at most one. Sum over the two sides and two cuts through the common intermediate forest. This proves a deterministic switch-Lipschitz bound, independent of $n$.

For two valid switches $ij\mapsto ik$ and $ab\mapsto ac$, suppose $d_T(\{j,k\},\{b,c\})\ge2h+4$, using the source's definition of separation. Their endpoint triples are at distance at least $2h+2$ in $T$, as $i$ is adjacent to $j$ and $a$ to $b$. To bound the distance in the union of the four switched trees, observe that every extra edge lies inside one of the triples, and a path between triples has a segment between its last visit to the first and first subsequent visit to the second, using only old edges. For any center, at least one triple is at distance greater than $h$ in the union. That switch has no effect on its ball in either state of the other switch. Summing the resulting zero mixed differences proves superposability with constant $2h+4$. [Isaev, Southwell and Zhukovskii, definitions before Theorem 1.1](https://doi.org/10.1017/S0963548321000523)

We also need a variance *limit* for $F$. There are finitely many rooted ball shapes of radius at most $h$ and size at most $D$. Count their occurrences by injective vertex maps, dividing by the appropriate rooted automorphisms. Impose all prescribed edges and forbid further edges at interior vertices. This characterizes a ball exactly in a tree. To compute its first or second moment, enumerate identifications in one or two copies. Only finitely many types occur at fixed $h,D$. A cycle or a prescribed extra edge at a forbidden vertex gives zero. In a consistent union, each vertex may receive extra edges only if all copies containing it permit them. If there are $m$ distinct vertices in $r$ connected components, $r=1$ or $2$, and $b_i$ permitted vertices in component $i$, the term including the choice of labels is a constant times

$$
(n)_m\,\frac{\left(\prod_i b_i\right)
 (n-m+\sum_i b_i)^{n-m+r-2}}{n^{n-2}} .
\tag{14}
$$

Here $n>2D$. A zero $b_i$ gives a zero term. The other components in the completion formula are the $n-m$ unused vertices, each a singleton. The same expression applies when the two marked centers coincide. For fixed parameters (14) has expansion $n^r(A_0+A_1/n+O(n^{-2}))$ by expanding the finite product $(n)_m/n^m$ and $(n-m+r-2)\log(1-(m-\sum b_i)/n)$ through order $1/n$. The latter uses the Taylor series through its quadratic term, with remainder $O(n^{-3})$ before multiplication by $n-m+r-2$. Subtracting the square of the mean therefore gives $A n^2+B n+O(1)$ for the variance. The bounded single-cut estimate and (11) give $O(n)$, hence $A=0$. Write $v_{h,D}=B\ge0$.

If $v_{h,D}>0$, the two smallness terms in the cited Theorem 1.1, with switch bound $L$ and separation $\rho$, are $n L^3/\operatorname{Var}(F)^{3/2}$ and $n^{1/4}L\rho/\sqrt{\operatorname{Var}(F)}$. They are $O(n^{-1/2})$ and $O(n^{-1/4})$ for fixed $h,D$, as required. Thus $(F-\mathbb EF)/\sqrt n$ converges to $\mathcal N(0,v_{h,D})$. When $v_{h,D}=0$, convergence to zero follows in $L^2$. [Isaev, Southwell and Zhukovskii, Theorem 1.1](https://doi.org/10.1017/S0963548321000523)

For any $\epsilon>0$, (13) supplies fixed $h,D$ so that the $L^2$ distance between the centered, $\sqrt n$-scaled $\mu$ and $F_{h,D}$ is at most $\epsilon$ for every $n$. Their standard deviations consequently differ by at most $\epsilon$. The limsup and liminf of $\sqrt{\operatorname{Var}(\mu)/n}$ differ by at most $2\epsilon$ and are finite by §4. Denote the common limit by $\sqrt v$. At each real $t$ the characteristic functions of the two centered variables differ by at most $|t|\epsilon$, by $|e^{ita}-e^{itb}|\le |t||a-b|$. Taking $n\to\infty$ and then $\epsilon\downarrow0$ gives the characteristic function $e^{-vt^2/2}$, hence the claimed normal limit for $\mu$. Equation (10) permits centering by $cn$. Equation (5) transfers both the distributional and the variance limit to every mode, since the $\sqrt n$-scaled error tends uniformly to zero. It remains to prove $v>0$.

## 6. A positive contribution to the variance

A rooted shape is a tuple of its child shapes. Write $L=()$, $P=(L)$ for a root with one leaf child, $Q=(P)$, and $Y=(L,L)$. Consider these ten-vertex shapes.

$$
A=(L,Y,(L,Q)),\qquad B=(P,P,(P,P)).
$$

Let $Z_0,Z_1$ be the numbers of sets excluding and including the root respectively. The corresponding sums of their sizes are $D_0,D_1$. The table follows by enumerating all $2^{10}$ subsets.

| Shape | $Z_0$ | $Z_1$ | $D_0$ | $D_1$ | Root-preserving automorphisms |
|---|---:|---:|---:|---:|---:|
| $A$ | 130 | 40 | 395 | 140 | 2 |
| $B$ | 117 | 36 | 336 | 120 | 4 |

The vectors $(Z_0,Z_1)$ are proportional. Replacing an $A$ by a $B$ attached at the same root preserves the joint law of the outside configuration and root state, since both state weights are multiplied by $9/10$. In *each* state the conditional mean size inside decreases by $1/6$, as one sees by subtracting $D_0/Z_0$ or $D_1/Z_1$. Thus the full mean decreases by exactly $1/6$ at any such replacement, independently of the other branches.

Root $T_n$ at label 1 and mark every descendant subtree of type $A$ or $B$. All marks are disjoint because the shapes have equal size. Reveal the outside tree, all attachment roots and their label blocks, retaining only the fact that each marked block has one of the two shapes. Any substitution preserves the marks because ancestors retain their sizes, and proper descendants of the blocks have size less than ten. Conditional on the revealed information, choices across blocks are independent. Indeed all completions are equiprobable and on each block of nine internal labels and one specified root the two multiplicities are $9!/a_A$ and $9!/a_B$, with $a_A=2$, $a_B=4$. With $R$ marks and $p_A=(1/a_A)/(1/a_A+1/a_B)$ the conditional variance of $\mu$ is $R p_A(1-p_A)/36$.

For $n\ge12$, the expected number of marks of shape $H$ divided by $n$ is

$$
\frac{(n-1)(n-2)_9(n-10)^{n-11}}{n\,a_H n^{\,n-2}}
\ \longrightarrow\ \frac{e^{-10}}{a_H}.
\tag{15}
$$

Choose the marked root outside label 1, the nine other labels outside both roots, a tree on the $n-10$ outside labels, and the outside attachment vertex. The count of outside trees is $(n-10)^{n-12}$ by Prüfer coding. Dividing by the total number of Cayley trees gives (15). The law of total variance yields $v\ge e^{-10}/(36(a_A+a_B))=e^{-10}/216$. Since $j!\ge2^{j-1}$ for $j\ge2$, strictly for $j\ge4$, the exponential series gives $e<3$ and the rational lower bound in (2). This completes the proof of Theorem 1.

## 7. A finite certificate for the limiting fraction

The calculation integrates distributions of messages in (6). For a height-$H$ truncation put message zero at depth $H$. At an ordinary vertex the incoming messages come from $N\sim\operatorname{Pois}(1)$ independent ordinary children. At a spine vertex there is one additional independent spine child. All Poisson variables and child subtrees in separate branches are independent. Update the two laws in parallel from the preceding height, and evaluate the expectation of $b(x)$ at the spine root. By (7) and $\mathbb E Z_H=H+1$ the boundary error is at most $(H+1)2^{-H}$.

**Interpolation error.** On a grid of mesh $1/J$ replace a value $y$ by one of its two neighboring grid points, with weights giving mean $y$. The computation propagates both weights exactly in the ideal grid model. For one replacement the change in expectation of a twice differentiable test function $\psi$ is bounded by $\|\psi''\|_\infty/(8J^2)$ by Taylor's formula about $y$, cancellation of the linear term, and variance at most $1/(4J^2)$. When all ancestors of a node at depth $d$ are still unrounded, the map from its message to root occupation has both derivative bounds at most $2^{-d}$. Start with $b$. Composing with $g(\cdot+s)$ for a fixed nonnegative sibling sum multiplies the first-derivative bound by at most $1/2$ and bounds the second by $1/4$ times the sum of the two preceding bounds. Apply the replacements starting at the deepest nodes. Conditional on their input and the other branches, each replacement thus has this expectation bound even if descendants have already been rounded. Summing over the expected $d+1$ vertices in generation $d$ bounds the total grid error by $1/(2J^2)$, since $\sum_{d\ge0}(d+1)2^{-d}=4$.

**Arithmetic and offspring tails.** Set $J=1024$, $H=40$, $K_0=17$, $S=2^{60}$. Use integer probability masses summing to $S$. For incoming grid sum $i$, the outgoing grid coordinate is $y_i=J\log_2(1+2^{-i/J})$. The verifier encloses these coordinates up to $i=(K_0+1)J$. It brackets $2^{-1/J}$ by integer bisection at precision $Q=2^{112}$, encloses its powers by directed multiplication, and uses

$$
\log(1+x)=2\sum_{r=0}^{R-1}\frac{w^{2r+1}}{2r+1}+E_R,\qquad
w=\frac{x}{2+x},\quad
0\le E_R\le\frac{2(1/3)^{2R+1}}{(2R+1)(1-1/9)}
$$

for $0\le x\le1$ and $R=44$. This follows by integrating the geometric series for $2/(1-w^2)$ from zero to $w$. The program checks that each grid coordinate differs from its stored value by at most $2/S$. The total variation distance between two adjacent-point splitting laws is at most the difference between their real grid coordinates. Moving the coordinate inside one grid cell transfers exactly that much mass between the two endpoints. Sum over intervening cells.

Use $N$ conditioned on $N\le K_0$, with exact integer weights proportional to $K_0!/d!$. The total variation cost is $\Pr(N>K_0)\le 2/(K_0+1)!$ by the factorial series and a geometric tail bound. Each convolution multiplies integer polynomials exactly (base $2^{128}$ exceeds the largest raw coefficient $S^2$). Restoring denominator $S$ floors the masses and returns the remainder to index zero. On a support of size $l$ this changes total variation by at most $l/S$, as the total discarded mass is less than $l/S$. The same estimate holds for the mixture normalization and the final split normalization. Convolution with a probability law and applying a Markov kernel do not increase total variation, by the triangle inequality after summing absolute differences. Thus the successive convolutions forming the powers up to $K_0$ accumulate at most $K_0$ such errors in any mixture term. The mixture, the extra spine convolution, the final split and the coordinate discrepancy together give the following one-update bound for either output, at fixed input laws.

$$
\epsilon=\frac{2}{(K_0+1)!}
 +\frac{(K_0+4)((K_0+1)J+2)}{S}.
\tag{16}
$$

To propagate errors compare first to an *exact Poisson* update on the computed input laws. If input ordinary and spine laws have total variation errors $a,b$, independent couplings of their children give output errors at most $a+\epsilon$ and $a+b+\epsilon$, since $\mathbb E N=1$. Iteration from zero errors bounds the spine error by $H(H+1)\epsilon/2$. Occupation lies in $[0,1]$, so the same error bounds its expectation. At the last step the powers of $2^{-1/J}$ already enclosed give outward rational bounds for the root average. Enlarging them on both sides by

$$
(H+1)2^{-H}+\frac1{2J^2}+\frac{H(H+1)}2\,\epsilon
$$

and rounding outward yields (2). For Corollary 2, the odd and even Taylor partial sums of $e^{-x}$ enclose it for $0\le x\le1$ because the terms decrease. Testing $x=0.56714329$ and $x=0.56714330$ with orders 17 and 18 respectively proves these endpoints bracket the unique root of $e^{-x}-x$. Division then gives the stated enclosure for $c/\omega$.

**Reproduction.** Run `python3 verify.py` in this directory. It computes the rational density certificate, verifies the two ten-vertex shapes by direct subset enumeration, compares the polynomial recurrence with exhaustive subset counts, and checks the exact laws of radius-one and radius-two balls for all 1,441 labelled trees of orders 2 through 6. As a separate, non-rigorous check, 48 seeded uniform Prüfer trees at each of $n=64,128,256$ gave mean mode fractions $0.30664,0.30762,0.30640$. The analytic arguments in §§2–6 establish the distributional limit. The remaining quantitative questions and a separate audit prompt accompany the manuscript. [Proof boundaries](GAPS.md) · [Independent review prompt](independent-review-prompt.md)

Work began at repository commit `528e7f6e8238ae59a172b2fb7da34c2e8762c793` and was integrated on top of `f5e976d`. OpenAI Codex prepared this research preprint at Jerald Hu Yuan’s request. Independent review of the tree-to-tree fluctuation proof is pending.

## References

Basit, A., and D. Galvin. *On the independent set sequence of a tree*. Electronic Journal of Combinatorics 28(3), P3.23. [Preprint](https://arxiv.org/abs/2006.12562)

Dembo, A., A. Montanari and N. Sun. *Factor models on locally tree-like graphs*. [Paper](https://arxiv.org/abs/1110.4821)

Devroye, L., M. K. Goh and R. Y. Zhao. *On the peel number and the leaf-height of Galton–Watson trees*. [Published paper](https://www.cambridge.org/core/journals/combinatorics-probability-and-computing/article/abs/on-the-peel-number-and-the-leafheight-of-galtonwatson-trees/631190C0A22C46DD848D98D26A5AC5CD)

Fang, E. X., J. Lu, E. Nevo, Y. Yao and H. Zheng. *Unimodality of independence polynomials for sufficiently large forests*. [arXiv:2609.20961](https://arxiv.org/abs/2609.20961)

Heilman, S. *Independent sets of random trees and sparse random graphs*. Journal of Graph Theory 109. [Published paper](https://doi.org/10.1002/jgt.23225)

Isaev, M., A. Southwell and M. Zhukovskii. *Distribution of tree parameters by martingale approach*. [Published paper](https://doi.org/10.1017/S0963548321000523)

Janson, S. *Simply generated trees, conditioned Galton–Watson trees, random allocations and condensation*. [Paper](https://arxiv.org/abs/1112.0510)

Li, W., K. Vallier and T. Zhang. *Unimodality of forest independence polynomials*. [arXiv:2610.07943](https://arxiv.org/abs/2610.07943)

Ralaivaosaona, D., M. Šileikis and S. Wagner. *A central limit theorem for almost local additive tree functionals*. [Paper](https://arxiv.org/abs/1810.00467)
