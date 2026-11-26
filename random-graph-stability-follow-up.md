# How often is a random graph's independence polynomial unstable?

Research draft · 26 November 2026 · Developed with OpenAI Codex · Unrefereed

Let $G\sim G(n,p)$, with independent edges of fixed probability $p$. Write $X_s$ for the number of independent $s$-sets, $I_G(z)=\sum_s X_sz^s$, and

$$
u_n(p)=\Pr\{I_G(z)=0\text{ for some }z\text{ with }\Re z>0\}.
$$

This is the stability convention in Brown and Cameron's Problem 2. At $p=1/2$ the distribution is uniform over all labelled graphs of the given order. An earlier note gave $\limsup_n u_n(1/2)\ge 1-e^{-1/4}$, using the second Hurwitz determinant. Here we determine the entire set of subsequential limits. [Original question](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v25i1p46/pdf/) · [Earlier note](https://github.com/jeraldhu-yuan/random-graph-independence-stability/blob/e49e2d90d2a43f936eb5714a6208967255799ff3/random-graph-stability.md)

**Theorem 1 (uniform labelled graphs).** The set of subsequential limits of $u_n(1/2)$ is the closed interval $[0,L]$, where

$$
0.27149223459<L=1-\exp(-1/\tau)<0.27149223460.
$$

The constant $\tau$ is defined by the convergent series in Section 3. Moreover,

$$
\frac1N\sum_{n=1}^N u_n(1/2)=\Theta\bigl(1/\log N\bigr).
$$

In particular, for every fixed $\varepsilon>0$, the proportion of orders $n\le N$ with instability probability above $\varepsilon$ is $O_\varepsilon(1/\log N)$. Choosing an order uniformly from $1,\dots,N$ and then drawing a uniform graph of that order gives stability with probability tending to one.

**Theorem 2 (smaller edge probabilities).** Let $s_*=0.3092493386\ldots$ be the first positive spectral value of the partial theta function. For every fixed

$$
0<p<1-s_*^{1/4}=0.254277589\ldots,\qquad \lim_{n\to\infty}u_n(p)=1.
$$

Section 6 gives an elementary proof on the subrange $p<1-2^{-1/4}$, which includes $p=0.1$. Section 5 also gives a sufficient condition at larger $p$ for the mean instability probability to decay as $1/\log N$. The behavior between these ranges remains to be classified.

The reason for the oscillation at $p=1/2$ is that the largest independent sets appear at widely spaced graph orders. Well below or above one of these orders, the coefficients force stability. During the short interval when a new size first appears, its count has a Poisson limit. A pair of roots crosses the imaginary axis when that count exceeds a computable multiple of its mean. Brown and Nowakowski previously proved that a uniform graph almost surely has a nonreal independence root as $n$ grows. The sign of its real part is the issue here. [Brown and Nowakowski, Proposition 1](https://doi.org/10.1016/j.jctb.2004.10.001)

## 1. A bound on the imaginary axis

Let $B(z)=\sum_{j=l}^d b_j z^j$, with consecutive positive coefficients. Suppose

$$
\frac{b_j^2}{b_{j-1}b_{j+1}}\ge \rho>1\qquad(l<j<d),
\qquad h(\rho)=1-\rho^{-1}-\rho^{-3}>0.
$$

For $y>0$ we have

$$
|B(iy)|\ \ge\ h(\rho)\max_j b_j y^j. \tag{1}
$$

To see this, put $T_j=b_jy^j$ and choose an index $m$ maximizing $T_j$. Let $v=T_{m-1}/T_m$ and $w=T_{m+1}/T_m$, using zero for missing neighbors. Then $0\le v,w\le1$ and $vw\le1/\rho$. In the real part of $i^{-m}B(iy)$, the terms on each side of $T_m$ alternate and decrease in magnitude. This follows from the decreasing consecutive ratios $T_{j+1}/T_j$. Each side is bounded below by its first negative term. Their magnitudes are at most $T_m v^2/\rho$ and $T_m w^2/\rho$. Finally $v^2+w^2\le1+v^2w^2\le1+\rho^{-2}$. This proves (1), including a maximum at an endpoint.

When $l=0$, all the roots lie in the **open** left half-plane. For completeness, multiply $b_j$ by $t^{\binom j2}$, $0<t\le1$. All the coefficient ratios in the hypothesis increase by a factor $1/t$, so (1) excludes imaginary-axis roots throughout the deformation. The constant term excludes zero. At $z=-t^{1/2-m}$ the power of $t$ in term $j$ is $(j-m)^2/2-m^2/2$. For sufficiently small $t$, term $m$ determines the sign for each $0\le m\le d$. The $d+1$ signs alternate in order on the negative real axis, giving $d$ negative real roots. The degree stays fixed and the roots vary continuously for $t$ between this positive value and one. The claim follows. Sufficient conditions for strict Hurwitz stability in terms of consecutive coefficient ratios are classical. [Lipatov and Sokolov](https://www.mathnet.ru/eng/at9834)

## 2. Uniform control of the coefficients below the random leading term

Put $q=1-p$, $b(j)=\binom j2$, $\mu_s=\mathbb E X_s$, and

$$
\mu_s=\binom ns q^{b(s)},\quad r_s=\frac{\mu_s}{\mu_{s-1}}=\frac{n-s+1}{s}q^{s-1},
\quad a_s=(s!)^{1/s}q^{-(s-1)/2}\quad(s\ge1).
$$

The symbol $a_s$ denotes a graph order scale. Stirling's formula gives $a_s\sim (s/e)q^{-(s-1)/2}$ and $a_{s+1}/a_s\to q^{-1/2}>1$. Choose $j$ with $a_j\le n<a_{j+1}$ for large $n$. Then $j\asymp\log n$, $n\asymp j q^{-j/2}$, and $r_{j+i}\asymp j/n$ for fixed $i$. Uniformly for $s=j+O(1)$,

$$
\mu_s=(n/a_s)^s\prod_{i=0}^{s-1}(1-i/n)
      =(n/a_s)^s(1+O(j^2/n)). \tag{2}
$$

In particular $\mu_j\ge1-o(1)$ and $\mu_{j+1}\le1$. Take $k=j$ if $\mu_j\le n^{1/3}$ and $k=j+1$ otherwise. In both cases $\mu_k\le n^{1/3}$, $\mu_{k-1}\ge n^{1/3}$ for all sufficiently large $n$, and $\mu_{k+1}=O(kn^{-2/3})$. Markov's inequality gives $X_{k+1}=0$ with probability tending to one. Any larger independent set would contain one of size $k+1$.

Write $V_s=\operatorname{Var}(X_s)/\mu_s^2$. For independent uniformly chosen $s$-subsets $S,T$, let $J_s=|S\cap T|$. Pair counting gives

$$
V_s=\mathbb E[q^{-b(J_s)}]-1.
$$

In fact, an ordered pair with overlap $i$ needs $2b(s)-b(i)$ absent edges. Coupling the two sets as the first $s$ vertices in two independent random permutations makes $J_s$ nondecreasing in $s$. Consequently $V_s$ is nondecreasing as well.

For $s=j+O(1)$ we need the uniform bound

$$
V_s=O(j^3/n+1/\mu_s). \tag{3}
$$

Here are the details. Put $h_i=\Pr(J_s=i)q^{-b(i)}$. The contributions at $i=0,1$ to $V_s$ vanish. For $2\le i\le s/2$, a union bound over the prescribed $i$ vertices gives $\Pr(J_s=i)\le(s^2/n)^i$. Thus $h_i\le A^i$ with $A=(s^2/n)q^{-s/4}=O(j^{3/2}/\sqrt n)$. For $i=s-t>s/2$, exact pair counting gives

$$
h_{s-t}=\mu_s^{-1}\binom st\binom{n-s}t q^{t(s-(t+1)/2)}
\le \mu_s^{-1}D^t,\qquad D=snq^{3s/4-1/2}=O(j^{5/2}/\sqrt n).
$$

This includes $t=0$. The two geometric sums give $V_s\le A^2/(1-A)+1/[\mu_s(1-D)]$, proving (3). The underlying overlap count is the classical one for cliques in the complement graph. [Bollobás and Erdős, “Cliques in random graphs”](https://doi.org/10.1017/S0305004100053056)

By monotonicity, Chebyshev's inequality and a union bound, for each $\delta>0$ the probability that any $s\le k-1$ violates $|X_s/\mu_s-1|\le\delta$ is at most $k V_{k-1}/\delta^2=o(1)$. Thus

$$
\max_{0\le s<k}|X_s/\mu_s-1|\ \longrightarrow 0
\quad\text{in probability}. \tag{4}
$$

Define $c_j=r_k^j X_{k-j}/\mu_k$ for $1\le j\le k$ and $c_j=0$ for $j>k$. At each fixed $j\ge1$, $c_j\to q^{b(j)}$ in probability. Indeed

$$
\frac{r_k^j\mu_{k-j}}{\mu_k}
=\prod_{\ell=0}^{j-1}\frac{r_k}{r_{k-\ell}}
\le q^{b(j)},
$$

and at fixed $j$ its ratio to the right side tends to one. Together with (4) this bounds every $c_j$ by $(1+o_{\Pr}(1))q^{b(j)}$, with one common error. Therefore $B_n(z)=\sum_{j\ge1}c_jz^j$ converges in probability, uniformly on compact sets together with its derivatives, to $B(z)=\sum_{j\ge1}q^{b(j)}z^j$. This follows by bounding the tails by the convergent Gaussian-decay series on each disk. Also, for $2\le j\le k-1$, the ratios $c_j^2/(c_{j-1}c_{j+1})$ are simultaneously at least $(1-o_{\Pr}(1))/q$. The ratios of the corresponding means are $r_s/r_{s+1}\ge1/q$.

On $X_{k+1}=0$, set $t=X_k/\mu_k$. We have the exact identity

$$
F_{n,t}(z):=t+B_n(z)=(r_kz)^k I_G(1/(r_kz))/\mu_k. \tag{5}
$$

The reciprocal map with positive scale preserves the sign of the real part of each nonzero root. If $X_k=0$, padding to degree $k$ contributes a zero at the origin and preserves this equivalence for instability.

## 3. The first imaginary-axis crossing for the uniform graph

Set $q=1/2$ throughout this section. The imaginary part of $B_n(iy)$ is $yO_n(y^2)$ where $O_n(v)=\sum_{m\ge0}(-1)^m c_{2m+1}v^m$. Its deterministic limit and the negative of the real part at $y=\sqrt v$ are

$$
O(v)=\sum_{m\ge0}(-1)^m2^{-b(2m+1)}v^m,\qquad
T(v)=\sum_{m\ge1}(-1)^{m+1}2^{-b(2m)}v^m .
$$

The function $O$ has one zero $\xi$ in $[0,81]$. On $[0,9]$ the alternating derivative satisfies $O'(v)\le-1/8+2v/1024<0$. On $[9,81]$ the alternating tail gives $O(v)\le 1-v/8+v^2/1024<0$, since this convex quadratic is negative at both endpoints. The ratio of successive absolute terms is $v/2^{4m+3}$ from index $m$ to $m+1$, which justifies the alternating estimates in their stated ranges. For the derivative, include the factor $(m+1)/m$ for $m\ge1$. The first four terms at $8.5$ and first three at $8.6$ give $O(8.5)>0>O(8.6)$. Hence $8.5<\xi<8.6$. Define $\tau=T(\xi)$. The first three terms of $T$ on this interval give $3<\tau<16/5$. Rational enclosures obtained by bisection and alternating remainders are

$$
8.57160089651<\xi<8.57160089652,\qquad
3.15699422295<\tau<3.15699422296.
$$

By compact convergence of $O_n$ and its derivative, with probability tending to one $O_n$ has a unique zero $\xi_n$ in $[0,81]$, lying in $(8.5,8.6)$, with $O'_n(\xi_n)<0$. On these events $\xi_n\to\xi$ and $\tau_n=-\Re B_n(i\sqrt{\xi_n})\to\tau$ in probability. Variables defined on such events may be assigned arbitrary values on their complements.

**Claim.** With probability tending to one, for *every* real $t\ge0$, the polynomial $F_{n,t}$ has a root with positive real part exactly when $t>\tau_n$.

On an event of probability tending to one we have $\tau_n<16/5$, $c_1>99/100$, all $c_1,\ldots,c_k>0$, and all their interior coefficient ratios exceed $\rho=99/50$. Equation (1) gives, uniformly for $|y|\ge9$,

$$
|B_n(iy)|\ge h(\rho)c_1|y|>9(99/100)h(99/50)>16/5.
$$

Thus an imaginary-axis zero for $0<t\le16/5$ must satisfy $|y|<9$. It cannot have $y=0$. For $0<y<9$ it requires $y=\sqrt{\xi_n}$ and $t=\tau_n$.

At $t=0$, $B_n(z)/z$ is strictly stable by Section 1, and the zero of $B_n$ at the origin is simple. For sufficiently small positive $t$, that zero moves to $-t/c_1+O(t^2)$ and the remaining roots stay in the open left half-plane. Continuity and the absence of a crossing give strict stability throughout $0<t<\tau_n$, and closed left-half-plane stability at $t=\tau_n$.

For any $t>\tau_n$, at $y=\sqrt{\xi_n}$ the value $F_{n,t}(iy)=t-\tau_n$ is positive real. The derivative of its imaginary part with respect to $y$ is $2\xi_n O'_n(\xi_n)<0$, so the derivative of its argument is negative there. For any polynomial $F$ whose roots $\zeta$ all have $\Re\zeta\le0$, at any $F(iy)\ne0$ one has

$$
\frac{d}{dy}\arg F(iy)
=\Re\frac{F'(iy)}{F(iy)}
=\sum_\zeta\frac{-\Re\zeta}{|iy-\zeta|^2}\ge0.
$$

This contradiction proves the claim for all $t>\tau_n$. In particular it covers every Hurwitz determinant, even when a lower-order determinant stays positive.

## 4. All subsequential limits

The preceding claim and (5) show that, up to an event of probability $o(1)$, instability is the event $X_k>\mu_k\tau_n$. If $\mu_k\to0$, Markov's inequality gives $u_n(1/2)\to0$. If $\mu_k\to\infty$, (3) gives $X_k/\mu_k\to1$ in probability. Since $\tau>3$, this also gives $u_n(1/2)\to0$.

If $\mu_k\to\lambda\in(0,\infty)$, the clique-count Poisson theorem applied to the complement graph gives $X_k\Rightarrow Y_\lambda\sim\operatorname{Poisson}(\lambda)$. Its hypotheses hold along this sequence of orders with fixed edge probability and mean tending to $\lambda$. Since $\mu_k\tau_n\to\lambda\tau$ in probability, the open and closed set bounds for convergence in distribution give

$$
\Pr(Y_\lambda>\lambda\tau)\le\liminf u_n(1/2)
\le\limsup u_n(1/2)\le\Pr(Y_\lambda\ge\lambda\tau). \tag{6}
$$

In particular the two bounds agree whenever $\lambda\tau$ lies strictly between consecutive integers. Convergence of the random threshold to a constant justifies this argument even when $\tau_n$ and $X_k$ are dependent. [Bollobás, *Random Graphs*, 2nd ed., Theorem 11.7](https://doi.org/10.1017/CBO9780511814068.013)

Here is a direct check of the Poisson step using Section 2. For fixed $0<x<1$, independently mark each of the $M=\binom nk$ candidate sets with probability $1-x$. Let $A_i$ be the event that candidate $i$ is marked and independent, with probability $\pi=(1-x)q^{b(k)}$. The probability that none occurs is $\mathbb E[x^{X_k}]$. Call distinct candidates neighbors if they share at least two vertices. Each $A_i$ is independent of the joint collection of events at its nonneighbors, since they use disjoint edge variables from candidate $i$ and independent marks.

Order the candidates and put $v_i=\Pr(\text{none of }A_1,\ldots,A_i)$, with $v_0=1$. In both $\Pr(A_i,\text{ no earlier event})$ and $\pi v_{i-1}$, discard the restrictions on earlier neighbors. The resulting expressions agree by independence, so the difference in absolute value is at most $\sum_{h}[\Pr(A_i\cap A_h)+\pi^2]$, summed over earlier neighbors $h$. This bounds $|v_i-(1-\pi)v_{i-1}|$. Summing over $i$ and using $\Pr(A_i\cap A_h)\ge\pi^2$ bounds $|v_M-(1-\pi)^M|$ by the ordered sum of joint probabilities of distinct neighboring candidates. In Section 2 notation this is

$$
(1-x)^2\mu_k^2\sum_{i=2}^{k-1} h_i
\ \le\ (1-x)^2\mu_k^2\left(\frac{A^2}{1-A}+\frac{D}{\mu_k(1-D)}\right)=o(1).
$$

Since $M\pi\to(1-x)\lambda$ and $\pi\to0$, the generating functions converge to $\exp[-(1-x)\lambda]$. Bounded means give tightness and the limiting generating function uniquely identifies the Poisson law.

Put $m=\lceil\lambda\tau\rceil$. If $m=1$, the upper bound is $1-e^{-\lambda}\le L$. If $m\ge2$, the second factorial moment gives

$$
\Pr(Y_\lambda\ge m)\le\frac{\lambda^2}{m(m-1)}
\le\frac{2}{\tau^2}<1-e^{-1/\tau}.
$$

The strict inequality follows from $\tau>3$ and $1-e^{-x}\ge x-x^2/2$. Every sequence has a further subsequence on which $\mu_k$ converges in $[0,\infty]$, so these cases give the global upper bound $L$.

Conversely, for any fixed $0<\lambda<1/\tau$, take $n$ to be the nearest integer to $a_d\lambda^{1/d}$ as $d\to\infty$. Rounding changes the $d$th power by a factor tending to one, since $d/n\to0$. Equation (2) gives $\mu_d\to\lambda$, and the adjacent means satisfy $\mu_{d-1}\asymp n/d$ and $\mu_{d+1}\to0$. The rule in Section 2 selects $k=d$ eventually because $j$ can only be $d-1$ or $d$ and in each case its mean selects $d$. Equation (6) gives the limit $1-e^{-\lambda}$. These limits fill $(0,L)$. The set of subsequential limits of a bounded sequence is closed, which supplies both endpoints and proves the first assertion of Theorem 1.

The difference from the four-coefficient test can be quantified. For $t>0$, the leading coefficients of the reciprocal of $F_{n,t}$, read in descending order, have limits $t,1,1/2,1/8,1/64,1/1024$. The second Hurwitz determinant is $1/2-t/8$. The third is $3/64-15t/1024$. They turn negative at $t>4$ and $t>16/5$, respectively. The full polynomial crosses at $t=\tau=3.156994\ldots$. These are thresholds in the asymptotic coefficient model.

## 5. Averaging over graph orders

Retain $q=1/2$. On a block $a_j\le n<a_{j+1}$, if $X_{j+1}=0$ and every coefficient through $j$ differs from its mean by at most the relative error $\delta=1/16$, the degree is $j$ and every interior coefficient ratio is at least $2(15/17)^2>3/2$. Section 1 applies since $h(3/2)>0$. Monotonicity of the normalized variances therefore gives

$$
u_n(1/2)\le \mu_{j+1}+\delta^{-2}(V_j+jV_{j-1})
\le C\bigl(\mu_{j+1}+\mu_j^{-1}+j^4/n\bigr). \tag{7}
$$

Here $\mu_j\ge1-o(1)$, $\mu_{j-1}\gtrsim n/j$, and (3) applies at both indices. On this block (2) bounds $\mu_{j+1}$ by $(n/a_{j+1})^{j+1}$ and $\mu_j^{-1}$ by $2(a_j/n)^j$. For an explicit sequence of orders with $u_n(1/2)\to0$, take $n$ nearest to $\sqrt{a_ja_{j+1}}$. Both powers then tend to zero. Summing the first power over integers in the entire block costs at most $1+a_{j+1}/(j+2)$, and summing the second power without the factor 2 costs at most $1+a_j/(j-1)$, by comparison with integrals. For $a_J\le N<a_{J+1}$ the block bounds sum to $O(a_{J+1}/J)=O(N/\log N)$ since $a_{j+1}/a_j\to\sqrt2$. The $j^4/n$ terms contribute $O((\log N)^5)$.

For a matching lower bound, consider integers in

$$
a_J(1/10)^{1/J}\le n\le a_J(1/5)^{1/J}.
$$

There are $\asymp a_J/J\asymp N/\log N$ such integers, all at most $N$. Their $\mu_J$ lie between $1/10-o(1)$ and $1/5+o(1)$, and $k=J$ eventually by the same adjacent-mean argument. The probabilities $u_n(1/2)$ are bounded below by a positive constant uniformly there. Otherwise one could extract a subsequence with $\mu_J\to\lambda\in[1/10,1/5]$ violating (6), which gives the positive limit $1-e^{-\lambda}$ throughout this interval. This proves the average-order assertion.

**Corollary for other edge probabilities.** For every fixed $0<p<1$ with $q=1-p$ satisfying $1-q-q^3>0$, the mean $N^{-1}\sum_{n\le N}u_n(p)$ is $\Theta(1/\log N)$. In particular the liminf is zero. For the upper bound choose a fixed $\delta>0$ small enough that $\rho=q^{-1}[(1-\delta)/(1+\delta)]^2$ satisfies $h(\rho)>0$. Equation (7) and the integral estimates then apply with constants depending on $p$. For the lower bound choose an interval of means $[\lambda_0,\lambda_1]\subset(0,q^2)$ with $\lambda_0<\lambda_1$. In the corresponding window $a_J\lambda_0^{1/J}\le n\le a_J\lambda_1^{1/J}$, $k=J$ eventually. On $X_k\ge1$, $t\ge 1/\mu_k$, and $c_1c_2-tc_3$ is negative with probability $1-o(1)$ apart from the possible failure of the simultaneous coefficient estimates. Indeed $c_1c_2/c_3\to q^{-2}<1/\lambda_1$. This second Hurwitz determinant of the reciprocal of $F_{n,t}$ excludes even closed left-half-plane stability, since the determinants of a closed stable polynomial are nonnegative by shifting its roots slightly to the left. Along every subsequence with $\mu_k\to\lambda$ in the given interval, $\Pr(X_k\ge1)\to1-e^{-\lambda}$ by the same Poisson calculation. This gives a uniform positive lower bound throughout a window of order $N/\log N$.

## 6. A regime with instability at every large order

Return to any fixed $0<q<1$. By Section 2 the odd part of $F_{n,t}$, viewed as a polynomial in $w=z^2$ after removing its factor $z$, converges uniformly on compact sets in probability to

$$
H_q(w)=\sum_{m\ge0}q^{b(2m+1)}w^m
      =\theta(q^4,w/q),\qquad
\theta(s,x)=\sum_{m\ge0}s^{m(m+1)/2}x^m .
$$

This odd part is independent of the choice of the constant term $t$ as a polynomial identity. By the weak Hermite–Biehler theorem, if $F_{n,t}$ has all its roots in the closed left half-plane, this odd-part polynomial has only real nonpositive roots. Kostov's Theorem 1(3) implies that $\theta(s,\cdot)$ has a nonreal zero for every $s\in(s_*,1)$. For fixed $q^4>s_*$, choose a closed disk around such a nonreal zero of $H_q$, disjoint from the real line and with no zero on its boundary. Uniform convergence on that boundary and Rouché's theorem put a nonreal zero of the random odd part in the disk with probability tending to one. Equivalence (5) proves Theorem 2. [Brown and Cameron, Theorem 2.1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v25i1p46/pdf/) · [Kostov, “A property of a partial theta function”](https://arxiv.org/abs/1504.01524)

There is a short independent check for $q>2^{-1/4}$. With high probability the first three odd-part coefficients are positive. If a polynomial with positive constant and nonnegative coefficients has only real nonpositive roots, its first three coefficients $A_0,A_1,A_2$ satisfy $A_1^2\ge2A_0A_2$, by expanding the product of its linear factors. Here $c_3^2/(c_1c_5)\to q^{-4}<2$. This excludes stability with probability tending to one.

## 7. Finite checks and the simulation at 86 vertices

- For a fixed graph on 86 vertices, the second Hurwitz determinant is $4530$ and the third is $-5\,819\,319$. The graph is specified by the SHAKE256 seed `BC:86:6`, with the edge encoding in the earlier note. Thus an unstable graph can pass the original four-coefficient test at precisely the order mentioned in the question.
- At $n=65\,600\,000$, the probability of instability is less than $0.009215$. Set $d=44$, $\delta=1/16$, and use $\mu_{d+1}+256\sum_{s=2}^{d} V_s$ as the upper bound, by the same event as in (7). The verifier evaluates the exact overlap sums as rational numbers.
- For the degree-16 truncation of $t+B(z)$ the exact Hurwitz test gives strict stability at $t=3.156994$. At $t=3.156995$ the first negative determinant is the sixth. This is a numerical cross-check using exact determinants of two finite polynomials. Section 3 controls all remaining degrees in the limit.

In a reproducible sample of 8,192 SHAKE256-generated graphs on 86 vertices (seeds 0 through 8191), exact Hurwitz determinants classified 123 graphs as unstable. The first negative determinant had order two in 14 cases, order three in 94, and order four in 15. Thus the four-coefficient test detected 14 of the 123 unstable graphs, or 11.38%. All 8,069 remaining graphs passed the strict Hurwitz test. A separate count of cliques in the complement and an exact rational Routh table reproduced all 8,192 classifications. Each unstable example had exactly two roots in the open right half-plane. No zero pivot occurred. The limiting probability in Theorem 1 follows from Sections 1–4.

Run `python3 -B verify_stability_followup.py` beside the original graph counter. The checker encloses $\xi,\tau,L$ by rational intervals, reproduces the finite probability bound, and counts the explicit graph in two ways. Add `--sample-86 8192` to repeat the simulation. [Follow-up checker](verify_stability_followup.py) · [Original graph counter](https://github.com/jeraldhu-yuan/random-graph-independence-stability/blob/e49e2d90d2a43f936eb5714a6208967255799ff3/verify_random_stability.py)
