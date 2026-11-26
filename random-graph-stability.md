# Instability of independence polynomials of random graphs

Research note · 26 November 2026 · Unrefereed draft

Brown and Cameron asked whether almost all graphs are stable [1, Problem 2]. In their model $G(n,p)$, every pair of distinct labelled vertices forms an edge independently with fixed probability $p$. Write $X_s$ for the number of independent sets of $s$ vertices. The independence polynomial is $I_G(z)=\sum_{s\ge0} X_s z^s$. A graph is stable if all zeros of this polynomial have nonpositive real part.

**Theorem.** Fix $0<p<1$ and put $q=1-p$. Let $U_n$ be the event that $I_G$ has a zero with strictly positive real part. Then

$$
\limsup_{n\to\infty} \Pr(U_n) \ge 1-\exp(-q^2).
$$

For $p=1/2$, every labelled graph on $n$ vertices is equally likely. For arbitrarily large $n$, more than 22% of these graphs have a zero in the open right half-plane. This gives a negative answer to Problem 2. The bound concerns a subsequence of graph orders. Whether the instability probability has a limit is a separate question.

Brown and Nowakowski showed that the average independence polynomial over all labelled graphs of a fixed order has only real, simple roots [2]. The proof here uses fluctuations in the highest coefficient of the individual graph's polynomial. Its probabilistic ingredients come from the classical theory of clique counts in the complement graph [3]. We give the estimates explicitly.

## A test using four coefficients

**Lemma 1.** Let $P(z)=a_d z^d+\cdots+a_0$ be real, with $a_d>0$ and $d\ge3$. If all its roots have nonpositive real part, then

$$
a_{d-1}a_{d-2}-a_d a_{d-3}\ge0.
$$

**Proof.** Divide by the leading coefficient. For a monic polynomial $f$ with first three coefficients below its leading term $u_1,u_2,u_3$, define $H(f)=u_1u_2-u_3$. Use zero for missing coefficients. If $g$ has corresponding coefficients $v_1,v_2,v_3$, multiplication gives

$$
H(fg)=H(f)+H(g)+u_1v_1(u_1+v_1).
$$

Factor a monic polynomial with roots in the closed left half-plane over the reals. The factors have the forms $z+a$ and $z^2+2a z+a^2+b^2$, with $a\ge0$. Their first coefficients below the leading term are nonnegative and their $H$ values are respectively $0$ and $2a(a^2+b^2)\ge0$. Induction gives $H(P/a_d)\ge0$. Multiplying by $a_d^2$ proves the claim. This is the second Hurwitz determinant.

<!-- page -->

## Choosing the graph orders

Write $C(n,s)$ for the binomial coefficient and $b(s)=s(s-1)/2$. Each candidate $s$-set needs $b(s)$ absent edges, so its expected count and the ratio of successive expected counts are

$$
\mu_s=C(n,s)q^{b(s)}, \qquad r_s=\mu_s/\mu_{s-1}=((n-s+1)/s)q^{s-1}.
$$

Fix $\lambda>0$. For each sufficiently large integer $k$, take $n=n_k$ to be the nearest integer to $(\lambda k!)^{1/k}q^{-(k-1)/2}$. Stirling's formula gives $n_k\sim (k/e)q^{-(k-1)/2}$. In particular $k^2/n_k\to0$. Rounding changes the $k$th power by a factor tending to one, and $C(n,k)=(1+o(1))n^k/k!$. Hence $\mu_k\to\lambda$. For every fixed nonnegative integer $j$,

$$
r_k\sim q^{(k-1)/2}/e\to0, \qquad r_{k-j}/r_k\to q^{-j}.
$$

Also $r_{k+1}/r_k\to q$. Thus $\mu_{k+1}\to0$ and $\mu_{k-j}\to\infty$ for $j=1,2,3$. By Markov's inequality, with probability tending to one there is no independent $(k+1)$-set. Every larger independent set would contain one.

**Lemma 2.** On these graph orders, $X_{k-j}/\mu_{k-j}\to1$ in probability for $j=1,2,3$, and $\Pr(X_k>0)\to1-\exp(-\lambda)$. The next two sections prove the lemma.

### Proof of the theorem

Choose $0<\lambda<q^2$. With probability tending to one, the denominator below is positive. Lemma 2 and the exact identity for the ratio of the means give convergence in probability to the indicated limit,

$$
X_{k-1}X_{k-2}/X_{k-3}=(1+o_{\Pr}(1))\mu_k r_{k-2}/r_k \longrightarrow \lambda/q^2 < 1.
$$

Thus the strict inequality $X_{k-1}X_{k-2}<X_{k-3}$ and the equality $X_{k+1}=0$ hold together with probability tending to one. Intersect this event with $X_k\ge1$. The intersection has probability $1-\exp(-\lambda)+o(1)$, even though the events may be dependent. On the intersection the degree is $k$ and

$$
X_{k-1}X_{k-2}<X_k X_{k-3}.
$$

Lemma 1 puts a zero in the open right half-plane. The limsup over all $n$ is at least $1-\exp(-\lambda)$ for every $\lambda<q^2$ with $\lambda>0$. Letting $\lambda$ increase to $q^2$ proves the theorem.

### A consequence for log-concavity

On the same construction, $X_{k-1}^2/X_{k-2}$ converges in probability to $\lambda/q$, since $\mu_{k-1}^2/\mu_{k-2}=\mu_k r_{k-1}/r_k$. Taking $0<\lambda<q$ gives $X_{k-1}^2<X_k X_{k-2}$ with probability at least $1-\exp(-\lambda)-o(1)$. The limsup probability of failure of log-concavity is at least $1-\exp(-q)$, approximately 39.35% at $p=1/2$. This also refutes the log-concavity conclusion of Cameron's Conjecture 6.1.1 [4].

<!-- page -->

## Concentration below the largest size

To prove the first part of Lemma 2, take $s=k-j$ with $j\in\{0,1,2,3\}$. Fix an $s$-set $S$ and choose a second $s$-set $T$ uniformly. The probability of an intersection of size $i$ and its contribution to the normalized second moment are

$$
w_i=C(s,i)C(n-s,s-i)/C(n,s), \qquad h_i=w_i q^{-b(i)}.
$$

An ordered pair with intersection of size $i$ needs precisely $2b(s)-b(i)$ absent edges. Summation gives $\mathbb{E}[X_s^2]/\mu_s^2=\sum_{i=0}^{s}h_i$. In particular $h_s=1/\mu_s$ and $h_0+h_1=w_0+w_1\le1$. We claim that the sum over the intermediate overlaps $2\le i<s$ tends to zero.

**Small overlaps.** Suppose $2\le i\le s/2$. If $T$ has at least $i$ vertices in $S$, it contains some prescribed $i$-subset of $S$. Such a subset is contained in $T$ with probability at most $(s/n)^i$, by sampling without replacement. A union bound gives $w_i\le C(s,i)(s/n)^i\le(s^2/n)^i$. Since $b(i)/i\le s/4$,

$$
h_i \le A^i, \qquad A=(s^2/n)q^{-s/4}.
$$

Here $A=O(k q^{k/4})\to0$ by our choice of $n$. The sum of these terms is bounded by $A^2/(1-A)$ once $A<1$.

**Large overlaps.** Suppose $s/2<i<s$. Put $t=s-i$, so $1\le t<s/2$. Dividing the expression for $h_{s-t}$ by $h_s$ gives

$$
h_{s-t}=\mu_s^{-1} C(s,t)C(n-s,t) q^{t(s-(t+1)/2)}.
$$

Use $C(s,t)C(n-s,t)\le(sn)^t$ and $s-(t+1)/2\ge3s/4-1/2$ to obtain

$$
h_{s-t}\le \mu_s^{-1} B^t, \qquad B=sn q^{3s/4-1/2}.
$$

Now $B=O(k^2 q^{k/4})\to0$. The factor $1/\mu_s$ stays bounded for all four choices of $s$. This range contributes at most $\mu_s^{-1} B/(1-B)$ once $B<1$. The two bounds prove

$$
\sum_{2\le i<s} h_i=o(1), \qquad \operatorname{Var}(X_s)/\mu_s^2\le1/\mu_s+o(1).
$$

For $s=k-1,k-2,k-3$, the means tend to infinity. Chebyshev's inequality gives the three relative concentration statements simultaneously. For $s=k$, the intermediate-overlap estimate still holds and the mean tends to $\lambda$. We use this case next.

<!-- page -->

## Probability that an independent set of size k occurs

List the $m=C(n,k)$ candidate sets in any fixed order. Let $A_i$ be the event that set $i$ is independent, with probability $\pi=q^{b(k)}$. Call two distinct candidates neighbors if they intersect in at least two vertices. Each $A_i$ is independent of the entire collection of its nonneighbors' events because those events depend only on edge indicators outside the internal edges of set $i$.

Write $u_i$ for the probability that none of the first $i$ events occurs, with $u_0=1$. Let $D_i$ be the event that some earlier neighbor occurs and $B_i$ the event that no earlier nonneighbor occurs. The event of no earlier occurrence is $B_i\setminus D_i$. Expanding it and using $\Pr(A_i\cap B_i)=\pi\Pr(B_i)$ gives

$$
u_i-(1-\pi)u_{i-1}=\Pr(A_i\cap B_i\cap D_i)-\pi\Pr(B_i\cap D_i).
$$

The absolute value on the right is at most the sum of $\Pr(A_i\cap A_j)+\pi^2$ over earlier neighbors $j$. For a pair overlapping in $\ell$ vertices, $\Pr(A_i\cap A_j)=\pi^2 q^{-b(\ell)}\ge\pi^2$. Summing the errors and iterating the recursion, whose multiplier $1-\pi$ is at most one, therefore gives

$$
|u_m-(1-\pi)^m|\le\mu_k^2 \sum_{2\le\ell<k} h_\ell=o(1).
$$

Indeed, twice the unordered sum of the joint probabilities is the ordered sum on the right. Since $m\pi=\mu_k\to\lambda$ and $\pi\to0$, we have $(1-\pi)^m\to\exp(-\lambda)$. Taking the complement proves the second part of Lemma 2. One can also obtain it directly from Bollobás's Poisson limit for the number of $k$-cliques in the complement graph [5, Theorem 11.7].

There is also a shorter route to the negative answer. Leaving out the diagonal term $h_k$ in the pair count gives $\mathbb{E}[X_k(X_k-1)]\le\lambda^2+o(1)$. Bonferroni's two-term bound $\Pr(X_k>0)\ge\mu_k-\mathbb{E}[X_k(X_k-1)]/2$ then gives the limsup lower bound $q^2-q^4/2$ by the same coefficient test. At $p=1/2$, this is $7/32$.

### A bound at one fixed graph order

The conclusion also admits an exact finite calculation. Set $p=q=1/2$, $n=260\,000$, $k=30$, and $\delta=1/8$. Let $v_s=\operatorname{Var}(X_s)/\mu_s^2$. If $|X_s-\mu_s|\le\delta\mu_s$ for $s=k-1,k-2,k-3$, then the ratio in the four-coefficient test is at most

$$
R=((1+\delta)^2/(1-\delta)) \mu_{k-1}\mu_{k-2}/\mu_{k-3}<0.74.
$$

From the two-term bound subtract Markov's bound for $X_{k+1}>0$ and three Chebyshev bounds. The resulting lower bound on the probability of a right-half-plane zero is

$$
L=\mu_k-\mathbb{E}[X_k(X_k-1)]/2-\mu_{k+1}-\delta^{-2}\sum_{j=1}^3 v_{k-j}>0.1016.
$$

The companion verifier checks the two strict inequalities using rational arithmetic and the pair-count formula on the preceding page.

<!-- page -->

## An 86-vertex example

Label the vertices $0,\ldots,85$. List the pairs $(u,v)$ with $u<v$ in lexicographic order. Generate 457 bytes of SHAKE256 output from the nine ASCII characters `BC:86:622`. Read each byte from its least significant bit to its most significant bit. Make the pair an edge if its corresponding bit is one. The first 3,655 bits specify the graph.

The entire coefficient sequence, starting at the constant term, is

$$
(1,\ 86,\ 1810,\ 12478,\ 31706,\ 31559,\ 12792,\ 2225,\ 205,\ 21,\ 2).
$$

There are two independent sets of maximum size ten. The second Hurwitz determinant is

$$
21\cdot205-2\cdot2225=-145.
$$

The verifier counts all independent sets by extending them in increasing vertex order. It checks the same polynomial by the separate recurrence $I_G(z)=I_{G-v}(z)+z I_{G-N[v]}(z)$, where $N[v]$ is the closed neighborhood of $v$.

A second exact root certificate locates one zero in $|z-(11/5+9i)|<1/100$. Put $Q(w)=5^{10} I_G((11+45i+w)/5)=\sum_{j=0}^{10}c_jw^j$ and $\rho=1/20$. Exact integer expansion gives

$$
\sum_{j\ne1} (|\Re c_j|+|\Im c_j|)\rho^j < 0.87\,\rho\,\max(|\Re c_1|,|\Im c_1|).
$$

Thus on $|w|=\rho$ the linear term exceeds the sum of the magnitudes of the remaining terms. Rouché's theorem gives exactly one zero in the disk. Its conjugate gives another right-half-plane zero.

Run `python3 verify_random_stability.py` to reproduce both counts, the root disk and the finite probability bound. The same script checks the first and second moments by exhaustive enumeration of all 32,768 labelled six-vertex graphs at $p=1/2$ and all 1,024 labelled five-vertex graphs at $p=1/3$ and $p=2/3$.

### References

[1] J. I. Brown and B. Cameron. On the Stability of Independence Polynomials. Electronic Journal of Combinatorics 25(1), P1.46 (2018), Problem 2. [Journal paper](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v25i1p46/pdf/) · [Cameron's problem list](https://www.smcs.upei.ca/~brcameron/graphpolys.html)

[2] J. I. Brown and R. J. Nowakowski. Average Independence Polynomials. Journal of Combinatorial Theory, Series B 93 (2005), 313–318. [Author's copy](https://jasonibrown.com/pdfs/AvgIndFinal3.pdf)

[3] B. Bollobás and P. Erdős. Cliques in random graphs. Mathematical Proceedings of the Cambridge Philosophical Society 80 (1976), 419–427. [Journal paper](https://doi.org/10.1017/S0305004100053056)

[4] B. Cameron. On the Roots of Independence Polynomials. PhD thesis, Dalhousie University (2019), Conjecture 6.1.1 and Problem 6.2.1. [Thesis](https://dalspaceb.library.dal.ca/server/api/core/bitstreams/c903bf10-142c-4ffe-965a-1927ad52197a/content)

[5] B. Bollobás. Random Graphs, 2nd edition. Cambridge University Press (2001), Chapter 11, Theorem 11.7. [Book chapter](https://doi.org/10.1017/CBO9780511814068.013)
