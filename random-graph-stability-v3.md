# A limit law and a phase diagram for the independence roots of random graphs

Research manuscript · 26 November 2026 · Developed with OpenAI Codex · Unrefereed

Starting revision 423f7949c981ca39a3721e6b3d98acfa94cedbe6. We use the checked coefficient estimates of the [follow-up, Sections 1–6](https://github.com/jeraldhu-yuan/random-graph-independence-stability/blob/423f7949c981ca39a3721e6b3d98acfa94cedbe6/random-graph-stability-follow-up.md), abbreviated **F2**. Let $G\sim G(n,p)$, with fixed edge probability $0<p<1$, and put $q=1-p$. Write $X_s$ for the number of independent sets of size $s$, $I_G(z)=\sum_s X_sz^s$, and $b(j)=j(j-1)/2$. Instability means a zero in the **open** right half-plane. Its probability is $u_n(p)$.

## 1. Main results

**Theorem A (random analytic function and its roots).** Use the deterministic index $k=k(n,q)$ of F2 §2, recalled in §3 below. Put $\mu_s=\binom ns q^{b(s)}$, $r_s=\mu_s/\mu_{s-1}$, $t_n=X_k/\mu_k$ and $B_n(z)=\sum_{j=1}^k r_k^j X_{k-j}z^j/\mu_k$. On every sequence on which $\mu_k\to\lambda\in[0,\infty]$,

$$
(t_n,B_n,t_n+B_n)\ \Longrightarrow\ (T,B_q,T+B_q),\qquad
B_q(z)=\sum_{j\ge1}q^{b(j)}z^j. \tag{1}
$$

The function coordinates have the topology of uniform convergence on compact sets. Here $T=0$ for $\lambda=0$, $T=1$ for $\lambda=\infty$, and $T=Y_\lambda/\lambda$ for $0<\lambda<\infty$, with $Y_\lambda$ Poisson of mean $\lambda$.

Let $Z(f)$ count zeros with multiplicity as a locally finite measure. The actual roots $\zeta$ of $I_G$, mapped to $1/(r_k\zeta)$, converge as a point measure in distribution to

$$
Z(T+B_q)-\mathbf1_{\{T=0\}}\delta_0. \tag{2}
$$

For the degree-$k$ padded polynomial the limit is $Z(T+B_q)$ on the whole plane. On the punctured plane the two limits agree. Counts on bounded Borel sets converge whenever the boundary almost surely contains no limiting roots. For example, a literal assertion about the count of actual roots on $\{0\}$ would be false when $\lambda=0$ if the padded limit were used.

**Theorem B (classification of all fixed edge probabilities).** Set $s=q^4$ and

$$
f_s(x)=\sum_{m\ge0}(-1)^m s^{b(m)}x^m .
$$

Let $Q=s_*^{1/4}$, where $s_*$ is the first positive spectral value of the partial theta function. For $q<Q$ let $x_1<x_2$ be the first two positive roots of $f_s$. At $q=Q$ write $x_1=x_2=x_*$. Define

$$
B=1-f_s(x_1/q^2),\qquad A=1-f_s(x_2/q^2).
\tag{3}
$$

Here $A,B$ are real numbers depending on $q$ and $B_q$ denotes the entire function in (1). There is exactly one $q_c\in(0,Q)$ at which $x_1/q^2=x_2$. At $q_c$ we have $(A,B)=(0,1)$. At $Q$ put $t_*=A=B$ and $\eta=1/t_*$. Write $\bar u_N=N^{-1}\sum_{n\le N}u_n(p)$. Constants in order estimates may depend on the fixed $q$.

| Parameter | Full set of subsequential limits of $u_n(p)$ | Average over orders |
| --- | --- | --- |
| $0<q<q_c$ | $[0,1-\exp(-1/B)]$ | $\bar u_N=\Theta(1/\log N)$ |
| $q=q_c$ | $[0,1-e^{-1}]$ | $\bar u_N=\Theta(\log\log N/\log N)$ |
| $q_c<q<Q$ | $[1-M(A,B),1]$ | $1-\bar u_N=\Theta(1/\log N)$ |
| $q=Q$ | $[1-\eta e^{-\eta},1]$ | $1-\bar u_N=\Theta((\log N)^{-3/2})$ |
| $Q<q<1$ | $\{1\}$ | $\bar u_N\to1$, with $1-\bar u_N=O_q(N^{-1/4})$ |

The middle regime has $0<A<B<1$. Its constant is

$$
M(A,B)=\max_{m\ge1} e^{-m/B}
     \sum_{j=\lceil Am/B\rceil}^{m}\frac{(m/B)^j}{j!}.
\tag{4}
$$

The expression at $m$ in this maximum is at most $(e^{1-1/B}/B)^m$, a geometric bound with base less than one. The formula includes endpoint atoms.

The stable constants of the entire function $t+B_q$, for $t\ge0$, form $[\max(0,A),B]$ for $q\le Q$ and the empty set for $q>Q$. Thus $\inf\{t\ge0:t+B_q\text{ is unstable}\}$ equals $B$ for $q\le q_c$ and zero for $q>q_c$. The onset of instability at $t=0$ precedes the loss of all stable constants.

## 2. Corollaries and numerical landmarks

**Instability at every large order.** The instability probability is bounded away from zero at every sufficiently large order exactly when $p<1-q_c$. It tends to one along all orders exactly when $p<1-Q$. To see strictness in the middle row, the $m=1$ term of (4) is positive, its maximum is attained at a finite $m$, and each included Poisson band omits the positive mass at zero. Thus $0<M(A,B)<1$. At $Q$ one has $0<\eta e^{-\eta}<1$ since $e^\eta>1+\eta$.

**Uniform labelled graphs.** At $q=1/2$ the first row recovers F2 Theorem 1. The substitution $x=q^3y^2$ identifies (3) with the first crossing in F2 §3. In particular the upper limit lies between $0.27149223459$ and $0.27149223460$. For $q>Q$, the last row recovers F2 Theorem 2.

**Four coefficients and log-concavity.** For every $0<q<1$ the limsup of the instability probability is at least $1-e^{-q^2}$. The limsup probability that the independence sequence fails log-concavity is at least $1-e^{-q}$. Choose orders with $\mu_k\to\lambda$ as in §3. The first four coefficients in ascending powers of the reciprocal variable are $t_n,c_1,c_2,c_3$, with $(c_1,c_2,c_3)\to(1,q,q^3)$ in probability. When $\lambda<q^2$, $X_k\ge1$ implies $t_nc_3>c_1c_2$ apart from an event of probability $o(1)$. This excludes closed stability by the [original note, Lemma 1](https://github.com/jeraldhu-yuan/random-graph-independence-stability/blob/423f7949c981ca39a3721e6b3d98acfa94cedbe6/random-graph-stability.md). For $\lambda<q$ the same event gives $t_nc_2>c_1^2$, violating log-concavity. Apply (1) and let $\lambda$ increase to the respective endpoints.

**Corollary C (Makowsky–Rakita Problem 1.1).** For each fixed $0<p<1$, the sequence $(X_s)_{0\le s\le\alpha(G)}$ is unimodal with probability tending to one. The proof follows Theorem B in §9.

**Certified examples.** We obtain $0.726<q_c<0.727$, $0.74572240<Q<0.74572242$, $0.3092493<s_*<0.3092494$ and $0.6816<t_*<0.6820$. Thus $0.273<1-q_c<0.274$ and $0.25427758<1-Q<0.25427760$. The former sufficient-condition endpoint $1-q-q^3=0$ has $0.6823<q<0.6824$ and falls in the first row. At $q=7/10$ the upper crossing is in $[1.1946219,1.1946220]$ and the upper limit is in $[0.5670291,0.5670292]$. At $q=37/50$ the endpoints of the stable interval are in $[0.4164890,0.4164891]$ and $[0.8600525,0.8600526]$. Its graph instability probability has lower limit in $[0.5084363,0.5084364]$. The lower limit at $Q$ lies in $[0.6616,0.6617]$. These bounds use rational interval arithmetic in [the verifier](verify_v3.py).

## 3. Proof of the analytic limit

Put $a_d=(d!)^{1/d}q^{-(d-1)/2}$ and choose $j$ by $a_j\le n<a_{j+1}$. Use $k=j$ if $\mu_j\le n^{1/3}$ and $k=j+1$ otherwise. F2 §2, equations (2)–(4), gives $k\asymp\log n$, $r_{k+i}\asymp k/n$ for fixed $i$, $\mu_k\le n^{1/3}$, $\mu_{k-1}\ge n^{1/3}$, and $\mu_{k+1}=O(kn^{-2/3})$. Define the good event

$$
X_{k+1}=0,\qquad \max_{s<k}|X_s/\mu_s-1|\le k^{-5}. \tag{5}
$$

Its complement has probability $O_q(k^{14}n^{-1/3})$ by Markov and the bound $k^{11}V_{k-1}$ from Chebyshev, where $V_s=\operatorname{Var}(X_s)/\mu_s^2=O(k^3/n+1/\mu_s)$ for $s=k+O(1)$. These bounds are uniform over orders for fixed $q$. The convergence $B_n\to B_q$ holds in probability on each closed disk with all fixed derivatives by F2 §2 and the common coefficient envelope $(1+k^{-5})q^{b(j)}$ on (5). This implies convergence in the compact-open topology by summing the truncated disk distances with weights $2^{-h}$ over integer radii $h$.

For $\lambda=0$, $\Pr(X_k\ne0)\le\mu_k$. For $\lambda=\infty$, $V_k\to0$. At finite positive $\lambda$, apply Bollobás, *Random Graphs*, Theorem 11.7, to cliques in the complement with edge probability $q$. Its hypotheses are a fixed edge probability and a convergent positive clique-count mean. F2 §4 also gives a direct generating-function proof. Replacing the function coordinate by its deterministic limit changes the product distance by $o_{\Pr}(1)$. Addition is continuous. This proves (1).

No limiting function vanishes identically since its linear coefficient is one. For locally uniform convergence to such a function, isolate its finitely many zeros in a bounded continuity domain in small disjoint disks. On each boundary the modulus has a positive minimum, so Rouché gives the multiplicity there. Uniform convergence on the compact remainder excludes extra zeros. Shrinking disks, or testing against compactly supported continuous functions, gives continuity of the zero measure in the vague topology. On (5), $t_n+B_n=(r_kz)^k I_G(1/(r_kz))/\mu_k$ and $X_{k-1}>0$. Exactly one padded zero occurs when $X_k=0$. Its indicator converges jointly to $\mathbf1_{\{T=0\}}$. At finite positive mean isolate the integer zero by a half-integer cutoff, and use the preceding probability bounds in the other cases. This proves (2).

For each fixed $\lambda>0$, orders nearest to $a_d\lambda^{1/d}$ give $k=d$ eventually and $\mu_d=\lambda+O_\lambda(d^2/n)$. Indeed $\mu_d=(n/a_d)^d\prod_{i=0}^{d-1}(1-i/n)$, rounding costs $O(d/n)$ relatively, and the product differs from one by at most $\binom d2/n$. Adjacent means have orders $n/d$ and $d/n$, respectively. The same index rule selects $k=d$ when $\mu_d$ stays in a fixed positive compact interval, or in $[d,d^{3/2}]$ at orders near $a_d$. In these cases $j=d$ or $d-1$ by the ratio of consecutive order scales, and the adjacent means force the stated choice. The mean mesh is multiplicatively $1+O(d/n)$ since $\mu_d(n+1)/\mu_d(n)=(n+1)/(n+1-d)$.


## 4. Controlling the whole right half-plane

**Exterior estimate.** For each fixed $q$ and $K<\infty$ there is an $R$ such that, on (5) for all sufficiently large $n$, none of the polynomials $t+B_n$ has a zero with $\Re z\ge0$ and $|z|\ge R$, simultaneously for $0\le t\le K$. The same assertion for $t+B_q$ is uniform over compact parameter intervals in $0<q<1$.

To prove it, use the **full** mean polynomial $P_n(w)=\sum_{s=0}^n\mu_s w^s$. Its roots are $-d_i$ with $0<d_1<\cdots<d_n$ and $d_i/q<d_{i+1}$. Here is the induction. The recurrence is $P_{n+1}(w)=P_n(w)+wP_n(qw)$ by Pascal's identity. At $x=d_i$ and $x=d_i/q$, $P_{n+1}(-x)$ has sign $(-1)^i$ by the induction hypothesis. Along with the signs at zero and infinity this gives one root in each of $(0,d_1)$, $(d_i/q,d_{i+1})$ and $(d_n/q,\infty)$. The resulting spacing is strict and the initial polynomial is $1+w$.

For $u\ge0$ and $|\phi|\le\pi/2$, $|1+ue^{i\phi}|\ge\max(1,u)$. The logarithmic loss from $1+u$ is therefore at most $\min(u,1/u)$ for $u>0$, by $\log(1+x)\le x$. The sums for $u=|w|/d_i$ on either side of one are each at most $1/(1-q)$ by the spacing. Thus, with $C_q=e^{-2/(1-q)}$,
$$
|P_n(w)|\ge C_q P_n(|w|)\qquad(\Re w\ge0). \tag{6}
$$

Define the Laurent polynomial $D_n(z)=(r_kz)^kP_n(1/(r_kz))/\mu_k$. Its nonnegative-power part is $1+\mathbb EB_n$. Its coefficient at $z^{-h}$ is at most $q^{h(h+1)/2}$, by multiplying $r_{k+i}/r_k\le q^i$ for $1\le i\le h$. For $r=|z|\ge1$, (5) yields
$$
|t+B_n(z)-D_n(z)|\le k^{-5}D_n(r)+K+1/(1-q).
$$
Equation (6) gives $|D_n(z)|\ge C_qD_n(r)$ on the right half-plane, and $D_n(r)\ge r$ since the mean linear coefficient is one. This proves the exterior estimate. At fixed $z\ne0$, $D_n$ converges to the bilateral series $\sum_{j\in\mathbb Z}q^{b(j)}z^j$ by the coefficient bounds on both tails and the exact ratios at fixed indices. Inequality (6) therefore passes to that series and gives the asserted deterministic uniformity, using $\inf C_q>0$ on compact parameter intervals. We shall restrict to bounded $t$ in the exterior estimate. For $t>q^{-2}+1$, finite polynomials are unstable on (5) for all large $n$ by $t>c_1c_2/c_3$.

**First degree correction.** Write $\bar B_n=\mathbb EB_n$. The coefficient of $z^j$, $1\le j\le k$, is exactly
$$
q^{b(j)}\prod_{h=0}^{j-1}(1-h/k)
             \prod_{h=0}^{j-1}(1+h/(n-k+1))^{-1}.
$$
The first product differs from $1-b(j)/k$ by at most $b(j)^2/(2k^2)$ by the two Bonferroni inequalities for a product of numbers in $[0,1]$. The second differs from one by at most $b(j)/(n-k+1)$ by the union bound. Summing against $r^j q^{b(j)}$ gives, on every fixed disk, also for each fixed derivative,
$$
\bar B_n(z)=B_q(z)-\frac{z^2}{2k}B_q''(z)+O_{q,r}(k^{-2}).
\tag{7}
$$
Polynomial factors in $j$ are summable against this envelope and the omitted tail $j>k$ decreases faster than any inverse power of $k$. On (5) replacing $\bar B_n$ by $B_n$ costs $O_{q,r}(k^{-5})$. In particular (7) controls shifts of order $1/k$ at a limiting imaginary-axis root.

## 5. Determining the stable constants of the entire function

Write $C_t=t+B_q$ and $x=q^3 y^2$. Direct separation of powers gives
$$
\Im C_t(iy)=y f_s(x),\qquad
\Re C_t(iy)=t-1+f_s(x/\sqrt s),\qquad
f_s(x)=1-xf_s(sx),\quad s=q^4. \tag{8}
$$
Here $f_s(x)=\theta(s,-x/s)$ for $\theta(s,w)=\sum_{m\ge0}s^{m(m+1)/2}w^m$. All these series have order zero: at radius $r\ge 2$, completing the square in $b(j)\log q+j\log r$ bounds the maximum modulus by $\exp[O_q(1+(\log r)^2)]$, the remaining translated Gaussian sum being uniformly bounded. We use the genus-zero case of Hadamard factorization (Stein–Shakarchi, *Complex Analysis*, Ch. 5, Thm 5.1). Thus the zeros give a convergent genus-zero product, up to a constant and any power at the origin. If $C_t$ is closed stable, wherever $C_t(iy)\ne0$ the derivative of its argument is
$$
\Re(C_t'(iy)/C_t(iy))=\sum_\zeta(-\Re\zeta)/|iy-\zeta|^2>0.
$$
The inequality is strict since there is a zero strictly to the left: otherwise the real entire product, after removing any power of $z$, is even, whereas $C_t$ has consecutive nonzero coefficients. Conjugate factors are paired. The same product argument and Brown–Cameron Thm 2.1 show that the odd part of a closed stable $C_t$ has only real nonpositive zeros in its squared variable: use finite conjugation-invariant subproducts, pass to their odd parts, and use Rouché. For $q>Q$ this is impossible by Kostov Thm 1(3), since the odd part in the squared variable is $\theta(q^4,w/q)$. On a fixed disk around a nonreal zero of this odd part, (5) and Rouché give the same obstruction for the random polynomial, for **all** $t$. This proves the last row of Theorem B, including its bound by (5).

**The higher odd zeros.** The following degree-independent estimates explain why only two crossings matter. Use $0<s\le5/16$. On $|x|=s^{-m-1/2}$ for $m\ge1$, divide $f_s$ by its term of index $m+1$. Writing $x=|x|e^{i\phi}$ gives the sum of $(-1)^\ell s^{\ell^2/2}e^{i\ell\phi}$ for $\ell\ge-(m+1)$. Its real part exceeds $1/20$. Indeed the central five terms are
$$
1-2\sqrt{s}\cos\phi+2s^2(2\cos^2\phi-1).
$$
This quadratic decreases as a function of $\cos\phi$ since $4s^{3/2}<1$. Its minimum is $1-2\sqrt{s}+2s^2$. The remaining absolute sum is at most $2s^{9/2}/(1-s^{7/2})$. The positive margin also follows from the stronger rational bound below. Homotopy to the central monomial on this circle preserves the winding number. There are exactly $m+1$ zeros inside. Thus there are two zeros in the first disk and one simple positive real zero in each subsequent annulus, by conjugation and positivity of $f_s$ on the negative real axis.

For $m\ge2$ that annular zero has the form $x=s^{-m}u$ with $19/20<u<21/20$. To prove this, divide $f_s(s^{-m}u)$ by its term at index $m$. The completed bilateral normalized sum pairs indices $-h,h+1$ and equals
$$
(1-u)\sum_{h\ge0}(-1)^h s^{h(h+1)/2}\sum_{j=-h}^h u^j .
$$
For $19/20\le u\le21/20$, successive inner sums grow by at most $1+u+u^{-1}$ (multiply the previous sum by this number). From $h=1$ the terms decrease. The factor multiplying $1-u$ is at least $1-s(1+u+u^{-1})$. The omitted terms at indices $\ell<-m$ have absolute sum at most $s^6v^3/(1-s^4v)$ where $v=20/19$. At the two proposed endpoints this is less than $|1-u|[1-s(1+v+v^{-1})]$. Both endpoints are in the annulus because $\sqrt s<19/20<21/20<1/\sqrt s$. Hence the sign changes and the localization follows.

At this zero evaluate $f_s(x/\sqrt s)$, now dividing by the term of index $m+1$. Its sign is $(-1)^{m+1}$ because the normalized real value is bounded below by
$$
1-(v+v^{-1})\sqrt s+2s^2-
 \frac{2v^3s^{9/2}}{1-vs^{7/2}}>1/20. \tag{9}
$$
Here the contributions at offsets $\pm2$ total at least $2s^2$ and the tail bound covers all available offsets $|\ell|\ge3$. In (9) the central lower bound is decreasing and the subtracted tail increasing in $s$ on the stated interval. The verifier checks (9) at $s=5/16$ using a rational upper bound for its square root. It checks the endpoint inequality for the omitted indices in the preceding paragraph in the same way. Consequently
$$
|f_s(x/\sqrt s)|\ge(19/20)^{m+1}s^{-(m+1)^2/2}/20>1+2/\sqrt s. \tag{10}
$$
The last inequality holds at $m=2$ after multiplication by $s^{9/2}$ (use $s\le5/16$), and the left side increases with $m$. These bounds are uniform even as the annular index tends to infinity.

**The low zeros and the two transitions.** For $s\le27/100$, $f_s(2)<0$ by the alternating upper bound $-1+4s-8s^3+16s^6$, increasing and negative at $27/100$. The tail decreases from the last retained term since the successive absolute ratio is $2s^j$. At $0$ and at $s^{-3/2}$ the function is positive by its constant and the circle estimate. Thus the first two zeros are distinct and simple. On $(18/25)^4\le s\le31/100$, the rational rectangle certificates of §10 give, for $8/5\le x\le7/2$,
$$
f_s''(x)>1/25,\qquad x f_s''(x)+f_s'(x)>0.
\tag{11}
$$
At the two endpoints $f_s>0$, and $f_s'$ is respectively negative and positive. The minimum in this interval increases strictly with $s$ because $2s\partial_s f_s=x^2f_s''$. Its sign changes once. More explicitly the certificates give a negative value at $s=(0.74572240)^4$ and a positive minimum at $s=(0.74572242)^4$. For the latter, use the quadratic lower tangent at $x=2.32038$ and subtract $|f_s'(x)|^2/(2/25)$ from $f_s(x)$. The two enclosed roots exhaust the first disk whenever the minimum is negative, including its overlap with the range $s\le27/100$ (there $2$ is inside the interval). The first multiple zero is therefore exactly this minimum at zero height, identifying $s_*$ and the stated bracket without assuming a published decimal approximation. The remaining annuli have simple roots by the preceding argument.


For $q\le17/25$ the positive-coefficient ratio bound in F2 §1 applies to all truncations of $C_1$ and $C_0/z$ since $1-q-q^3>0$ there. It gives closed stability in the limit by Rouché, and its imaginary-axis lower bound excludes nonzero boundary roots. At the first positive zero of the imaginary part the derivative of that part is negative, and at the second it is positive. Strictly positive argument derivative at both, applied to $C_1$ and $C_0$, implies $B>1$ and $A<0$. This also implies $x_1/q^2<x_2$: $A<0$ and (8) give $f_s(q^2x_2)<0$, hence $q^2x_2>x_1$. On $17/25\le q\le18/25$ the verifier bounds both $f_s(10/3)$ and $f_s(10q^2/3)$ strictly below zero. These arguments lie in the first disk and hence between its two roots. This gives $x_1/q^2<x_2$ directly.

On $18/25\le q<Q$, implicit differentiation using (11) shows $x_2$ decreases and $x_1/\sqrt s$ increases. Indeed $dx_i/ds=-x_i^2f_s''(x_i)/(2sf_s'(x_i))$ and $x_1 f_s''(x_1)+f_s'(x_1)>0$. The verifier brackets equality between $q=0.726$ and $q=0.727$. Before equality, $x_1/q^2$ and $q^2x_2$ both lie strictly between $x_1,x_2$. Equation (8) gives $B>1$, $A<0$. At equality it gives $B=1,A=0$. After equality $q^2x_2<x_1$, and $x_2<x_1/q^2<x_3$: this last upper bound follows from $x_1\le7/2$, $q\ge18/25$ and $x_3>(19/20)s^{-2}$. Consequently $B<1$ and $A>0$. In the range $0.725\le q\le Q$, both $x_i/q^2$ lie in $[3,13/2]$: at $0.725$ we certify $x_1>1.68$, $x_2<3.3$, then use the monotonicity of the two individual roots and the bounds on $q$. Throughout this rectangle $f'_s>0$ by the verifier. Thus $A<B$ after equality up to coalescence. At $Q$ they agree and $0<t_*<1$. For every $q\le Q$ one also has $0<B<q^{-2}$: $B=x_1 f_s(\sqrt s\,x_1)/\sqrt s$, $1=x_1f_s(sx_1)$, and $f_s$ decreases on $[0,x_1)$ by its product over positive zeros. The same argument includes the double zero.

**From signs to stability.** For $q<Q$ and $\max(0,A)<t<B$, (8)–(10) exclude imaginary-axis zeros. For $t>B$ or $0\le t<A$, the argument derivative at the first or second odd zero respectively is negative, so $C_t$ has a strictly right-half-plane zero. At $q=Q$, $t\ne t_*$ would give argument derivative zero at the double odd zero, also impossible for a closed stable function. To prove sufficiency, start at a small fixed $q$ with $t=1$, where stability was established above. For any destination $q<Q$, follow the continuous midpoint of $(\max(0,A),B)$, connecting to it and to the destination $t$ by segments inside the open interval. For destinations below the starting point use the same path in reverse parameter order. The nonzero boundary values and the uniform exterior estimate of §4 permit the argument principle on one fixed right half-disk for this compact path. Its right-half-plane zero count stays zero. Taking limits supplies all closed endpoints, including $q=Q,t=t_*$ and $t=0$ when allowed. At $t=0$ the origin is simple and is removed when counting actual graph roots. The only other possible boundary roots in the stable cases are the conjugate imaginary pair belonging to each active endpoint, by (8)–(10). They are simple at distinct low odd zeros. At the double odd zero of $q=Q$ they are still simple since the derivative of the **real** part in $x$ is $f_s'(x_*/q^2)/q^2>0$.

## 6. Perturbing boundary roots at finite degree

Work on the event (5). All error statements below hold for every coefficient array on that event. Values $t>q^{-2}+1$ are already excluded. Compact convergence and the exterior estimate give the following test away from a boundary: limiting open right-half-plane roots persist by Rouché, and a closed stable limit without imaginary roots remains strictly stable for all large $n$. If $t\to0$, its simple root near zero is real by conjugation and uniqueness, is zero at $t=0$, and is negative for small $t>0$ since the coefficients are positive. No instability can arise there. The same compactness argument excludes other right-half-plane roots near a boundary, except in small disjoint disks about its active imaginary pair. These conclusions are uniform on compact $t$-sets separated from the active endpoints: otherwise extract a violating convergent subsequence of constants.

At a simple odd crossing, the finite odd-part zero varies by $O_q(k^{-1})$ by (7), uniform derivative convergence and the nonzero odd slope. It gives one and only one finite crossing in its disk. Its real-part threshold is the negative real part of $B_n$ there. The side of stability follows by continuity from the interior, and the opposite side by the negative argument derivative at this finite odd zero. At a nonnegative simple active endpoint $h=A$ or $B$, the new threshold actually expands the stable interval by an amount between $c/k$ and $C/k$, with positive constants for the fixed $q$. For $A=0$ this means the local lower crossing moves to negative $t$; only nonnegative $t$ are tested. Here is the sign and strict size calculation. For the simple zero $\zeta=iy$ of $C_h$, (7) and Taylor expansion give the zero at $t=h$ as
$$
\zeta+\frac{\zeta^2}{2k}\frac{C_h''(\zeta)}{C_h'(\zeta)}+O_q(k^{-2}). \tag{12}
$$
Factoring off $z-\zeta$ in the genus-zero product shows that the real part of $C_h''(\zeta)/(2C_h'(\zeta))$ equals $\sum_{\omega\ne\zeta}(-\Re\omega)/|\zeta-\omega|^2>0$, by the same strictness argument as in §5. Thus the root moves left by a positive multiple of $1/k$. Its real velocity as $t$ varies tends to $\Re(-1/C_h'(\zeta))\ne0$ since the odd slope is nonzero. The implicit function theorem gives a simple crossing displaced by the claimed amount. The local argument applies to real $t$ just below zero as well for the nonzero imaginary pair. At $q_c$, therefore, $t=0$ is stable and the upper finite threshold lies in $[1+c/k,1+C/k]$. For $q<q_c$ the only active threshold on $t\ge0$ is $B+O_q(k^{-1})$ and $t=0$ is stable.

At $q=Q$ there are constants $0<c<C$ such that the finite stable constants **contain** $|t-t_*|\le c/\sqrt k$ and **are contained in** $|t-t_*|\le C/\sqrt k$. To see this, locally solve the real-part equation at $iy$ for $x=q^3 y^2$, denoting its root by $x_e$. Its $x$-derivative is positive and bounded away from zero by §5, so $x_e-x_*$ equals a nonzero constant times $t-t_*$ plus $O((t-t_*)^2+k^{-1})$. Let $O_n(x)=\Im B_n(iy)/y$. At the real-part zero the sign of the real part of the nearby polynomial zero is exactly the sign of $O_n(x_e)$. In detail, if that zero is $w$, let $D$ be the integral of the complex derivative on the segment from $w$ to $iy$, parametrized by unit time. Then $w=iy-iyO_n(x_e)/D$ and $\Im D<0$ since in the limit $-\Im C'_{t_*}(iy)$ is the positive $y$-derivative of the real part. This also treats equality. By (7), uniformly on a neighborhood of $x_*$,
$$
O_n(x)=f_s(x)-[2x^2 f_s''(x)+3xf_s'(x)]/k+O_q(k^{-2}). \tag{13}
$$
Indeed $b(2m+1)=2m(m-1)+3m$. The first term has a quadratic zero with positive second derivative, and the bracket is strictly positive there. Taylor's formula gives the two bands. Outside a fixed neighborhood use compactness as above. This explains why the stable constants at a double odd zero have positive finite-degree width, although the limiting set is a singleton.


## 7. Extremal probabilities and all intermediate limits

Write $\mu=\mu_k$. In addition to (1) we use the quantitative part of Bollobás, *Random Graphs*, Thm 11.9: at fixed edge probability, for clique-count means bounded away from zero and $o(n^2/(\log n)^4)$ the total variation distance from Poisson of the same mean is $O(\mu(\log n)^4/n^2+(\log n)^3/n)$. We use it along subsequences where the mean tends to a positive number or to infinity, as in the theorem. All uses for $c\le\mu_k\le n^{1/3}$ are uniform for fixed $c>0$: a sequence contradicting the bound would have a subsequence with a limit in $[c,\infty]$ contradicting that estimate. A Poisson variable $Y_\mu$, centered and scaled by $\sqrt\mu$, tends to the standard normal as $\mu\to\infty$, by the characteristic function $\exp[\mu(e^{iu/\sqrt\mu}-1-iu/\sqrt\mu)]$.

**The first two rows.** If $\mu\to0$, stability at $t=0$ and Markov give limit zero. If $\mu\to\infty$ and $B>1$, the probability again tends to zero by Theorem A and §6. At $q=q_c$ and $\mu\to\infty$, it is at most $1/2+o(1)$ by the threshold above one, the Poisson approximation and the normal limit. For $\mu\to\lambda\in(0,\infty)$ the upper bound is $\Pr(Y_\lambda\ge\lambda B)$, while the strict inequality event gives a lower bound, by §6 and convergence of $t_n$. If $\lambda B$ is nonintegral the bounds coincide.

For $a=1/B\le1$ and $m=\lceil\lambda B\rceil$, monotonicity in the Poisson mean bounds the upper probability by $G_m(a)=\Pr(Y_{ma}\ge m)$. In fact $G_m(a)\le G_1(a)$. Differentiating the Poisson series gives $G'_m(a)=m^m a^{m-1}e^{-ma}/(m-1)!$. The ratio $G'_{m+1}/G'_m$ is a constant times $ae^{-a}$, increasing on $[0,1]$, so the ratio of their integrals from zero is increasing there too (it is a weighted average of this increasing ratio). At $a=1$ compare
$$
G_{m+1}(1)-G_m(1)=\int_m^{m+1}e^{-x}x^m/m!\,dx-e^{-m}m^m/m!\le0
$$
since the integrand is maximized at $m$. This proves the bound for all $a\le1$. In particular the upper limit at $q_c$ also dominates $1/2$, since $e>2$ by its power series. Every value in $(0,1-e^{-1/B})$ is attained by taking $0<\lambda<1/B$ in §3, when the limit is $1-e^{-\lambda}$. The closedness of a bounded sequence's set of subsequential limits gives both endpoints.

**The middle and coalescence rows.** Here both $t=0$ and $t=1$ are unstable. At $\mu\to0,\infty$ the instability probability tends to one. For finite positive $\lambda$ the upper bound on *stability* is $\Pr(A\lambda\le Y_\lambda\le B\lambda)$, by closed-set convergence and §6. This applies also when $A=B$. For $\lambda$ with $\lfloor B\lambda\rfloor=m\ge1$, moving the mean down to $m/B$ can only increase this bound: the upper index stays $m$, the lower index cannot increase, and every included mass $e^{-\lambda}\lambda^j/j!$ decreases with $\lambda$ because $j\le m<\lambda$. If $m=0$ the probability is zero since $A>0$. Markov applied to $B^{Y_{m/B}}$ gives the tail majorant in (4); its base is less than one since $\log x<x-1$ for $x=1/B>1$. Hence the maximum exists. At $Q$ the mass at $m$ is $e^{-m\eta}(m\eta)^m/m!$, decreasing with $m$: the successive ratio is $e^{-\eta}\eta(1+1/m)^m<e^{1-\eta}\eta<1$. This gives $\eta e^{-\eta}$.

The upper bounds on stability are attained even for coincident endpoint atoms. Take the tuned orders of §3 with $\lambda=m/B$. The deterministic error in $\mu$ is $O_q(k^2/n)$. For each integer count in the closed band, the scaled constant differs from its target by $O_q(k^2/n)$, so all such counts are stable on (5) for large $n$, by the inward correction (12) or the inner band from (13). Counts outside the band have a strictly unstable limit for each fixed count. Tightness of the counts gives the claimed attainment.

**No gaps between the extrema.** This step is needed because a Poisson law can have an atom at a crossing. We show $|u_{n+1}-u_n|\to0$ for $q_c<q\le Q$. Couple the two graphs by deleting the last vertex. First take $\mu\to\lambda\in(0,\infty)$; the same $k$ applies at both orders for large $n$. Put $\sigma=\mu_{k-1}^{-1/2}\asymp_q\sqrt{k/n}$ along this subsequence (constants can depend on $\lambda$). On compact disks, with each fixed derivative,
$$
B_n(z)=\bar B_n(z)+(c_1-1)z+o_{\Pr}(\sigma),\qquad
(c_1-1)/\sigma\ \Longrightarrow\ N(0,1). \tag{14}
$$
Indeed $c_1=X_{k-1}/\mu_{k-1}$, and $\mu_{k-1}\asymp n/k$ satisfies Thm 11.9. For the other coefficients refine the overlap estimate of F2 §2 at $s'=k-2$. Overlap two costs $O_q(k^4/n^2)$ by the union bound for two prescribed vertices. For $3\le i\le s'/2$ its geometric majorant sums to $O_q(k^{9/2}/n^{3/2})$ (F2's $A^3/(1-A)$). The upper-overlap sum is $O_q(1/\mu_{k-2})=O_q(k^2/n^2)$ since the adjacent ratios are of order $k/n$. Therefore $kV_{k-2}=o(\sigma^2)$. The monotonicity of $V_s$ and Chebyshev give simultaneous relative errors $o_{\Pr}(\sigma)$ through $k-2$. The common Gaussian coefficient envelope proves (14).

For bounded $X_k=m$, only a limiting boundary value $m/\lambda$ requires a comparison. In a neighborhood of its active root solve $m/\mu+\Re\bar B_n(iy)=0$ for the **deterministic** abscissa $\bar x=q^3y^2$. Its derivative in $x$ is bounded below by a positive constant, by the rectangle $3\le x_i/q^2\le13/2$ in §5. For the actual even part the root differs by $o_{\Pr}(\sigma)$, since the random linear term in (14) is purely imaginary. Let $\bar O_n$ denote the mean odd part divided by $y$ on this axis. The sign deciding instability at that root, by the segment argument in §6, is the sign of
$$
\bar O_n(\bar x)+(c_1-1)+o_{\Pr}(\sigma). \tag{15}
$$
This covers a double odd zero as well. There is only one conjugate pair to test at a given endpoint in these regimes. The estimates can be used simultaneously for finitely many specified $m$, without conditioning on the count.

Under the deletion coupling, $X_k$ changes with probability $O(k/n)$ since its expected increment is $O(k/n)$. The expected increment of $X_{k-1}$ is $O(1)$ by the exact adjacent-order mean ratio in §3. Thus $c_1$ changes by $O_{\Pr}(k/n)=o_{\Pr}(\sigma)$ including the change in its normalizer. The deterministic term in (15) changes by $O_q(k/n)$: at fixed $k$ the only $n$-dependent product in §4 differs from one by at most $b(j)/(n-k+1)$ at either order, so the mean analytic functions and derivatives differ by $O(1/n)$ on this disk. The change of $m/\mu$ is $O(k/n)$, and the even derivative has a fixed positive lower bound. By (14) the probability that $c_1-1$ lies within $o_{\Pr}(\sigma)$ of **any deterministic moving center** tends to zero. To see this, bound the random width by $\varepsilon\sigma$ with probability tending to one, take convergent subsequences of the scaled centers in the extended reals, use continuity and tightness of the normal law, and let $\varepsilon\downarrow0$. No independence of $c_1$ and $X_k$ is assumed. Thus (15) changes sign with probability $o(1)$. Truncating $m$ by tightness proves the adjacent-order claim at finite mean. If both means escape finite positive limits, both orders have instability probability tending to one. A finite positive subsequential mean at either order forces the same index and mean limit at its neighbor by §3. This exhausts subsequences. Finally, a bounded real sequence whose consecutive differences tend to zero takes subsequences tending to each value strictly between its lower and upper limits, by first passage between successive excursions. This proves the interval assertions including the singleton-band case.


## 8. Averaging over graph orders

Here are the counting details, including the possibly incomplete last block. If $k=d$ then $a_{d-1}\le n<a_{d+1}$, so $n\asymp_q a_d$. By the exact mesh in §3, an interval of length $w$ in $\log\mu_d$ contributes at most $C_q(1+a_dw/d)$ orders with this index. Conversely, for log-means between fixed constants or of magnitude $O(\log d)$, the interval corresponds to orders $n/a_d\to1$. In the windows used below $k=d$ by §3; if $w\gg d/a_d$ their counts are $\asymp_q a_dw/d$. To verify the inverse statement, use monotonicity of the binomial mean, the product estimate in §3 at $n=a_d\exp(O(\log d)/d)$ and the upper and lower mesh bounds there. Errors from (5) contribute $O_q(N^{3/4})$ in total since powers of $\log n$ are absorbed into $n^{1/12}$. The approximation of $X_k$ by Poisson when $\mu_k\ge c>0$ costs at most the same total by Thm 11.9 and $\mu_k\le n^{1/3}$.

For $q<q_c$, instability requires at least one set when $\mu\le1$, at a cost at most $\mu$ plus the exceptional probability. For $\mu\ge1$ the threshold is at least $1+\epsilon$ for some fixed $\epsilon>0$, so the cost is $O_q(1/\mu)$ by the variance of a Poisson variable. Sum over dyadic intervals of $\mu$ on both sides of one. Each costs at most $C_q(1+a_d/d)$ times a summable geometric weight, giving $O_q(a_d/d)$. For a lower bound choose any positive closed interval of means of positive length inside $(0,1/B)$. There the instability probability is uniformly bounded below, by the finite-mean limit in §7 and compactness.

At $q=q_c$ the same small-mean bound holds. On $1\le\mu\le d^2$ use probability at most one, costing $O_q(1+a_d\log d/d)$. Beyond $d^2$, Chebyshev at the threshold $1+c/d$ costs $O_q(d^2/\mu)$, summable dyadically starting there. Conversely on $d\le\mu\le d^{3/2}$ the upper threshold $1+C/d$ is less than $1+o(1/\sqrt\mu)$ uniformly. The normal limit from §7 gives a uniform positive lower bound on instability by subsequence compactness. The logarithmic width of these means is $(\log d)/2$, proving the extra factor.

For $q_c<q<Q$ estimate stability instead. On (5) it requires $t_n$ to lie in a fixed interval $[\alpha,\beta]\subset(0,1)$ containing $[A,B]$ in its interior, so $\mu\ge1/\beta$. For a Poisson variable the bound at most $\beta\mu$ is $\exp[-\mu(1-\beta+\beta\log\beta)]$ by Markov applied to $\beta^{Y_\mu}$. The exponent coefficient is positive by differentiation or strict convexity. Summing over dyadic intervals gives $O_q(a_d/d)$. A fixed positive-length compact interval inside $(1/B,1/A)$ gives the lower bound by $X_k=1$.

At $Q$ let $\epsilon=C/\sqrt d$. For large $d$, §6 restricts each stable integer count $m$ to $m\ge1$ and
$$
m/(t_*+\epsilon)\le\mu\le m/(t_*-\epsilon).
$$
Each such interval has logarithmic width $O_q(\epsilon)$ by the mean value theorem. The Poisson probability of this count is at most $e^{-c_qm}$ by the same lower-tail bound with any fixed $\beta\in(t_*,1)$ and large $d$. Sum the mesh estimate weighted by $e^{-c_qm}$ to obtain $O_q(1+a_d/d^{3/2})$. The errors of (5) and Poisson approximation are charged once per order, before summing the exact Poisson masses. Conversely use $m=1$ and the narrower band $|1/\mu-t_*|\le c/\sqrt d$ in §6. Its logarithmic width is $\asymp_q d^{-1/2}$ and the Poisson probability of one stays bounded below.

For $a_J\le N<a_{J+1}$ at most indices through $J+1$ contribute. Because $a_{d+1}/a_d\to q^{-1/2}>1$, the upper block estimates sum to their last order of magnitude: split at $J/2$, bound the early part geometrically and sum a geometric progression times comparable powers of $J$ in the later part. A lower window with $d=J-1$ lies entirely below $N$, has $a_d\asymp_q N$ and $d\asymp_q\log N$. Divide by $N$ to obtain the first four rates. For $q>Q$ summing the exceptional probability from §5 gives the last bound directly. This completes Theorem B.

## 9. Unimodality

On (5), for $1\le s\le k-2$ the neighboring coefficient ratio $X_s^2/(X_{s-1}X_{s+1})$ is at least $q^{-1}[(1-k^{-5})/(1+k^{-5})]^2>1$ eventually, since $r_s/r_{s+1}\ge q^{-1}$. Thus the positive prefix through $k-1$ is log-concave, hence unimodal by monotonicity of the consecutive ratios. Moreover $r_{k-1}\asymp k/n\to0$, so on the same event $X_{k-1}<X_{k-2}$ for large $n$. Markov gives $\Pr(X_k>\mu_{k-1}/2)\le2r_k=O_q(k/n)$, whereas $X_{k-1}>\mu_{k-1}/2$ on (5). Appending the last term (or discarding a terminal zero) preserves unimodality. This proves Corollary C, while the comparison of the last *ratios* can still fail as in §2.

## 10. Certificates, feasible experiments and prior work

The standard-library verifier uses outward integer intervals on a grid of spacing $10^{-34}$ for the compact rectangles and exact rational alternating remainders for point bisections. For $s\le5/16$, $0\le x\le7$ and derivative orders through two, the first omitted term after degree 15 is bounded by $16^2 7^{16}(5/16)^{120}$. The next ratios are at most $7(5/16)^{16}(17/15)<1/2$, so twice the first bound covers the whole tail. The signs in (11), in the $3\le x\le13/2$ derivative rectangle, and at the transition brackets are strict interval assertions on closed boxes whose union covers the specified sets. In (9) the term $s^{7/2}$ controls each ratio after offset three. For (4) the code encloses all included Poisson masses by directed bounds at both ends of the mean interval, allows outward rounding of both integer cutoffs, and bounds the omitted maximum geometrically. The checked F2 verifier supplies the rational uniform-graph constant. Run **python3 -B verify_v3.py** from the repository directory.

The verifier also uses a second upper bound on omitted indices in (4). For $a\ge1$ let $H_m(a)=\Pr(Y_{ma}\le m)$. Then $H_{m+1}(a)\le H_m(a)$. Indeed, $H_m(a)$ is the integral from $a$ to infinity of $h_m(x)=m^{m+1}x^me^{-mx}/m!$. The density ratio $h_{m+1}/h_m$ is a constant times $xe^{-x}$, decreasing on $[1,\infty)$, so the ratio of the tail integrals decreases there by weighted averaging. At $a=1$ their difference is
$$
H_{m+1}(1)-H_m(1)=e^{-m}\frac{m^{m+1}}{(m+1)!}
 -\int_m^{m+1}\frac{e^{-x}x^{m+1}}{(m+1)!}\,dx\le0
$$
because the integrand increases on that interval. If $B\le b_u<1$, every term of (4) with $m\ge K$ is at most $H_K(1/b_u)$. The verifier takes the smaller of this and the geometric tail bound.

**Non-rigorous Monte Carlo.** At each pair of parameters we generated 128 labelled graphs by SHAKE256, counted every independent set and computed integer Hurwitz determinants. At 36 vertices all 128 graphs were strictly stable at each of $p=17/64$ and $p=19/64$. At 72 vertices, 2 of 128 were unstable at $p=17/64$ and 6 of 128 at $p=19/64$. All remaining graphs in those samples were strictly stable. The 72-vertex samples use seeds 1 through 128, excluding a timing pilot at seed 0. The byte thresholds give the specified edge probabilities for uniform input bytes. The labels are exact for the sampled graphs. The samples provide no estimate of an asymptotic extremum.

**Prior work.** Bollobás and Erdős studied clique sizes and overlaps in the independent-edge model. Bollobás's *Random Graphs*, second edition, Ch. 11, Thms 11.7 and 11.9, pp. 291–292, give respectively the finite-mean law and the growing-mean total variation bound used here. The random coefficient at the top of (1) is exactly that classical clique count in the complement. [Bollobás–Erdős](https://doi.org/10.1017/S0305004100053056) · [Bollobás, Chapter 11](https://doi.org/10.1017/CBO9780511814068.013)

Kostov's *A property of a partial theta function*, Thm 1(3), puts a nonreal conjugate pair in the odd part for every $q^4>s_*$, including later spectral values. Brown and Cameron's *On the Stability of Independence Polynomials*, Thm 2.1, gives the weak Hermite–Biehler equivalence for standard real polynomials, with nonpositive squared-variable zeros and the odd part interlacing or alternating with the even part. Their Problem 2 asks about almost-sure stability. The two cited theorems, combined with coefficient concentration, yield the last phase. The first two odd zeros and the motion of the even part in (8) are needed for the other phases. [Kostov, Thm 1](https://arxiv.org/pdf/1504.01524) · [Brown–Cameron, Thm 2.1 and Problem 2](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v25i1p46/pdf/)

Brown and Nowakowski's *Average Independence Polynomials*, Prop. 1 and Thm 2, concerns nonreal roots of individual uniform graphs and real zeros of their average. Brown, Dilcher and Manna studied real zeros of the expected polynomial at other fixed edge probabilities. Section 4 proves the separation estimate it needs directly. Makowsky and Rakita's Problem 1.1 asks about unimodality for uniform graphs, in the sense of their Definition 1.5. Section 9 covers every fixed edge probability. [Brown–Nowakowski](https://jasonibrown.com/pdfs/AvgIndFinal3.pdf) · [Brown–Dilcher–Manna](https://doi.org/10.1002/jgt.21678) · [Makowsky–Rakita](https://arxiv.org/pdf/2102.00268)

Among these sources, the scalar clique law alone gives no random analytic function or global zero classification. Theorem A is a joint-limit consequence of the already checked F2 §§2,4 and local zero counting, with the padding correction (2). F2 already proves two special regimes of Theorem B. We found no prior result giving the full table, its two different critical averaging rates, or the intermediate stable band. This is a scoped literature assessment. The analytic factorization used in §5 is the classical finite-order product theorem. [Stein–Shakarchi, *Complex Analysis*, Ch. 5, Thm 5.1](https://cs.mcgill.ca/~akroit/math/analysis/Stein%20and%20Shakarchi%20Complex%20Analysis.pdf)
