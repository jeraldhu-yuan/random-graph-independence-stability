# Instability of independence polynomials of random graphs

Brown and Cameron asked whether the independence polynomial of a random graph has all its roots in the closed left half-plane with probability tending to one. These unrefereed drafts were developed with OpenAI Codex. [Brown and Cameron, Problem 2](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v25i1p46/pdf/)

## How often is a random graph unstable?

Let $u_n$ be the probability that the independence polynomial of a uniformly random labelled graph on $n$ vertices has a root with strictly positive real part. The follow-up derives $\liminf u_n=0$ and $\limsup u_n=0.27149223459\ldots$. Every value between the two endpoints is a subsequential limit. The mean of $u_n$ over orders $1$ through $N$ is of order $1/\log N$.

At smaller fixed edge probabilities, including $p=0.1$, the follow-up proves that the probability of instability tends to one. Theorem 2 states the range of $p$.

[Read the follow-up](random-graph-stability-follow-up.md) · [Download the PDF](random-graph-stability-follow-up.pdf) · [Independent referee prompt](follow-up-review-prompt.txt)

Run the exact checks with Python 3 and its standard library. The second command also reproduces the sample of 8,192 graphs at 86 vertices.

```sh
python3 -B verify_stability_followup.py
python3 -B verify_stability_followup.py --sample-86 8192
```

## Earlier four-coefficient bound

For a uniformly random labelled graph, the manuscript gives a lower bound of $1-e^{-1/4}$ on the limsup probability of a zero with strictly positive real part. Consequently, more than 22% of graphs are unstable at arbitrarily large graph orders. The manuscript also gives a bound for every fixed edge probability between zero and one.

[Read the proof](random-graph-stability.md) · [Download the PDF](random-graph-stability.pdf)

### Reproduce the earlier finite checks

Run the verifier from the repository directory with Python 3.8 or later. It uses the standard library.

```sh
python3 verify_random_stability.py
```

The script counts the independent sets of the specified 86-vertex graph by two methods, checks the Hurwitz determinant and a disk containing a zero, and verifies a finite probability bound at 260,000 vertices. The manuscript gives the proof of the asymptotic claim. [Verifier source](verify_random_stability.py)
