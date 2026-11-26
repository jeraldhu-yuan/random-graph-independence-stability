# Instability of independence polynomials of random graphs

Brown and Cameron asked whether the independence polynomial of a random graph has all its zeros in the closed left half-plane with probability tending to one. This unrefereed draft, developed with OpenAI Codex, proposes a negative answer. [Original paper, Problem 2](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v25i1p46/pdf/)

For a uniformly random labelled graph, the manuscript gives a lower bound of $1-e^{-1/4}$ on the limsup probability of a zero with strictly positive real part. Consequently, more than 22% of graphs are unstable at arbitrarily large graph orders. The manuscript also gives a bound for every fixed edge probability between zero and one.

[Read the proof](random-graph-stability.md) · [Download the PDF](random-graph-stability.pdf)

## Reproduce the finite checks

Run the verifier from the repository directory with Python 3.8 or later. It uses the standard library.

```sh
python3 verify_random_stability.py
```

The script counts the independent sets of the specified 86-vertex graph by two methods, checks the Hurwitz determinant and a disk containing a zero, and verifies a finite probability bound at 260,000 vertices. The manuscript gives the proof of the asymptotic claim. [Verifier source](verify_random_stability.py)
