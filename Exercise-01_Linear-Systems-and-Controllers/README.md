# Exercise 1 — Linear Systems and Linear Controllers

Problem sheet: Exercise 1 of the course (keep a local copy in [`../course-material/`](../course-material/); it is not published here).

## Files

| File | What it contains |
|---|---|
| [`01_Concepts.md`](01_Concepts.md) · [PDF](pdf/01_Concepts.pdf) | Linear algebra toolkit, eigenvalues, discrete-time LTI systems, stability, reachability/controllability, observability, state feedback and tracking — all derived step by step |
| [`02_Solutions.md`](02_Solutions.md) · [PDF](pdf/02_Solutions.pdf) | Eigenvalue practice (8 problems), Task 0 and Task 1, fully worked |
| [`python/task0_reachability_observability.py`](python/task0_reachability_observability.py) | Verifies Task 0: $W_c$, $W_o$, ranks, PBH tests, subspaces, input sequence |
| [`python/task1_tracking_error.py`](python/task1_tracking_error.py) | Verifies Task 1: simulated vs. predicted $x_\infty$ and tracking error |

## Progress

- [x] Eigenvalue practice
- [x] Task 0 — Controllable and observable subspaces
- [x] Task 1 — Tracking error for a state-feedback controller
- [ ] Task 2 — State feedback + reference prefilter
- [ ] Task 3 — State feedback + integrator
- [ ] Numerical experiment — controller comparison

## Key results so far

- **Task 0:** $A$ has eigenvalues $1,1,4$. $\operatorname{rank}W_c = \operatorname{rank}W_o = 2$, so the system is neither reachable/controllable nor observable; $\lambda=1$ is the bad mode in both cases. Reachable subspace: $x_3 = x_1 + x_2$. $[6,2,8]$ is reachable in 2 steps with $u = (1,5)$.
- **Task 1:** $x_\infty - x_{\text{ref}} = (I - A + BK)^{-1}(A - I)x_{\text{ref}}$, which is zero iff $Ax_{\text{ref}} = x_{\text{ref}}$, independent of $K$.
