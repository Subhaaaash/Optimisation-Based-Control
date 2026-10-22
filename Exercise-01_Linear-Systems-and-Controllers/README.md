# Exercise 1 — Linear Systems and Linear Controllers

Problem sheet: Exercise 1 of the course (keep a local copy in [`../course-material/`](../course-material/); it is not published here).

## Files

| File | What it contains |
|---|---|
| [`01_Concepts.md`](01_Concepts.md) · [PDF](pdf/01_Concepts.pdf) | Linear algebra toolkit, eigenvalues, discrete-time LTI systems, stability, reachability/controllability, observability, state feedback and tracking — all derived step by step |
| [`02_Solutions.md`](02_Solutions.md) · [PDF](pdf/02_Solutions.pdf) | Eigenvalue practice (8 problems), Task 0 and Task 1, fully worked |
| [`03_Lab.md`](03_Lab.md) · [PDF](pdf/03_Lab.pdf) | Hands-on lab: every concept shown on a simulated planar quadcopter and robot arm, with plots and animations |
| [`python/task0_reachability_observability.py`](python/task0_reachability_observability.py) | Verifies Task 0: $W_c$, $W_o$, ranks, PBH tests, subspaces, input sequence |
| [`python/task1_tracking_error.py`](python/task1_tracking_error.py) | Verifies Task 1: simulated vs. predicted $x_\infty$ and tracking error |

## Lab experiments ([`lab/`](lab/))

Run all with `python lab/run_all.py`, or one at a time (add `--show` for interactive windows).

| Script | Concept | Main figure |
|---|---|---|
| [`exp1_eigenvalues_and_modes.py`](lab/exp1_eigenvalues_and_modes.py) | eigenvalues, modes, $\lambda^k$ | [png](lab/figures/exp1_eigenvalue_modes.png) |
| [`exp2_drone_open_loop.py`](lab/exp2_drone_open_loop.py) | stability, repeated eigenvalues | [png](lab/figures/exp2_drone_open_loop.png) · [gif](lab/figures/exp2_drone_drift.gif) |
| [`exp3_controllability.py`](lab/exp3_controllability.py) | $W_c$, reachable subspace, PBH, minimum steps | [png](lab/figures/exp3_controllability.png) · [gif](lab/figures/exp3_steering.gif) |
| [`exp4_observability.py`](lab/exp4_observability.py) | $W_o$, unobservable subspace, state reconstruction | [png](lab/figures/exp4_observability.png) · [gif](lab/figures/exp4_same_sensor_different_reality.gif) |
| [`exp5_state_feedback.py`](lab/exp5_state_feedback.py) | $u=-Kx$, pole placement, closed-loop stability | [png](lab/figures/exp5_state_feedback.png) · [gif](lab/figures/exp5_feedback_recovery.gif) |
| [`exp6_tracking_error.py`](lab/exp6_tracking_error.py) | Task 1 tracking error | [png](lab/figures/exp6_tracking_error.png) · [gif](lab/figures/exp6_arm_sag.gif) |

## Progress

- [x] Eigenvalue practice
- [x] Task 0 — Controllable and observable subspaces
- [x] Task 1 — Tracking error for a state-feedback controller
- [x] Lab — drone and robot-arm experiments for Tasks 0–1
- [ ] Task 2 — State feedback + reference prefilter
- [ ] Task 3 — State feedback + integrator
- [ ] Numerical experiment — controller comparison

## Key results so far

- **Task 0:** $A$ has eigenvalues $1,1,4$. $\operatorname{rank}W_c = \operatorname{rank}W_o = 2$, so the system is neither reachable/controllable nor observable; $\lambda=1$ is the bad mode in both cases. Reachable subspace: $x_3 = x_1 + x_2$. $[6,2,8]$ is reachable in 2 steps with $u = (1,5)$.
- **Task 1:** $x_\infty - x_{\text{ref}} = (I - A + BK)^{-1}(A - I)x_{\text{ref}}$, which is zero iff $Ax_{\text{ref}} = x_{\text{ref}}$, independent of $K$.
