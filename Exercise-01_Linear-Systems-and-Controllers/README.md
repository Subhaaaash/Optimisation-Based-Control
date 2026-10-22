# Exercise 1 — Linear Systems and Linear Controllers

Problem sheet: Exercise 1 of the course (keep a local copy in [`../course-material/`](../course-material/); it is not published here).

## Files

| File | What it contains |
|---|---|
| [`01_Concepts.md`](01_Concepts.md) · [PDF](pdf/01_Concepts.pdf) | Linear algebra toolkit, eigenvalues, discrete-time LTI systems, stability, reachability/controllability, observability, state feedback and tracking — all derived step by step |
| [`02_Solutions.md`](02_Solutions.md) · [PDF](pdf/02_Solutions.pdf) | Eigenvalue practice (8 problems), Tasks 0–3 and the numerical experiment, fully worked (checked against the official solution) |
| [`03_Lab.md`](03_Lab.md) · [PDF](pdf/03_Lab.pdf) | Hands-on lab: every concept shown on a simulated planar quadcopter and robot arm, with plots and animations |
| [`python/task0_reachability_observability.py`](python/task0_reachability_observability.py) | Verifies Task 0: $W_c$, $W_o$, ranks, PBH tests, subspaces, input sequence |
| [`python/task1_tracking_error.py`](python/task1_tracking_error.py) | Verifies Task 1: simulated vs. predicted $x_\infty$ and tracking error |
| [`python/task2_prefilter.py`](python/task2_prefilter.py) | Verifies Task 2: prefilter formulas (1-D and matrix), model sensitivity |
| [`python/task3_integrator.py`](python/task3_integrator.py) | Verifies Task 3: augmented model, zero error, the Lyapunov-stable cases of 3.3* |
| [`python/numerical_experiment.py`](python/numerical_experiment.py) | Numerical experiment Q1–Q4: analysis, designs of $K, F, K_p, K_I$, step responses, robustness, windup → [`python/figures/`](python/figures/) |

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
| [`exp7_prefilter_vs_integrator.py`](lab/exp7_prefilter_vs_integrator.py) | Tasks 2 and 3: prefilter vs integrator, unknown payload | [png](lab/figures/exp7_prefilter_vs_integrator.png) · [gif](lab/figures/exp7_payload.gif) |

## Progress

- [x] Eigenvalue practice
- [x] Task 0 — Controllable and observable subspaces
- [x] Task 1 — Tracking error for a state-feedback controller
- [x] Lab — drone and robot-arm experiments for Tasks 0–3
- [x] Task 2 — State feedback + reference prefilter
- [x] Task 3 — State feedback + integrator (incl. 3.3*)
- [x] Numerical experiment — controller comparison

## Key results

- **Task 0:** $A$ has eigenvalues $1,1,4$. $\operatorname{rank}W_c = \operatorname{rank}W_o = 2$, so the system is neither reachable/controllable nor observable; $\lambda=1$ is the bad mode in both cases. Reachable subspace: $x_3 = x_1 + x_2$. $[6,2,8]$ is reachable in 2 steps with $u = (1,5)$.
- **Task 1:** $x_\infty - x_{\text{ref}} = (I - A + BK)^{-1}(A - I)x_{\text{ref}}$, which is zero iff $Ax_{\text{ref}} = x_{\text{ref}}$, independent of $K$.
- **Task 2:** zero error iff $(A - BK + BKF)x_{\text{ref}} = x_{\text{ref}}$; 1-D: $f = (1 - a + bk)/(bk)$; matrix: $F = (BK)^{-1}(I - A) + I$ (needs $\operatorname{rank}B = n$). Works only with an exact model.
- **Task 3:** augmented matrix $\tilde A = \begin{bmatrix} A - BK_p & -BK_I \\ T_sI & I\end{bmatrix}$; if it is asymptotically stable, the integrator row forces $x_\infty = x_{\text{ref}}$, robust to model errors. With $m < n$, $\lambda = 1$ is always an eigenvalue (windup).
- **Numerical experiment:** a), b): Tasks 2 and 3 reach zero error, Task 1 does not. c) DC motor: $[1, 0]$ is not holdable by any constant input; $x_1$ can be made exact, $x_2$ stays $9\cdot10^{-5}$ off.
