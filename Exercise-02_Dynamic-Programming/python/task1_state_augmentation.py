"""
Optimization-Based Control Methods - Exercise 2, Task 1
(c)(iv) Slew-rate penalty  g = u_k^2 + rho (u_k - u_{k-1})^2  handled by state augmentation,
(c)(i)/(ii) hard barrier vs soft dead-zone penalty under disturbances,
(b) causality: which cost terms a disturbance can influence.

Run:  python task1_state_augmentation.py
"""
import itertools

import numpy as np

# ---------------------------------------------------------------------------
# (c)(iv) augmented-state DP == brute force
# ---------------------------------------------------------------------------
print("=== (c)(iv) slew-rate penalty via augmented state x~ = (x, u_prev) ===")
a, b, rho, N = 1.1, 1.0, 4.0, 4
X = np.arange(-6, 7)            # integer state grid
U = np.arange(-3, 4)            # integer inputs
x0, u_prev0 = 5, 0              # u_{-1} = 0 (actuator at rest before start)


def f(x, u):
    return int(np.clip(round(a * x + b * u), X.min(), X.max()))   # deterministic, stays on grid


def stage(x, u, up):
    return x**2 + u**2 + rho * (u - up) ** 2


# DP over the augmented state (x, u_prev): this is an ordinary Markov DP problem
J = {(x, up): x**2 for x in X for up in U}                 # terminal cost g_N(x_N) = x_N^2
for k in range(N - 1, -1, -1):
    J = {(x, up): min(stage(x, u, up) + J[(f(x, u), u)] for u in U) for x in X for up in U}
J_dp = J[(x0, u_prev0)]


# brute force over all input sequences
def total(useq):
    x, up, c = x0, u_prev0, 0.0
    for u in useq:
        c += stage(x, u, up)
        x, up = f(x, u), u
    return c + x**2


J_bf = min(total(s) for s in itertools.product(U, repeat=N))
print(f"augmented DP J0 = {J_dp},  brute force = {J_bf}  ->  equal: {J_dp == J_bf}")

# what goes wrong if u_{k-1} is ignored: a DP on x alone can only use rho*(u - ?)^2 with a guess
print("A DP on x alone cannot even evaluate the stage cost: u_{k-1} is not part of its state.\n")

# ---------------------------------------------------------------------------
# (c)(i)/(ii) barrier vs soft penalty when a disturbance pushes x past xbar
# ---------------------------------------------------------------------------
print("=== (c)(i)/(ii) hard barrier vs soft dead-zone penalty ===")
xbar, lam = 1.0, 50.0
for xk in (0.5, 1.0, 1.05, 1.2, 2.0):
    hard = xk**2 if abs(xk) <= xbar else np.inf
    soft = xk**2 + lam * max(0.0, abs(xk) - xbar) ** 2
    print(f"x = {xk:4.2f}:  barrier cost = {hard:>6}   soft cost = {soft:8.4f}")
print("A tiny violation (x = 1.05) makes the barrier cost infinite; the soft penalty grows smoothly.\n")

# ---------------------------------------------------------------------------
# (b) causality: perturb w0 or w_{N-1} and see which terms change
# ---------------------------------------------------------------------------
print("=== (b) which terms does a disturbance influence? (closed loop u = -0.5 x) ===")
N2 = 5


def run(w):
    x = [1.0]
    g = []
    for k in range(N2):
        u = -0.5 * x[k]
        g.append(x[k] ** 2 + u**2)
        x.append(x[k] + u + w[k])
    return np.array(x), np.array(g + [x[-1] ** 2])


x_ref, g_ref = run(np.zeros(N2))
for j, name in ((0, "w_0"), (N2 - 1, "w_{N-1}")):
    w = np.zeros(N2)
    w[j] = 0.3
    x_new, g_new = run(w)
    print(f"perturb {name:8}: states changed  x_k for k in {np.nonzero(~np.isclose(x_new, x_ref))[0].tolist()}, "
          f"cost terms changed g_k for k in {np.nonzero(~np.isclose(g_new, g_ref))[0].tolist()} (index {N2} = g_N)")
