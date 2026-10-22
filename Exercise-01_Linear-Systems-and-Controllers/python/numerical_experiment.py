"""
Optimization-Based Control Methods - Exercise 1, Numerical Experiment
Tracking control: comparison of the controllers from Task 1, Task 2 and Task 3

    Task 1:  u_k = -K (x_k - x_ref)                       (P state feedback)
    Task 2:  u_k = -K (x_k - F x_ref)                     (P + reference prefilter)
    Task 3:  u_k = -Kp e_k - KI e^i_k,  e^i_{k+1} = e^i_k + Ts e_k     (MIMO PI)

Systems a), b) (two inputs) and c) DC motor (one input); step reference
x_ref = [1, 0], initial state x_0 = [0, 0].

Assumption: the sampling time Ts of the systems is not given on the sheet; we use
Ts = 1. Ts only scales the integrator (KI*Ts is what matters), so no conclusion
depends on it.

Run:  python numerical_experiment.py        -> prints every step, saves figures/
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import place_poles

np.set_printoptions(precision=5, suppress=True)
FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)

Ts = 1.0
I2 = np.eye(2)
x_ref = np.array([1.0, 0.0])
x0 = np.zeros(2)
N = 40

SYSTEMS = {
    "a": (np.array([[0.9, 0.5], [0.1, 0.3]]), np.array([[0.0, 1.0], [1.0, 0.0]])),
    "b": (np.array([[1.0, 0.5], [0.1, 0.7]]), np.array([[0.0, 1.0], [1.0, 0.0]])),
    "c": (np.array([[0.9759, 0.3452], [-0.005502, -0.001944]]), np.array([[0.561], [0.1304]])),
}


def title(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


def rank(M):
    return np.linalg.matrix_rank(M, tol=1e-9)


# =============================================================================
# Q1: single- or multiple-input?
# =============================================================================
title("Q1  Number of inputs m = number of columns of B")
for name, (A, B) in SYSTEMS.items():
    n, m = B.shape
    print(f"  system {name}): B is {n}x{m}  ->  m = {m}  "
          f"({'multiple-input' if m > 1 else 'single-input'}),  rank B = {rank(B)}")

# =============================================================================
# Q2: can the controllers of Tasks 1-3 reach zero tracking error?
# =============================================================================
title("Q2  Conditions for zero asymptotic tracking error (x_ref = [1, 0])")
for name, (A, B) in SYSTEMS.items():
    n, m = B.shape
    print(f"\n  system {name}):  eig(A) = {np.linalg.eigvals(A).real}")
    # Task 1: A x_ref = x_ref ?
    print(f"    Task 1: A x_ref = {A @ x_ref}  ==  x_ref ? {np.allclose(A @ x_ref, x_ref)}")
    # Is x_ref an equilibrium that SOME constant input can hold?  (I - A) x_ref in im(B)?
    v = (I2 - A) @ x_ref
    feasible = rank(np.column_stack([B, v])) == rank(B)
    print(f"    (I - A) x_ref = {v};  in im(B)?  {feasible}   "
          f"(x_ref can be held by a constant input iff this is true)")
    # Task 2 general formula needs BK invertible -> rank B = n
    print(f"    Task 2: F = (BK)^-1 (I - A) + I needs BK invertible -> rank B = n ? {rank(B) == n}")
    # Task 3: augmented system stabilizable?  needs rank [A - I, B; Ts I, 0] = 2n -> m >= n
    Aa = np.block([[A, np.zeros((n, n))], [Ts * np.eye(n), np.eye(n)]])
    Ba = np.vstack([B, np.zeros((n, m))])
    Wa = np.hstack([np.linalg.matrix_power(Aa, i) @ Ba for i in range(2 * n)])
    pbh1 = rank(np.hstack([np.eye(2 * n) - Aa, Ba]))
    print(f"    Task 3: augmented (A_aug, B_aug): rank W_c = {rank(Wa)} of {2*n};  "
          f"PBH at lambda=1: rank = {pbh1} of {2*n}  -> "
          f"{'stabilizable' if pbh1 == 2*n else 'lambda = 1 can NOT be moved'}")

# best achievable steady state for c)
A, B = SYSTEMS["c"]
d = np.linalg.solve(I2 - A, B).ravel()           # all equilibria: x = d * u_bar
print(f"\n  system c): every equilibrium has the form x = (I - A)^-1 B u_bar = u_bar * {d}")
print(f"            the x-direction is {d / d[0]}  ->  x2 = {d[1]/d[0]:.3e} * x1")
print(f"            so x_ref = [1, 0] is NOT exactly holdable; closest: [1, {d[1]/d[0]:.2e}]"
      f" with u_bar = {1/d[0]:.5f}")

# =============================================================================
# Q3: designs
# =============================================================================
title("Q3  Controller design")
designs = {}

# ---- systems a) and b): B is invertible (B^-1 = B), design by hand ------------
for name in "ab":
    A, B = SYSTEMS[name]
    Binv = np.linalg.inv(B)
    D = 0.5 * I2                                # desired closed loop A - BK = 0.5 I
    K = Binv @ (A - D)
    F = np.linalg.inv(B @ K) @ (I2 - A) + I2
    alpha, beta = 0.1, 0.2                      # A - B Kp = alpha I,  B KI = beta I
    Kp = Binv @ (A - alpha * I2)
    KI = beta * Binv
    At = np.block([[A - B @ Kp, -B @ KI], [Ts * I2, I2]])
    designs[name] = dict(K=K, F=F, Kp=Kp, KI=KI)
    print(f"\n  system {name}):")
    print(f"    K  = B^-1 (A - 0.5 I) =\n{K}\n    A - BK =\n{A - B @ K}"
          f"   eig = {np.linalg.eigvals(A - B @ K).real}")
    print(f"    F  = (BK)^-1 (I - A) + I =\n{F}")
    print(f"    Kp = B^-1 (A - 0.1 I) =\n{Kp}\n    KI = 0.2 B^-1 =\n{KI}")
    print(f"    eig(A_tilde) = {np.sort(np.linalg.eigvals(At).real)}")

# ---- system c): one input -----------------------------------------------------
A, B = SYSTEMS["c"]
a11, a12, a21, a22 = A.ravel()
b1, b2 = B.ravel()
# tr(A - BK) = tr A - b1 k1 - b2 k2 ;  det(A - BK) = det A - k1 (b1 a22 - b2 a12) - k2 (b2 a11 - b1 a21)
Mlin = np.array([[b1, b2], [b1 * a22 - b2 * a12, b2 * a11 - b1 * a21]])
rhs = np.array([np.trace(A) - (0.5 + 0.6), np.linalg.det(A) - 0.5 * 0.6])
K = np.linalg.solve(Mlin, rhs).reshape(1, 2)
Mcl = A - B @ K
g = np.linalg.solve(I2 - Mcl, B).ravel()        # x_inf = g * (K F x_ref)
f = (x_ref[0] / g[0]) / (K @ x_ref).item()        # F = f I makes x1 exact
F = f * I2
# PI: integrate x1-error only (KI = [kI, 0]); design the 3-state x1-integrator loop
Ao = np.block([[A, np.zeros((2, 1))], [Ts * np.array([[1.0, 0.0]]), np.eye(1)]])
Bo = np.vstack([B, [[0.0]]])
G = place_poles(Ao, Bo, [0.5, 0.6, 0.7]).gain_matrix
Kp = G[:, :2]
KI = np.array([[G[0, 2], 0.0]])
At = np.block([[A - B @ Kp, -B @ KI], [Ts * I2, I2]])
designs["c"] = dict(K=K, F=F, Kp=Kp, KI=KI)
print("\n  system c):")
print(f"    pole placement equations  [b1 b2; b1a22-b2a12  b2a11-b1a21] k = [trA-1.1; detA-0.3]:")
print(f"      {Mlin[0]} k = {rhs[0]:.5f}\n      {Mlin[1]} k = {rhs[1]:.5f}")
print(f"    K = {K}   eig(A - BK) = {np.linalg.eigvals(Mcl).real}")
print(f"    F = f I with f = {f:.5f}  (puts x1 exactly on 1; x2 cannot be fixed)")
print(f"    Kp = {Kp},  KI = {KI}   (integrator acts on the x1-error only)")
print(f"    eig(A_tilde) = {np.sort(np.linalg.eigvals(At).real)}   <- 1 is always there (m < n)")

# =============================================================================
# Q4: simulate and compare
# =============================================================================
title("Q4  Simulation: step x_ref = [1, 0], x_0 = [0, 0]")


def run(A, B, ctrl, des, x_ref=x_ref, N=N):
    n = A.shape[0]
    X = np.zeros((N + 1, n))
    EI = np.zeros((N + 1, n))
    U = np.zeros((N, B.shape[1]))
    X[0] = x0
    for k in range(N):
        e = X[k] - x_ref
        if ctrl == "P":
            U[k] = -des["K"] @ e
        elif ctrl == "P+F":
            U[k] = -des["K"] @ (X[k] - des["F"] @ x_ref)
        else:
            U[k] = -des["Kp"] @ e - des["KI"] @ EI[k]
        X[k + 1] = A @ X[k] + B @ U[k]
        EI[k + 1] = EI[k] + Ts * e
    return X, U, EI


CTRL = [("P", "Task 1: P"), ("P+F", "Task 2: P + prefilter F"), ("PI", "Task 3: PI")]
COL = {"P": "#d64545", "P+F": "#e8833a", "PI": "#2a6fdb"}
results = {}
for name, (A, B) in SYSTEMS.items():
    for c, _ in CTRL:
        X, U, EI = run(A, B, c, designs[name])
        results[(name, c)] = (X, U, EI)
        e = X[-1] - x_ref
        print(f"  system {name}) {c:4s}: x_{N} = {X[-1]},  error = {e},  |e| = {np.linalg.norm(e):.2e}")

# ---- model mismatch: the true plant is 3 % different -------------------------
title("Extra: robustness — true A is 3 % larger than the model used for the design")
mismatch = {}
for name, (A, B) in SYSTEMS.items():
    At_true = 1.03 * A
    for c, _ in CTRL:
        X, U, EI = run(At_true, B, c, designs[name], N=60)
        mismatch[(name, c)] = X
        print(f"  system {name}) {c:4s}: final error = {X[-1] - x_ref}")

# ---- windup in system c) with an unholdable reference ------------------------
title("Extra: system c) with x_ref = [1, 0.1] (clearly not holdable) -> integrator windup")
x_ref_bad = np.array([1.0, 0.1])
Xw, Uw, EIw = run(*SYSTEMS["c"], "PI", designs["c"], x_ref=x_ref_bad, N=60)
print(f"  x_60 = {Xw[-1]}   e^i_60 = {EIw[-1]}  (second component grows by ~0.1 every step)")

# =============================================================================
# Figures
# =============================================================================
k = np.arange(N + 1)
fig, axs = plt.subplots(3, 3, figsize=(15, 11))
for r, name in enumerate("abc"):
    for c, lab in CTRL:
        X, U, EI = results[(name, c)]
        err = np.linalg.norm(X - x_ref, axis=1)
        axs[r, 0].plot(k, X[:, 0], "o-", ms=3, color=COL[c], label=lab)
        axs[r, 1].plot(k, X[:, 1], "o-", ms=3, color=COL[c], label=lab)
        axs[r, 2].semilogy(k, np.maximum(err, 1e-16), "o-", ms=3, color=COL[c], label=lab)
    axs[r, 0].axhline(1, color="k", ls="--", lw=1)
    axs[r, 1].axhline(0, color="k", ls="--", lw=1)
    axs[r, 0].set_ylabel(f"system {name})\n\nx₁", fontsize=11)
    axs[r, 1].set_ylabel("x₂")
    axs[r, 2].set_ylabel("‖x_k − x_ref‖ (log)")
    axs[r, 2].set_ylim(1e-16, 10)
    for a in axs[r]:
        a.grid(alpha=0.3)
axs[0, 0].set_title("first state x₁ (reference 1)")
axs[0, 1].set_title("second state x₂ (reference 0)")
axs[0, 2].set_title("tracking error norm")
axs[0, 2].legend(fontsize=9)
for a in axs[2]:
    a.set_xlabel("step k")
fig.suptitle("Numerical experiment — step response of the three tracking controllers "
             "(x_ref = [1, 0], x₀ = 0)", fontsize=13, weight="bold")
fig.tight_layout()
fig.savefig(FIG / "numexp_step_responses.png", dpi=110)
plt.close(fig)

fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
k60 = np.arange(61)
for a, name in zip(axs, "abc"):
    for c, lab in CTRL:
        X = mismatch[(name, c)]
        a.semilogy(k60, np.maximum(np.linalg.norm(X - x_ref, axis=1), 1e-16), color=COL[c],
                   lw=2, label=lab)
    a.set_title(f"system {name}) with a 3 % model error")
    a.set_xlabel("step k")
    a.set_ylim(1e-12, 10)
    a.grid(alpha=0.3)
axs[0].set_ylabel("‖x_k − x_ref‖ (log)")
axs[0].legend(fontsize=9)
fig.suptitle("Robustness — the prefilter relies on the exact model, the integrator does not",
             fontsize=13, weight="bold")
fig.tight_layout()
fig.savefig(FIG / "numexp_model_mismatch.png", dpi=110)
plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(13, 4.3))
Xc, Uc, EIc = results[("c", "PI")]
axs[0].plot(k60, Xw[:, 0], color="#2a6fdb", lw=2, label="x₁ (ref 1)")
axs[0].plot(k60, Xw[:, 1], color="#2e9e5b", lw=2, label="x₂ (ref 0.1)")
axs[0].axhline(1, color="#2a6fdb", ls="--", lw=1)
axs[0].axhline(0.1, color="#2e9e5b", ls="--", lw=1)
axs[0].set_title("system c), PI, x_ref = [1, 0.1]: x₁ exact, x₂ never reaches 0.1")
axs[0].set_xlabel("step k")
axs[0].legend()
axs[1].plot(k60, EIw[:, 0], color="#2a6fdb", lw=2, label="eⁱ₁ (acts through K_I)")
axs[1].plot(k60, EIw[:, 1], color="#d64545", lw=2, label="eⁱ₂ (in ker K_I: grows forever)")
axs[1].set_title("integrator states: the λ = 1 mode winds up")
axs[1].set_xlabel("step k")
axs[1].legend()
for a in axs:
    a.grid(alpha=0.3)
fig.suptitle("System c) — one input, two states: eigenvalue 1 of the augmented matrix "
             "cannot be moved (Task 3.3*)", fontsize=12, weight="bold")
fig.tight_layout()
fig.savefig(FIG / "numexp_dc_motor_windup.png", dpi=110)
plt.close(fig)
print("\nFigures saved in", FIG)
