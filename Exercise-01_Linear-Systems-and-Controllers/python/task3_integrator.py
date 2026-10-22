"""
Optimization-Based Control Methods - Exercise 1, Task 3
MIMO PI controller:  u_k = -Kp e_k - KI e^i_k,   e_k = x_k - x_ref,
                     e^i_{k+1} = e^i_k + Ts e_k

Augmented model (Task 3.1):
   [x_{k+1}  ]   [A - B Kp   -B KI] [x_k  ]   [ B Kp  ]
   [e^i_{k+1}] = [  Ts I       I  ] [e^i_k] + [-Ts I  ] x_ref

Task 3.2: asymptotically stable  ->  x_inf = x_ref.
Task 3.3*: only Lyapunov stable  ->  no convergence in general (oscillation / windup).
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

np.set_printoptions(precision=5, suppress=True)
FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)


def augmented(A, B, Kp, KI, Ts):
    n = A.shape[0]
    At = np.block([[A - B @ Kp, -B @ KI], [Ts * np.eye(n), np.eye(n)]])
    Bt = np.vstack([B @ Kp, -Ts * np.eye(n)])
    return At, Bt


def sim_direct(A, B, Kp, KI, Ts, x_ref, N):
    """Simulate the controller exactly as written in the task."""
    n = A.shape[0]
    x, ei = np.zeros(n), np.zeros(n)
    X, EI = [x.copy()], [ei.copy()]
    for _ in range(N):
        e = x - x_ref
        u = -Kp @ e - KI @ ei
        x, ei = A @ x + B @ u, ei + Ts * e
        X.append(x.copy())
        EI.append(ei.copy())
    return np.array(X), np.array(EI)


# ---------------------------------------------------------------- 3.1 check
print("--- 3.1: augmented model gives the same trajectory as the controller itself ---")
A = np.array([[0.9, 0.5], [0.1, 0.3]])
B = np.array([[0.0, 1.0], [1.0, 0.0]])
Kp = np.array([[0.1, 0.2], [0.8, 0.5]])
KI = np.array([[0.0, 0.2], [0.2, 0.0]])
Ts, x_ref, N = 1.0, np.array([1.0, 0.0]), 40
At, Bt = augmented(A, B, Kp, KI, Ts)
print("A_tilde =\n", At, "\nB_tilde =\n", Bt)
X, EI = sim_direct(A, B, Kp, KI, Ts, x_ref, N)
z = np.zeros(4)
Z = [z]
for _ in range(N):
    z = At @ z + Bt @ x_ref
    Z.append(z)
Z = np.array(Z)
print("max difference direct vs augmented:", np.abs(Z - np.hstack([X, EI])).max())

# ---------------------------------------------------------------- 3.2 check
print("\n--- 3.2: asymptotically stable A_tilde -> zero error ---")
print("eig(A_tilde) =", np.linalg.eigvals(At).real, "  all |lambda| < 1:",
      bool(np.all(np.abs(np.linalg.eigvals(At)) < 1)))
z_bar = np.linalg.solve(np.eye(4) - At, Bt @ x_ref)
print("equilibrium z_bar = (I - A_tilde)^-1 B_tilde x_ref =", z_bar,
      "\n  -> x_inf =", z_bar[:2], " e^i_inf =", z_bar[2:])
print("simulated x_40 =", X[-1])

# ---------------------------------------------------------------- 3.3* demos
print("\n--- 3.3*: merely Lyapunov stable A_tilde (1-D plant a = 0.5, b = 1, Ts = 1) ---")
a, b = 0.5, 1.0
cases = {
    "asymptotically stable (kp=0.5, kI=0.25)": (0.5, 0.25),
    "complex pair ON the unit circle (kp=0.5, kI=1)": (0.5, 1.0),
    "lambda = 1, no integral gain (kp=0.5, kI=0)": (0.5, 0.0),
}
N = 40
runs = {}
for label, (kp, kI) in cases.items():
    At1, _ = augmented(np.array([[a]]), np.array([[b]]), np.array([[kp]]), np.array([[kI]]), 1.0)
    ev = np.linalg.eigvals(At1)
    X1, EI1 = sim_direct(np.array([[a]]), np.array([[b]]), np.array([[kp]]),
                         np.array([[kI]]), 1.0, np.array([1.0]), N)
    runs[label] = (X1[:, 0], EI1[:, 0], ev)
    print(f"  {label}\n     eig(A_tilde) = {np.round(ev, 4)}, |lambda| = {np.round(np.abs(ev), 4)}"
          f"\n     x_40 = {X1[-1,0]:+.4f},  e^i_40 = {EI1[-1,0]:+.3f}")

print("\n--- 3.3*: with fewer inputs than states (m < n), lambda = 1 is ALWAYS an eigenvalue ---")
A_c = np.array([[0.9759, 0.3452], [-0.005502, -0.001944]])
B_c = np.array([[0.561], [0.1304]])
rng = np.random.default_rng(1)
for _ in range(3):
    Kp_r, KI_r = rng.normal(size=(1, 2)), rng.normal(size=(1, 2))
    At_r, _ = augmented(A_c, B_c, Kp_r, KI_r, 1.0)
    w = np.linalg.svd(KI_r)[2][-1]           # w with KI w = 0
    v = np.concatenate([np.zeros(2), w])     # eigenvector [0; w]
    print(f"  random Kp, KI: eig = {np.round(np.linalg.eigvals(At_r), 3)};"
          f"  A_tilde [0; w] - [0; w] = {np.round(At_r @ v - v, 12)}")

# ---------------------------------------------------------------- figure
k = np.arange(N + 1)
fig, axs = plt.subplots(1, 2, figsize=(13, 4.3))
cols = ["#2e9e5b", "#e8833a", "#d64545"]
for (label, (x1, ei, ev)), c in zip(runs.items(), cols):
    axs[0].plot(k, x1, "o-", ms=3, color=c, label=label)
    axs[1].plot(k, ei, "o-", ms=3, color=c, label=label)
axs[0].axhline(1, color="k", ls="--", lw=1)
axs[0].set_title("state x_k (reference 1)")
axs[1].set_title("integrator state eⁱ_k")
for a_ in axs:
    a_.set_xlabel("step k")
    a_.grid(alpha=0.3)
axs[0].legend(fontsize=8)
fig.suptitle("Task 3.3* — only an asymptotically stable augmented matrix guarantees zero error",
             fontsize=12, weight="bold")
fig.tight_layout()
fig.savefig(FIG / "task3_lyapunov_cases.png", dpi=110)
print("\nsaved figures/task3_lyapunov_cases.png")
