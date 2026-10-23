"""
Optimization-Based Control Methods - Exercise 2, Numerical experiment
Effect of the terminal weight q_N on the finite-horizon LQ policy of Task 2(b).

  x_{k+1} = a x_k + b u_k,  a = 1.2, b = 1.0, x0 = 10, N = 4
  cost   q_N x_4^2 + sum_{k=0}^{3} (x_k^2 + u_k^2),   q_N in {0, 0.1, 1, 10, 100}

  p_4 = q_N,  p_k = 1 + a^2 p_{k+1} / (1 + b^2 p_{k+1})
  u_k = -L_k x_k,  L_k = a b p_{k+1} / (1 + b^2 p_{k+1})

Prints the p- and L-sequences, the trajectories and the costs, and saves
figures/numexp_states.png, figures/numexp_controls.png, figures/numexp_gains.png

Run:  python numerical_experiment.py   [--show]
"""
import sys
from pathlib import Path

import matplotlib
if "--show" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

a, b, x0, N = 1.2, 1.0, 10.0, 4
qNs = [0, 0.1, 1, 10, 100]
FIG = Path(__file__).parent / "figures"
FIG.mkdir(exist_ok=True)


def backward(qN):
    p = np.zeros(N + 1)
    L = np.zeros(N)
    p[N] = qN
    for k in range(N - 1, -1, -1):
        L[k] = a * b * p[k + 1] / (1 + b**2 * p[k + 1])
        p[k] = 1 + a**2 * p[k + 1] / (1 + b**2 * p[k + 1])
    return p, L


def simulate(L):
    x = np.zeros(N + 1)
    u = np.zeros(N)
    x[0] = x0
    for k in range(N):
        u[k] = -L[k] * x[k]
        x[k + 1] = a * x[k] + b * u[k]
    return x, u


results = {}
print(f"{'q_N':>6} | {'p4':>7} {'p3':>7} {'p2':>7} {'p1':>7} {'p0':>7} | "
      f"{'L0':>6} {'L1':>6} {'L2':>6} {'L3':>6} | {'x4':>8} {'u3':>8} | {'J':>9} {'p0*x0^2':>9}")
for qN in qNs:
    p, L = backward(qN)
    x, u = simulate(L)
    J = qN * x[N] ** 2 + np.sum(x[:N] ** 2 + u**2)
    assert abs(J - p[0] * x0**2) < 1e-8
    results[qN] = (p, L, x, u, J)
    print(f"{qN:>6} | " + " ".join(f"{p[k]:7.4f}" for k in range(N, -1, -1)) + " | "
          + " ".join(f"{l:6.4f}" for l in L) + f" | {x[N]:8.4f} {u[N-1]:8.4f} | {J:9.4f} {p[0]*x0**2:9.4f}")

print("\nstate trajectories x_k:")
for qN in qNs:
    print(f"  q_N={qN:>5}: ", np.round(results[qN][2], 4))
print("control trajectories u_k:")
for qN in qNs:
    print(f"  q_N={qN:>5}: ", np.round(results[qN][3], 4))

# limit q_N -> infinity: p_3 -> 1 + a^2/b^2, L_3 -> a/b (deadbeat on the last step)
print(f"\nlimit q_N -> inf: L3 -> a/b = {a/b}, p3 -> 1 + a^2/b^2 = {1 + a**2/b**2}")

colors = plt.cm.viridis(np.linspace(0, 0.9, len(qNs)))
k_x, k_u = np.arange(N + 1), np.arange(N)

fig, ax = plt.subplots(figsize=(7, 4.2))
for c, qN in zip(colors, qNs):
    ax.plot(k_x, results[qN][2], "o-", color=c, label=f"$q_N$ = {qN}")
x_open = x0 * a ** k_x
ax.plot(k_x, x_open, "k:", label="open loop ($u=0$)")
ax.axhline(0, color="grey", lw=0.6)
ins = ax.inset_axes([0.5, 0.2, 0.3, 0.36])
for c, qN in zip(colors, qNs):
    ins.plot(k_x[2:], results[qN][2][2:], "o-", color=c, ms=4)
ins.axhline(0, color="grey", lw=0.6); ins.set_xticks([2, 3, 4]); ins.set_title("zoom $k = 2..4$", fontsize=8)
ins.tick_params(labelsize=7); ins.grid(alpha=0.3)
ax.set_xlabel("stage $k$"); ax.set_ylabel("state $x_k$"); ax.set_xticks(k_x)
ax.set_title("State trajectories for different terminal weights")
ax.legend(loc="upper left", fontsize=8); ax.grid(alpha=0.3); fig.tight_layout()
fig.savefig(FIG / "numexp_states.png", dpi=150)

fig, ax = plt.subplots(figsize=(7, 4.2))
for c, qN in zip(colors, qNs):
    ax.plot(k_u, results[qN][3], "s-", color=c, label=f"$q_N$ = {qN}")
ax.axhline(0, color="grey", lw=0.6)
ax.set_xlabel("stage $k$"); ax.set_ylabel("control $u_k$"); ax.set_xticks(k_u)
ax.set_title("Control trajectories for different terminal weights")
ax.legend(); ax.grid(alpha=0.3); fig.tight_layout()
fig.savefig(FIG / "numexp_controls.png", dpi=150)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.8))
for c, qN in zip(colors, qNs):
    p, L = results[qN][0], results[qN][1]
    a1.plot(k_u, p[:N], "o-", color=c, label=f"$q_N$ = {qN}")
    a2.plot(k_u, L, "s-", color=c, label=f"$q_N$ = {qN}")
a1.set_xlabel("stage $k$"); a1.set_ylabel("$p_k$"); a1.set_xticks(k_u)
a1.set_title("Cost-to-go weights $p_k$, $k=0..3$  ($p_4 = q_N$)")
a2.set_xlabel("stage $k$"); a2.set_ylabel("$L_k$  ($u_k=-L_k x_k$)"); a2.set_xticks(k_u)
a2.axhline(a / b, color="k", ls=":", lw=1); a2.text(0.05, a / b + 0.02, "$a/b$ (deadbeat)", fontsize=8)
a2.set_title("Feedback gains $L_k$")
for axx in (a1, a2):
    axx.grid(alpha=0.3)
a2.legend(fontsize=8); fig.tight_layout()
fig.savefig(FIG / "numexp_gains.png", dpi=150)

print("\nsaved figures to", FIG)
if "--show" in sys.argv:
    plt.show()
