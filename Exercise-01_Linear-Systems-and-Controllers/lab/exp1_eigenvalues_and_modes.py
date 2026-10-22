"""
Experiment 1 — What an eigenvalue DOES to a discrete-time system.

For a single mode, x_k = lambda^k x_0. We place eight different eigenvalues in the
complex plane and plot the resulting sequences side by side, so you can see:
    |lambda| < 1  -> decays        |lambda| = 1 -> stays the same size
    |lambda| > 1  -> grows         negative / complex lambda -> flips sign / oscillates

Run:  python exp1_eigenvalues_and_modes.py        (add --show to open windows)
"""
import numpy as np
import matplotlib.pyplot as plt
from lab_common import (unit_circle, save, header, C_BLUE, C_GREEN, C_RED,
                        C_ORANGE, C_PURPLE, C_GREY)

header("Experiment 1: eigenvalues and modes")

cases = [
    ("A", 0.5,                       "λ = 0.5: fast decay"),
    ("B", 0.9,                       "λ = 0.9: slow decay"),
    ("C", -0.7,                      "λ = −0.7: decay, sign flips"),
    ("D", 0.85 * np.exp(0.5j),       "λ = 0.85·e^{±0.5i}: damped oscillation"),
    ("E", 1.0,                       "λ = 1: stays constant"),
    ("F", np.exp(1.2j),              "|λ| = 1, complex: endless oscillation"),
    ("G", 1.08,                      "λ = 1.08: grows"),
    ("H", 1.05 * np.exp(2.2j),       "|λ| = 1.05, complex: growing oscillation"),
]


def colour(lam):
    r = abs(lam)
    return C_GREEN if r < 1 - 1e-9 else (C_ORANGE if abs(r - 1) < 1e-9 else C_RED)


fig = plt.figure(figsize=(15, 8.2))
gs = fig.add_gridspec(4, 6, width_ratios=[1.25, 1.25, 1, 1, 1, 1], hspace=0.75, wspace=0.45)

# --- left: eigenvalue map --------------------------------------------------
ax0 = fig.add_subplot(gs[:, :2])
unit_circle(ax0, "Where each eigenvalue sits")
for tag, lam, _ in cases:
    pts = [lam] if abs(np.imag(lam)) < 1e-12 else [lam, np.conj(lam)]
    for p in pts:
        ax0.scatter(p.real, p.imag, s=110, marker="x", lw=2.5, color=colour(lam), zorder=5)
    off = {"B": (-6, -20), "E": (-4, 10), "G": (4, -20)}.get(tag, (6, 6))
    ax0.annotate(tag, (np.real(lam), abs(np.imag(lam))), xytext=off,
                 textcoords="offset points", weight="bold", fontsize=12, color=colour(lam))
ax0.set_xlim(-1.35, 1.35)
ax0.set_ylim(-1.35, 1.35)
ax0.text(-1.3, 1.22, "green: |λ|<1 (decays)\norange: |λ|=1 (persists)\nred: |λ|>1 (grows)",
         fontsize=9, va="top")

# --- right: the sequence lambda^k for each case ------------------------------
k = np.arange(0, 31)
for i, (tag, lam, title) in enumerate(cases):
    ax = fig.add_subplot(gs[i // 2, 2 + 2 * (i % 2): 4 + 2 * (i % 2)])
    seq = np.real(lam ** k)                 # real part = response of the real 2x2 system
    ax.stem(k, seq, linefmt=colour(lam), markerfmt="o", basefmt=" ")
    ax.axhline(0, color=C_GREY, lw=0.6)
    ax.set_title(f"{tag}: {title}", fontsize=9.5, color=colour(lam))
    ax.set_xlabel("step k", fontsize=8)
    ax.tick_params(labelsize=8)
    lim = max(1.1, np.max(np.abs(seq)) * 1.1)
    ax.set_ylim(-lim, lim)
    print(f"  {tag}: |lambda| = {abs(lam):.3f}  ->  lambda^30 has size {abs(lam**30):.3g}")

fig.suptitle("Experiment 1 — a mode behaves like λᵏ: the distance from the origin decides "
             "decay or growth, the angle decides oscillation", fontsize=13, weight="bold")
save(fig, "exp1_eigenvalue_modes.png")

# ----------------------------------------------------------------------------
# Second figure: a 2-state system = sum of its modes (eigenvector picture)
# ----------------------------------------------------------------------------
A = np.array([[0.5, 0.2], [0.1, 0.4]])        # practice problem 5: eigenvalues 0.6, 0.3
lam, V = np.linalg.eig(A)
order = np.argsort(-lam)
lam, V = lam[order].real, V[:, order].real
print(f"\n  Practice problem 5 matrix: eigenvalues {lam}, eigenvectors (columns)\n{V}")

x0 = np.array([0.2, 1.0])
c = np.linalg.solve(V, x0)                     # x0 = c1 v1 + c2 v2
N = 12
X = np.array([np.linalg.matrix_power(A, j) @ x0 for j in range(N + 1)])
M1 = np.array([c[0] * lam[0] ** j * V[:, 0] for j in range(N + 1)])
M2 = np.array([c[1] * lam[1] ** j * V[:, 1] for j in range(N + 1)])

fig, axs = plt.subplots(1, 2, figsize=(13, 5))
ax = axs[0]
for v, l, col in zip(V.T, lam, (C_BLUE, C_PURPLE)):
    ax.plot([-1.3 * v[0], 1.3 * v[0]], [-1.3 * v[1], 1.3 * v[1]], "--", color=col, lw=1.2,
            label=f"eigenvector for λ = {l:.1f}")
ax.plot(X[:, 0], X[:, 1], "o-", color=C_ORANGE, label="state x_k")
for j in range(3):
    ax.annotate(f"x{j}", X[j], xytext=(6, -12), textcoords="offset points", fontsize=9)
ax.set_aspect("equal")
ax.axhline(0, color=C_GREY, lw=0.5)
ax.axvline(0, color=C_GREY, lw=0.5)
ax.set_xlabel("x₁")
ax.set_ylabel("x₂")
ax.set_title("State plane: x_k slides in along the slow eigenvector")
ax.legend(fontsize=9)

ax = axs[1]
kk = np.arange(N + 1)
ax.plot(kk, X[:, 0], "o-", color=C_ORANGE, label="x₁ (total)")
ax.plot(kk, M1[:, 0], "s--", color=C_BLUE, label="part from λ = 0.6 mode")
ax.plot(kk, M2[:, 0], "^--", color=C_PURPLE, label="part from λ = 0.3 mode")
ax.axhline(0, color=C_GREY, lw=0.5)
ax.set_xlabel("step k")
ax.set_title("x_k = c₁·0.6ᵏ·v₁ + c₂·0.3ᵏ·v₂  (the fast mode dies first)")
ax.legend(fontsize=9)
fig.suptitle("Experiment 1b — a system is a sum of modes", fontsize=13, weight="bold")
save(fig, "exp1b_sum_of_modes.png")
