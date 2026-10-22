"""
Experiment 3 — Controllability / reachability: where can the motors take the drone?

1. Rank of W_c(N) = [B, AB, ..., A^{N-1}B] as N grows, for
      - the healthy drone        (inputs: thrust dF and torque tau)
      - torque actuator failed   (only thrust dF)
      - thrust control failed    (only torque tau)
2. PBH test for each case.
3. Steering the healthy drone from hover at the origin to a target in exactly N steps
   with the minimum-energy input  u = pinv(W_c(N)) (x_f - A^N x_0),
   for several N: the fewer steps, the more violent the input.
4. The broken drone (thrust only) tries to reach the same target: it can only reach
   the part of the target that lies in its reachable subspace.

Run:  python exp3_controllability.py   (add --show to open windows)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import orth
from lab_common import (planar_quadcopter, ctrb, rank, pbh_controllability, simulate,
                        draw_quad, save, save_anim, header,
                        C_BLUE, C_ORANGE, C_RED, C_GREEN, C_PURPLE, C_GREY)

header("Experiment 3: controllability")
d = planar_quadcopter()
A, B, Ts = d["A"], d["B"], d["Ts"]
n = A.shape[0]
names = d["states"]

variants = {
    "healthy drone (thrust + torque)": B,
    "torque failed (thrust only)":     B[:, [0]],
    "thrust control failed (torque only)": B[:, [1]],
}
colors = [C_GREEN, C_RED, C_PURPLE]


def min_energy_steer(A, B, x0, xf, N):
    """Minimum-norm input sequence reaching xf in N steps (exactly, if reachable)."""
    m = B.shape[1]
    W = ctrb(A, B, N)                        # columns: B, AB, ..., A^{N-1}B
    ubar = np.linalg.pinv(W) @ (xf - np.linalg.matrix_power(A, N) @ x0)
    # ubar = [u_{N-1}; u_{N-2}; ...; u_0]  -> reverse the blocks
    U = ubar.reshape(N, m)[::-1]
    return U


# ----------------------------------------------------------------------------
# 1 + 2: rank growth and PBH
# ----------------------------------------------------------------------------
Ns = np.arange(1, 9)
ranks = {}
for name, Bv in variants.items():
    ranks[name] = [rank(ctrb(A, Bv, N)) for N in Ns]
    Wc = ctrb(A, Bv)
    print(f"\n{name}:  rank W_c = {rank(Wc)} of {n}")
    for lam, r in pbh_controllability(A, Bv):
        print(f"    PBH  lambda = {lam.real:.4f}: rank[lambda I - A, B] = {r}"
              f"  -> {'controllable' if r == n else 'UNCONTROLLABLE'} mode")
    R = orth(Wc)
    basis_txt = ", ".join(names[i] for i in range(n) if np.linalg.norm(R[i]) > 1e-6)
    print(f"    reachable subspace has dimension {R.shape[1]}; states it involves: {basis_txt}")

# ----------------------------------------------------------------------------
# 3: steer the healthy drone in N steps
# ----------------------------------------------------------------------------
x0 = np.zeros(n)
xf = np.array([2.0, 1.0, 0, 0, 0, 0])        # 2 m right, 1 m up, level, at rest
N_list = [4, 10, 30]
runs = {}
for N in N_list:
    U = min_energy_steer(A, B, x0, xf, N)
    X, _ = simulate(A, B, x0, N, u=U)
    runs[N] = (X, U)
    energy = np.sum(U ** 2)
    print(f"\n  N = {N:2d} steps ({N*Ts:.1f} s): final state {X[-1]}, "
          f"max |tau| = {np.abs(U[:,1]).max():8.3f} N m, energy sum(u^2) = {energy:10.3f}")

N_energy = np.arange(4, 61, 2)
energies = []
for N in N_energy:
    U = min_energy_steer(A, B, x0, xf, N)
    energies.append(np.sum(U ** 2))

# ----------------------------------------------------------------------------
# 4: the broken drone tries
# ----------------------------------------------------------------------------
Bb = B[:, [0]]
Ub = min_energy_steer(A, Bb, x0, xf, 30)
Xb, _ = simulate(A, Bb, x0, 30, u=Ub)
print(f"\n  Broken drone (thrust only), best attempt at x_f = {xf}:\n    reaches {Xb[-1]}"
      "\n    -> altitude is fine, the sideways part of the target is unreachable")

# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(14, 9.5))

ax = axs[0, 0]
for (name, r), col, mk in zip(ranks.items(), colors, "osd"):
    ax.plot(Ns, r, mk + "-", color=col, lw=2, ms=8, label=name)
ax.axhline(n, color=C_GREY, ls="--")
ax.text(8, n + 0.15, "n = 6 (everything reachable)", ha="right", fontsize=9, color=C_GREY)
ax.set_xlabel("number of steps N")
ax.set_ylabel("rank [B, AB, …, A^(N−1)B]")
ax.set_title("How many directions can N steps of input reach?")
ax.set_ylim(0, 7)
ax.legend(fontsize=9, loc="lower right")

ax = axs[0, 1]
for N, col in zip(N_list, (C_RED, C_ORANGE, C_BLUE)):
    X, U = runs[N]
    ax.plot(X[:, 0], X[:, 1], "o-", ms=3, color=col, label=f"N = {N} steps ({N*Ts:.1f} s)")
ax.plot(Xb[:, 0], Xb[:, 1], "s-", ms=3, color=C_PURPLE, label="broken drone (thrust only)")
for X, col in [(runs[30][0], C_BLUE)]:
    for k in range(0, 31, 6):
        draw_quad(ax, X[k, 0], X[k, 1], X[k, 2], color=col, alpha=0.35, arm=0.15)
draw_quad(ax, 0, 0, 0, color=C_GREY, arm=0.15)
ax.scatter([xf[0]], [xf[1]], marker="*", s=300, color=C_GREEN, zorder=6, label="target")
ax.set_xlabel("y [m]")
ax.set_ylabel("z [m]")
ax.set_aspect("equal")
ax.set_xlim(-0.6, 3.0)
ax.set_ylim(-0.5, 1.6)
ax.set_title("Steering from hover to the target in exactly N steps")
ax.legend(fontsize=8, loc="upper left")

ax = axs[1, 0]
for N, col in zip(N_list, (C_RED, C_ORANGE, C_BLUE)):
    X, U = runs[N]
    tt = np.arange(N) * Ts
    ax.step(tt, U[:, 1], where="post", color=col, lw=2, label=f"torque τ, N = {N}")
ax.axhline(0, color=C_GREY, lw=0.6)
ax.set_yscale("symlog", linthresh=0.01)
ax.set_xlabel("time [s]")
ax.set_ylabel("torque τ [N m]  (symlog scale)")
ax.set_title("The input sequences (torque shown)")
ax.legend(fontsize=9)

ax = axs[1, 1]
ax.semilogy(N_energy * Ts, energies, "o-", color=C_BLUE)
ax.set_xlabel("time allowed  N·Ts  [s]")
ax.set_ylabel("input energy  Σ‖u_k‖²  (log scale)")
ax.set_title("Reachable ≠ cheap: rushing costs huge inputs")
ax.grid(alpha=0.3, which="both")

fig.suptitle("Experiment 3 — controllability: what the motors can (and cannot) reach",
             fontsize=13, weight="bold")
fig.tight_layout()
save(fig, "exp3_controllability.png")

# ----------------------------------------------------------------------------
# Animation: healthy (N = 30) vs broken drone
# ----------------------------------------------------------------------------
Xh = runs[30][0]
fig, ax = plt.subplots(figsize=(8, 4.6))


def update(k):
    ax.clear()
    ax.scatter([xf[0]], [xf[1]], marker="*", s=300, color=C_GREEN, zorder=6)
    ax.plot(Xh[:k + 1, 0], Xh[:k + 1, 1], ":", color=C_BLUE)
    ax.plot(Xb[:k + 1, 0], Xb[:k + 1, 1], ":", color=C_PURPLE)
    draw_quad(ax, Xh[k, 0], Xh[k, 1], Xh[k, 2], color=C_BLUE)
    draw_quad(ax, Xb[k, 0], Xb[k, 1], Xb[k, 2], color=C_PURPLE)
    ax.text(Xh[k, 0], Xh[k, 1] + 0.25, "healthy", color=C_BLUE, ha="center", fontsize=9)
    ax.text(Xb[k, 0], Xb[k, 1] - 0.3, "torque motor failed", color=C_PURPLE, ha="center",
            fontsize=9)
    ax.set_xlim(-0.8, 2.8)
    ax.set_ylim(-0.6, 1.5)
    ax.set_aspect("equal")
    ax.set_xlabel("y [m]")
    ax.set_ylabel("z [m]")
    ax.set_title(f"Steering to ★ with the minimum-energy input,  t = {k*Ts:.1f} s\n"
                 "rank W_c = 6 (healthy) vs rank W_c = 2 (thrust only)", fontsize=10)


save_anim(fig, update, list(range(0, 31)) + [30] * 10, "exp3_steering.gif", fps=10)
