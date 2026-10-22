"""
Experiment 2 — The drone WITHOUT a controller (open loop).

We build the linearised planar quadcopter, look at the eigenvalues of A, and fly it
with zero input (u = 0, i.e. exactly hover thrust and no torque) from three
slightly disturbed starting states. Concepts: eigenvalues, stability,
lambda_discrete = e^{lambda_continuous * Ts}, repeated eigenvalue 1 (Jordan block).

Run:  python exp2_drone_open_loop.py   (add --show to open windows)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import null_space
from lab_common import (planar_quadcopter, simulate, unit_circle, plot_eigs, draw_quad,
                        save, save_anim, header, C_BLUE, C_ORANGE, C_RED, C_GREEN, C_GREY)

header("Experiment 2: the drone in open loop")
sys_ = planar_quadcopter()
A, B, Ts = sys_["A"], sys_["B"], sys_["Ts"]
print("Discrete A (Ts = 0.1 s):\n", A)
print("Discrete B:\n", B)

eig_c = np.linalg.eigvals(sys_["Ac"])
eig_d = np.linalg.eigvals(A)
print("\nContinuous-time eigenvalues  :", np.sort(eig_c.real))
print("Discrete-time eigenvalues    :", np.sort(eig_d.real))
print("Check lambda_d = e^{lambda_c Ts}:", np.sort(np.exp(eig_c.real * Ts)))

# algebraic vs geometric multiplicity of lambda = 1
n = A.shape[0]
alg = int(np.sum(np.isclose(eig_d, 1.0)))
geo = null_space(np.eye(n) - A, rcond=1e-9).shape[1]
print(f"\nlambda = 1: algebraic multiplicity {alg}, geometric multiplicity {geo}"
      f"  ->  only {geo} independent eigenvectors, so some states GROW like k, k^2, ...")
print("Asymptotically stable?", bool(np.all(np.abs(eig_d) < 1)))

# ----------------------------------------------------------------------------
# Three open-loop flights, u = 0
# ----------------------------------------------------------------------------
N = 50                       # 5 seconds
t = np.arange(N + 1) * Ts
cases = {
    "a) tilted 3° at start":        np.array([0, 0, np.deg2rad(3), 0, 0, 0]),
    "b) small pitch-rate kick":     np.array([0, 0, 0, 0, 0, np.deg2rad(2)]),
    "c) small upward velocity":     np.array([0, 0, 0, 0, 0.3, 0]),
}
cols = [C_BLUE, C_RED, C_GREEN]
sims = {name: simulate(A, B, x0, N)[0] for name, x0 in cases.items()}
for name, X in sims.items():
    print(f"  {name:28s} after 5 s: y = {X[-1,0]:+7.2f} m, z = {X[-1,1]:+5.2f} m, "
          f"theta = {np.rad2deg(X[-1,2]):+6.1f}°")

fig = plt.figure(figsize=(15, 8))
gs = fig.add_gridspec(3, 3, width_ratios=[1.25, 1, 1], hspace=0.55, wspace=0.45)

ax = fig.add_subplot(gs[:2, 0])
unit_circle(ax, "Open-loop eigenvalues of the drone")
plot_eigs(ax, eig_d, color=C_RED, label="eigenvalues of A", annotate=False)
ax.set_xlim(-1.3, 1.3)
ax.set_ylim(-1.3, 1.3)
ax.legend(loc="lower left", fontsize=9)
axz = ax.inset_axes([0.58, 0.62, 0.38, 0.33])
t_ = np.linspace(-0.25, 0.25, 50)
axz.plot(np.cos(t_), np.sin(t_), color=C_GREEN)
plot_eigs(axz, eig_d, color=C_RED, size=60, jitter=0.004)
axz.set_xlim(0.95, 1.02)
axz.set_ylim(-0.04, 0.04)
axz.set_title("zoom near 1", fontsize=8)
axz.tick_params(labelsize=7)

axt = fig.add_subplot(gs[2, 0])
axt.axis("off")
axt.text(-0.05, 1.05, "What the eigenvalues say\n"
         "• 4 × λ = 1  (y, z, θ, ω):\n"
         "    integrators, e^(0·Ts) = 1\n"
         f"• 2 × λ = e^(−c·Ts) = {np.exp(-0.2*Ts):.3f}\n"
         "    velocities, damped by drag\n"
         f"• λ = 1: {alg} repeats, only {geo} eigenvectors\n"
         "    → growth like k, k², k³\n"
         "• Not asymptotically stable:\n"
         "    hovering needs a controller", va="top", fontsize=9.5, family="monospace")

labels = [("y  horizontal position [m]", 0, 1), ("θ  tilt [deg]", 2, 180 / np.pi),
          ("z  altitude [m]", 1, 1)]
for r, (lab, idx, scale) in enumerate(labels):
    a = fig.add_subplot(gs[r, 1:])
    for (name, X), col in zip(sims.items(), cols):
        a.plot(t, X[:, idx] * scale, color=col, lw=2, label=name)
    a.axhline(0, color=C_GREY, lw=0.6)
    a.set_ylabel(lab, fontsize=9)
    if r == 0:
        a.legend(fontsize=9, loc="lower left")
    if r == 2:
        a.set_xlabel("time [s]")
fig.suptitle("Experiment 2 — without feedback, tiny disturbances never die out: "
             "the drone drifts away", fontsize=13, weight="bold")
save(fig, "exp2_drone_open_loop.png")

# ----------------------------------------------------------------------------
# Animation of case a)
# ----------------------------------------------------------------------------
X = sims["a) tilted 3° at start"]
fig, ax = plt.subplots(figsize=(8, 4))


def update(k):
    ax.clear()
    ax.plot(X[:k + 1, 0], X[:k + 1, 1], ":", color=C_GREY)
    draw_quad(ax, 0, 0, 0, color=C_GREEN, alpha=0.25)
    draw_quad(ax, X[k, 0], X[k, 1], X[k, 2] * 4, color=C_BLUE, thrust=0.4)  # tilt drawn ×4
    ax.set_xlim(X[:, 0].min() - 0.6, 0.8)
    ax.set_ylim(-1, 1)
    ax.set_aspect("equal")
    ax.set_xlabel("y [m]")
    ax.set_ylabel("z [m]")
    ax.set_title(f"Open loop, started 3° tilted   t = {k*Ts:.1f} s   y = {X[k,0]:+.2f} m\n"
                 "(green ghost = hover target; tilt drawn ×4 so you can see it)", fontsize=10)


save_anim(fig, update, range(0, N + 1), "exp2_drone_drift.gif", fps=12)
