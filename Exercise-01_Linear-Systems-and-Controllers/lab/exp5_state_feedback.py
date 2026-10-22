"""
Experiment 5 — State feedback u = -K x: moving the eigenvalues.

The drone is reachable (Experiment 3), so K can place the eigenvalues of A - BK anywhere
(pole placement). We design three controllers whose closed-loop eigenvalues sit at
different distances from the origin, plus one controller with too much gain, and
compare how they recover the drone from a bad starting state.

Run:  python exp5_state_feedback.py   (add --show to open windows)
"""
import numpy as np
import matplotlib.pyplot as plt
from lab_common import (planar_quadcopter, drone_pole_placement, simulate, unit_circle, plot_eigs, draw_quad,
                        save, save_anim, header,
                        C_BLUE, C_ORANGE, C_RED, C_GREEN, C_PURPLE, C_GREY)

header("Experiment 5: state feedback and pole placement")
d = planar_quadcopter()
A, B, Ts = d["A"], d["B"], d["Ts"]


def design(radius):
    return drone_pole_placement(A, B, radius)


designs = {
    "slow    (|λ| ≈ 0.92)": (design(0.92), C_BLUE),
    "medium  (|λ| ≈ 0.80)": (design(0.80), C_GREEN),
    "fast    (|λ| ≈ 0.55)": (design(0.55), C_ORANGE),
}
K_med = designs["medium  (|λ| ≈ 0.80)"][0]
designs["too much gain (3 × medium)"] = (3 * K_med, C_RED)

x0 = np.array([-1.5, -0.5, np.deg2rad(10), 0, 0, 0])   # left, low, tilted 10°
N = 60
t = np.arange(N + 1) * Ts
results = {}
for name, (K, col) in designs.items():
    eig = np.linalg.eigvals(A - B @ K)
    X, U = simulate(A, B, x0, N, K=K)
    results[name] = (X, U, eig, col)
    stable = np.all(np.abs(eig) < 1)
    print(f"\n{name}\n  K =\n{K}\n  eig(A-BK) = {np.round(eig, 3)}\n  max|λ| = {np.abs(eig).max():.3f}"
          f" -> {'stable' if stable else 'UNSTABLE'};  |y| after 6 s = {abs(X[-1,0]):.2e} m;"
          f"  max |tau| = {np.abs(U[:,1]).max():.3g} N m")

# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(15, 9))
gs = fig.add_gridspec(4, 2, width_ratios=[1, 1.5], hspace=0.6, wspace=0.25)

ax = fig.add_subplot(gs[:, 0])
unit_circle(ax, "Closed-loop eigenvalues of A − BK")
plot_eigs(ax, np.linalg.eigvals(A), color=C_GREY, label="open loop (A)", annotate=False)
for name, (X, U, eig, col) in results.items():
    plot_eigs(ax, eig, color=col, label=name.split("(")[0].strip(), marker="o", size=55,
              jitter=0, annotate=False)
ax.set_xlim(-1.9, 1.4)
ax.set_ylim(-1.65, 1.65)
ax.annotate("λ = −1.66\n(outside → unstable)", (-1.655, 0), xytext=(-1.75, 0.55),
            fontsize=9, color=C_RED, arrowprops=dict(arrowstyle="->", color=C_RED))
ax.legend(fontsize=9, loc="lower left")

rows = [("y [m]", lambda X, U: X[:, 0]), ("z [m]", lambda X, U: X[:, 1]),
        ("θ [deg]", lambda X, U: np.rad2deg(X[:, 2])), ("torque τ [N m]", None)]
for r, (lab, f) in enumerate(rows):
    a = fig.add_subplot(gs[r, 1])
    for name, (X, U, eig, col) in results.items():
        bad = np.where(np.abs(X[:, 0]) > 3)[0]          # stop drawing once diverged
        cut = bad[0] + 1 if len(bad) else len(X)
        X, U = X[:cut], U[:cut - 1] if cut < len(X) else U
        if f is None:
            a.step(t[:len(U)], U[:, 1], where="post", color=col, lw=1.6)
        else:
            a.plot(t[:len(X)], f(X, U), color=col, lw=2, label=name)
    a.axhline(0, color=C_GREY, lw=0.6)
    a.set_ylabel(lab, fontsize=9)
    lim = {0: (-2.5, 1.5), 1: (-0.8, 0.6), 2: (-25, 25), 3: (-0.3, 0.3)}[r]
    a.set_ylim(*lim)
    if r == 0:
        a.legend(fontsize=8, ncol=2, loc="lower right")
        a.text(2.0, 1.0, "red: diverges (curve stopped at |y| > 3 m)", color=C_RED, fontsize=9)
    if r == 3:
        a.set_xlabel("time [s]")
fig.suptitle("Experiment 5 — feedback moves the eigenvalues: closer to 0 = faster, "
             "outside the circle = unstable", fontsize=13, weight="bold")
save(fig, "exp5_state_feedback.png")

# ----------------------------------------------------------------------------
# Animation: medium design vs too-much-gain design
# ----------------------------------------------------------------------------
Xm = results["medium  (|λ| ≈ 0.80)"][0]
Xr = results["too much gain (3 × medium)"][0]
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))


def update(k):
    for a, X, col, title in ((axs[0], Xm, C_GREEN, "medium gain: max|λ| = 0.80"),
                             (axs[1], Xr, C_RED, "3× gain: max|λ| > 1")):
        a.clear()
        draw_quad(a, 0, 0, 0, color=C_GREY, alpha=0.3)
        a.plot(X[:k + 1, 0], X[:k + 1, 1], ":", color=col)
        if np.all(np.abs(X[k, :3]) < 5):
            draw_quad(a, X[k, 0], X[k, 1], X[k, 2], color=col)
        else:
            a.text(-0.5, 0.2, "diverged!", color=C_RED, fontsize=16, weight="bold")
        a.set_xlim(-2.2, 1.2)
        a.set_ylim(-1.0, 0.8)
        a.set_aspect("equal")
        a.set_title(f"{title}\nt = {k*Ts:.1f} s", fontsize=10)
        a.set_xlabel("y [m]")


save_anim(fig, update, range(0, 41), "exp5_feedback_recovery.gif", fps=10)
