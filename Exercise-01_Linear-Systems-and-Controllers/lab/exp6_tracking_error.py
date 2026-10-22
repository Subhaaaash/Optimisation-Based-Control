"""
Experiment 6 — Task 1 live: the tracking error of u = -K (x - x_ref).

Task 1 result:   x_inf - x_ref = (I - A + BK)^{-1} (A - I) x_ref
                 zero  <=>  A x_ref = x_ref   (x_ref is an equilibrium WITHOUT input)

Two systems, two outcomes:
  * Drone, target "hover 2 m right, 1 m up":  positions are integrators, so a hovering
    drone stays where it is with u = 0  ->  A x_ref = x_ref  ->  ZERO error.
  * Robot arm, target "hold the link at 30°": gravity pulls the link back, so holding it
    needs a constant torque, but the P-controller gives zero torque at the target
    ->  A x_ref != x_ref  ->  the arm SAGS below the target, for every gain.

Run:  python exp6_tracking_error.py   (add --show to open windows)
"""
import numpy as np
import matplotlib.pyplot as plt
from lab_common import (planar_quadcopter, robot_arm, drone_pole_placement, simulate, draw_quad, save, save_anim,
                        header, C_BLUE, C_ORANGE, C_RED, C_GREEN, C_PURPLE, C_GREY)

header("Experiment 6: tracking error (Task 1)")

# ----------------------------------------------------------------------------
# Part A: drone -> zero error
# ----------------------------------------------------------------------------
dq = planar_quadcopter()
Aq, Bq, Ts_q = dq["A"], dq["B"], dq["Ts"]
Kq = drone_pole_placement(Aq, Bq, 0.80)          # the "medium" design from Experiment 5
xq_ref = np.array([2.0, 1.0, 0, 0, 0, 0])
Xq, _ = simulate(Aq, Bq, np.zeros(6), 80, K=Kq, x_ref=xq_ref)
I6 = np.eye(6)
err_q = np.linalg.solve(I6 - Aq + Bq @ Kq, (Aq - I6) @ xq_ref)
print(f"\nDRONE: A x_ref == x_ref ? {np.allclose(Aq @ xq_ref, xq_ref)}")
print(f"  predicted error (I-A+BK)^-1 (A-I) x_ref = {err_q}")
print(f"  simulated final state                   = {Xq[-1]}")

# ----------------------------------------------------------------------------
# Part B: robot arm -> steady-state error
# ----------------------------------------------------------------------------
arm = robot_arm()
A, B, Ts, m, l, g = arm["A"], arm["B"], arm["Ts"], arm["m"], arm["l"], arm["g"]
b_damp = 0.15
theta_ref = np.deg2rad(30)
x_ref = np.array([theta_ref, 0.0])
I2 = np.eye(2)
print(f"\nARM: A x_ref = {A @ x_ref}  vs  x_ref = {x_ref}  -> equal? "
      f"{np.allclose(A @ x_ref, x_ref)}  (gravity pulls it back)")


def arm_gain(kp, zeta=0.7):
    """K = [kp, kd] with kd chosen for damping ratio zeta (kd does not affect the error)."""
    J = m * l**2
    wn = np.sqrt(g / l + kp / J)
    kd = 2 * zeta * wn * J - b_damp
    return np.array([[kp, kd]])


def predicted(K):
    x_inf = np.linalg.solve(I2 - (A - B @ K), B @ K @ x_ref)
    err = np.linalg.solve(I2 - A + B @ K, (A - I2) @ x_ref)
    return x_inf, err


kps = [2, 5, 15, 50]
cols = [C_RED, C_ORANGE, C_BLUE, C_GREEN]
N = 60
t = np.arange(N + 1) * Ts
arm_runs = {}
for kp in kps:
    K = arm_gain(kp)
    X, U = simulate(A, B, np.zeros(2), N, K=K, x_ref=x_ref)
    x_inf, err = predicted(K)
    eig = np.linalg.eigvals(A - B @ K)
    arm_runs[kp] = (X, U, x_inf)
    print(f"  kp = {kp:3d}: max|eig(A-BK)| = {np.abs(eig).max():.3f}  "
          f"theta_inf: formula {np.rad2deg(x_inf[0]):6.2f}°, sim {np.rad2deg(X[-1,0]):6.2f}°,"
          f"  error {np.rad2deg(err[0]):+6.2f}°,  "
          f"kp/(kp+mgl) = {kp/(kp+m*g*l):.3f}")

kp_sweep = np.logspace(-0.5, 3, 80)
err_sweep, stable_sweep = [], []
for kp in kp_sweep:
    K = arm_gain(kp)
    err_sweep.append(np.rad2deg(predicted(K)[1][0]))
    stable_sweep.append(np.all(np.abs(np.linalg.eigvals(A - B @ K)) < 1))
err_sweep, stable_sweep = np.array(err_sweep), np.array(stable_sweep)
unstable_from = kp_sweep[~stable_sweep].min() if (~stable_sweep).any() else None
print(f"  gain sweep: closed loop becomes unstable from kp ≈ "
      f"{unstable_from:.0f}" if unstable_from else "  gain sweep: stable throughout")


def draw_arm(ax, theta, color, alpha=1.0, lw=6, label=None):
    tip = np.array([l * np.sin(theta), -l * np.cos(theta)])
    ax.plot([0, tip[0]], [0, tip[1]], color=color, lw=lw, alpha=alpha,
            solid_capstyle="round", label=label)
    ax.scatter([tip[0]], [tip[1]], s=80, color=color, alpha=alpha, zorder=5)


# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(15, 10))
gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.3)

ax = fig.add_subplot(gs[0, 0])
ax.scatter([0], [0], s=150, color="k", zorder=6)
draw_arm(ax, theta_ref, C_GREY, alpha=0.35, lw=10, label="target 30°")
for kp, col in zip(kps, cols):
    draw_arm(ax, arm_runs[kp][0][-1, 0], col, lw=3,
             label=f"kp = {kp}: settles at {np.rad2deg(arm_runs[kp][0][-1,0]):.1f}°")
ax.plot([0, 0], [0, -0.6], ":", color=C_GREY)
ax.annotate("gravity ↓", (-0.13, -0.62), fontsize=9, color=C_GREY)
ax.set_aspect("equal")
ax.set_xlim(-0.75, 0.4)
ax.set_ylim(-0.68, 0.1)
ax.set_title("Robot arm: where it actually ends up")
ax.legend(fontsize=8, loc="lower left", bbox_to_anchor=(0.0, 0.12))

ax = fig.add_subplot(gs[0, 1:])
for kp, col in zip(kps, cols):
    X, U, x_inf = arm_runs[kp]
    ax.plot(t, np.rad2deg(X[:, 0]), color=col, lw=2, label=f"kp = {kp}")
    ax.axhline(np.rad2deg(x_inf[0]), color=col, ls=":", lw=1.2)
ax.axhline(30, color="k", ls="--", lw=1.5, label="θ_ref = 30°")
ax.set_xlabel("time [s]")
ax.set_ylabel("θ [deg]")
ax.set_title("Arm angle: dotted lines = Task 1 prediction  x∞ = (I − A + BK)⁻¹ BK x_ref")
ax.legend(fontsize=9, ncol=3)

ax = fig.add_subplot(gs[1, 0])
ax.semilogx(kp_sweep[stable_sweep], err_sweep[stable_sweep], color=C_BLUE, lw=2,
            label="formula (stable K)")
if unstable_from:
    ax.axvspan(unstable_from, kp_sweep.max(), color=C_RED, alpha=0.12)
    ax.text(unstable_from * 1.1, -25, "A − BK\nunstable", color=C_RED, fontsize=9)
for kp, col in zip(kps, cols):
    X = arm_runs[kp][0]
    ax.scatter(kp, np.rad2deg(X[-1, 0]) - 30, color=col, s=70, zorder=5)
ax.axhline(0, color=C_GREY, lw=0.6)
ax.set_xlabel("proportional gain kp  (log scale)")
ax.set_ylabel("steady-state error θ∞ − θ_ref [deg]")
ax.set_title("More gain → smaller error, never zero")
ax.legend(fontsize=9)

ax = fig.add_subplot(gs[1, 1])
needed = m * g * l * theta_ref            # linearised holding torque m g l theta
labels = ["torque needed\nto hold 30°", "P-controller torque\nAT the target", ]
ax.bar(labels, [needed, 0], color=[C_GREY, C_RED])
ax.text(1, 0.15, "u = −K(x_ref − x_ref) = 0", ha="center", color=C_RED, fontsize=9)
ax.set_ylabel("torque [N m]")
ax.set_title("Why: at the target the controller pushes with zero torque")

ax = fig.add_subplot(gs[1, 2])
tq = np.arange(81) * Ts_q
ax.plot(tq, Xq[:, 0], color=C_BLUE, lw=2, label="drone y")
ax.plot(tq, Xq[:, 1], color=C_GREEN, lw=2, label="drone z")
ax.axhline(2, color=C_BLUE, ls="--", lw=1)
ax.axhline(1, color=C_GREEN, ls="--", lw=1)
ax.set_xlabel("time [s]")
ax.set_ylabel("position [m]")
ax.set_title("Drone: A x_ref = x_ref  →  zero error ✓")
ax.legend(fontsize=9, loc="lower right")

fig.suptitle("Experiment 6 — Task 1: zero tracking error only if x_ref is an equilibrium "
             "without input (A x_ref = x_ref)", fontsize=13, weight="bold")
save(fig, "exp6_tracking_error.png")

# ----------------------------------------------------------------------------
# Animation: arm with kp = 5 and kp = 50 trying to reach 30°
# ----------------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(9, 4.4))


def update(k):
    for a, kp, col in ((axs[0], 5, C_ORANGE), (axs[1], 50, C_GREEN)):
        a.clear()
        X = arm_runs[kp][0]
        a.scatter([0], [0], s=120, color="k", zorder=6)
        draw_arm(a, theta_ref, C_GREY, alpha=0.3, lw=10)
        draw_arm(a, X[k, 0], col, lw=5)
        e = np.rad2deg(X[k, 0] - theta_ref)
        a.set_aspect("equal")
        a.set_xlim(-0.35, 0.45)
        a.set_ylim(-0.6, 0.15)
        a.set_title(f"kp = {kp}   t = {k*Ts:.2f} s\nθ = {np.rad2deg(X[k,0]):5.1f}°   "
                    f"error = {e:+5.1f}°", fontsize=10)
        a.text(0.25, -0.12, "grey = 30° target", fontsize=8, color=C_GREY)


save_anim(fig, update, list(range(0, N + 1)) + [N] * 10, "exp6_arm_sag.gif", fps=15)
