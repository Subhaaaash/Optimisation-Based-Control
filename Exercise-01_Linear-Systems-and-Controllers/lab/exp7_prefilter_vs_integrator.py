"""
Experiment 7 — Tasks 2 and 3 live: removing the robot arm's sag.

Same arm and target (30°) as Experiment 6. Three controllers:
  P      u = -K (x - x_ref)                         (Task 1: sags)
  P + F  u = -K (x - F x_ref),  F = diag(f, 1)      (Task 2: aims higher on purpose)
  PI     u = -Kp e - KI e^i,    e^i_{k+1} = e^i_k + Ts e_k   (Task 3)
and two situations:
  nominal  : the real arm is exactly the model
  payload  : the arm carries 30 % extra mass that the controller does not know about

Run:  python exp7_prefilter_vs_integrator.py   (add --show to open windows)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import place_poles
from lab_common import (robot_arm, save, save_anim, header,
                        C_BLUE, C_ORANGE, C_RED, C_GREEN, C_GREY)

header("Experiment 7: prefilter vs integrator on the robot arm")
nom = robot_arm()                      # model used for the design
heavy = robot_arm(m=1.3)               # real arm with payload
A, B, Ts, m, l, g = nom["A"], nom["B"], nom["Ts"], nom["m"], nom["l"], nom["g"]
b_damp = 0.15
theta_ref = np.deg2rad(30)
x_ref = np.array([theta_ref, 0.0])

# ---- P gain (as in Experiment 6, kp = 15) -------------------------------------
kp, zeta = 15.0, 0.7
J = m * l**2
wn = np.sqrt(g / l + kp / J)
K = np.array([[kp, 2 * zeta * wn * J - b_damp]])

# ---- Task 2 prefilter ---------------------------------------------------------
# Holding torque: x_ref = A x_ref + B u_bar  ->  u_bar = m g l theta_ref (exact for ZOH).
# Need K F x_ref = u_bar + K x_ref  ->  kp f theta_ref = m g l theta_ref + kp theta_ref
f = 1 + m * g * l / kp
F = np.diag([f, 1.0])
u_bar = np.linalg.lstsq(B, (np.eye(2) - A) @ x_ref, rcond=None)[0][0]
print(f"holding torque u_bar = {u_bar:.4f} N m  (m g l theta_ref = {m*g*l*theta_ref:.4f})")
print(f"prefilter f = 1 + m g l / kp = {f:.4f}  -> the controller aims at "
      f"{np.rad2deg(f*theta_ref):.1f}° to end at 30°")

# ---- Task 3 PI (integrate the angle error) -------------------------------------
# With one motor and two states, a full-state integrator always keeps eigenvalue 1
# (Task 3.3*), so we integrate only the angle error: KI = [kI, 0].
Ao = np.block([[A, np.zeros((2, 1))], [Ts * np.array([[1.0, 0.0]]), np.eye(1)]])
Bo = np.vstack([B, [[0.0]]])
p_P = np.linalg.eigvals(A - B @ K)
G = place_poles(Ao, Bo, np.concatenate([p_P, [0.85]])).gain_matrix
Kp, KI = G[:, :2], np.array([[G[0, 2], 0.0]])
print(f"P : K = {K}, eig(A-BK) = {np.round(p_P, 3)}")
print(f"PI: Kp = {Kp}, KI = {KI}")


def run(plant, ctrl, N=100):
    Ap, Bp = plant["A"], plant["B"]
    X = np.zeros((N + 1, 2))
    ei = np.zeros(2)
    U = np.zeros(N)
    for k in range(N):
        e = X[k] - x_ref
        if ctrl == "P":
            u = -K @ e
        elif ctrl == "P+F":
            u = -K @ (X[k] - F @ x_ref)
        else:
            u = -Kp @ e - KI @ ei
        U[k] = u[0]
        X[k + 1] = Ap @ X[k] + Bp @ u
        ei = ei + Ts * e
    return X, U


CTRL = [("P", "Task 1: P", C_RED), ("P+F", "Task 2: P + prefilter", C_ORANGE),
        ("PI", "Task 3: PI", C_BLUE)]
res = {}
for plant_name, plant in (("nominal", nom), ("payload +30 %", heavy)):
    for c, lab, col in CTRL:
        X, U = run(plant, c)
        res[(plant_name, c)] = (X, U)
        print(f"  {plant_name:14s} {c:4s}: final angle {np.rad2deg(X[-1,0]):6.2f}°"
              f"  error {np.rad2deg(X[-1,0]) - 30:+6.2f}°   final torque {U[-1]:.3f} N m")

# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
N = 100
t = np.arange(N + 1) * Ts
fig, axs = plt.subplots(2, 3, figsize=(15, 8.5), gridspec_kw=dict(width_ratios=[1.6, 1, 1]))
for r, plant_name in enumerate(("nominal", "payload +30 %")):
    ax = axs[r, 0]
    for c, lab, col in CTRL:
        X, U = res[(plant_name, c)]
        ax.plot(t, np.rad2deg(X[:, 0]), color=col, lw=2, label=lab)
    ax.axhline(30, color="k", ls="--", lw=1.2, label="target 30°")
    ax.set_ylabel("θ [deg]")
    ax.set_title(f"{plant_name}: arm angle")
    ax.set_ylim(0, 40)
    ax.grid(alpha=0.3)
    if r == 0:
        ax.legend(fontsize=9, loc="lower right")
    else:
        ax.set_xlabel("time [s]")

    ax = axs[r, 1]
    for c, lab, col in CTRL:
        X, U = res[(plant_name, c)]
        ax.plot(t[:-1], U, color=col, lw=2)
    mass = 1.0 if r == 0 else 1.3
    hold = mass * g * l * theta_ref
    ax.axhline(hold, color="k", ls=":", lw=1.2)
    ax.text(t[-1], hold + 0.1, "torque needed to hold 30°", ha="right", fontsize=8)
    ax.set_title("motor torque [N m]")
    ax.set_ylim(-4, 14)
    ax.grid(alpha=0.3)
    if r == 1:
        ax.set_xlabel("time [s]")

    ax = axs[r, 2]
    errs = [np.rad2deg(res[(plant_name, c)][0][-1, 0]) - 30 for c, _, _ in CTRL]
    bars = ax.bar([lab.split(":")[0] for _, lab, _ in CTRL], errs, color=[col for *_, col in CTRL])
    for bar_, e in zip(bars, errs):
        ax.text(bar_.get_x() + bar_.get_width() / 2, e + (0.4 if e >= 0 else -0.9),
                f"{e:+.2f}°", ha="center", fontsize=10)
    ax.axhline(0, color=C_GREY, lw=0.8)
    ax.set_ylim(-10, 3)
    ax.set_title("steady-state error")
fig.suptitle("Experiment 7 — the prefilter removes the sag only if the model is exact; "
             "the integrator removes it anyway", fontsize=13, weight="bold")
fig.tight_layout()
save(fig, "exp7_prefilter_vs_integrator.png")

# ----------------------------------------------------------------------------
# Animation: payload case, three arms side by side
# ----------------------------------------------------------------------------


def draw_arm(ax, theta, color, alpha=1.0, lw=6):
    tip = np.array([l * np.sin(theta), -l * np.cos(theta)])
    ax.plot([0, tip[0]], [0, tip[1]], color=color, lw=lw, alpha=alpha, solid_capstyle="round")
    ax.scatter([tip[0]], [tip[1]], s=80, color=color, alpha=alpha, zorder=5)


fig, axs = plt.subplots(1, 3, figsize=(12, 4.3))


def update(k):
    for ax, (c, lab, col) in zip(axs, CTRL):
        ax.clear()
        X = res[("payload +30 %", c)][0]
        ax.scatter([0], [0], s=120, color="k", zorder=6)
        draw_arm(ax, theta_ref, C_GREY, alpha=0.3, lw=10)
        draw_arm(ax, X[k, 0], col, lw=5)
        ax.set_aspect("equal")
        ax.set_xlim(-0.3, 0.45)
        ax.set_ylim(-0.6, 0.12)
        ax.set_title(f"{lab}\nθ = {np.rad2deg(X[k,0]):5.1f}°  error {np.rad2deg(X[k,0])-30:+5.1f}°",
                     fontsize=10)
    fig.suptitle(f"Arm with unknown +30 % payload,  t = {k*Ts:.2f} s   (grey = 30° target)",
                 fontsize=11)


save_anim(fig, update, list(range(0, 101, 2)) + [100] * 8, "exp7_payload.gif", fps=12)
