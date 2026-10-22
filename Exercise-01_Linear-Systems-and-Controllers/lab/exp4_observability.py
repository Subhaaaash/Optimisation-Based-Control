"""
Experiment 4 — Observability: can the sensors tell us the full state?

Sensor options for the drone (C matrices):
    altimeter only         y_k = z
    GPS (horizontal) only  y_k = y
    gyroscope only         y_k = omega
    GPS + altimeter        y_k = [y, z]

1. rank of W_o and the PBH test for each sensor set.
2. "Two drones, same sensor readings": with the altimeter only, two initial states that
   differ by a vector in ker(W_o) give IDENTICAL measurements, yet fly completely
   differently.
3. Reconstructing x_0 from the first n measurements:  x0_hat = pinv(W_o) [y_0; ...; y_{n-1}].
   Exact when observable; otherwise only the observable component is recovered.

Run:  python exp4_observability.py   (add --show to open windows)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import null_space
from lab_common import (planar_quadcopter, obsv, rank, pbh_observability, simulate,
                        draw_quad, save, save_anim, header,
                        C_BLUE, C_ORANGE, C_RED, C_GREEN, C_PURPLE, C_GREY)

header("Experiment 4: observability")
d = planar_quadcopter()
A, B, Ts = d["A"], d["B"], d["Ts"]
n = A.shape[0]
names = d["states"]


def sel(*idx):
    C = np.zeros((len(idx), n))
    for r, i in enumerate(idx):
        C[r, i] = 1
    return C


sensors = {
    "altimeter (z)":         sel(1),
    "gyroscope (ω)":         sel(5),
    "GPS horizontal (y)":    sel(0),
    "GPS + altimeter (y,z)": sel(0, 1),
}
ranks = {}
for name, C in sensors.items():
    Wo = obsv(A, C)
    ranks[name] = rank(Wo)
    ker = null_space(Wo, rcond=1e-9)
    hidden = [names[i] for i in range(n) if ker.shape[1] and np.linalg.norm(ker[i]) > 1e-6]
    print(f"\n{name:24s} rank W_o = {ranks[name]} of {n}"
          f"   unobservable states: {hidden if hidden else 'none'}")
    for lam, r in pbh_observability(A, C):
        print(f"    PBH lambda = {lam.real:.4f}: rank = {r} -> "
              f"{'observable' if r == n else 'UNOBSERVABLE'} mode")

# ----------------------------------------------------------------------------
# Two drones, identical altimeter readings
# ----------------------------------------------------------------------------
C_alt = sensors["altimeter (z)"]
Wo_alt = obsv(A, C_alt)
x0_a = np.array([0, 0.3, 0, 0, 0.1, 0])                   # drone A: just a bit high
hidden_dir = np.array([1.0, 0, np.deg2rad(4), 0.4, 0, 0])  # lies in ker(W_o): no z part
print("\n  W_o (altimeter) @ hidden direction =", Wo_alt @ hidden_dir, " (zero -> invisible)")
x0_b = x0_a + hidden_dir                                    # drone B: same + hidden part
N = 40
t = np.arange(N + 1) * Ts
Xa, _ = simulate(A, B, x0_a, N)
Xb, _ = simulate(A, B, x0_b, N)
Ya, Yb = Xa @ C_alt.T, Xb @ C_alt.T
print(f"  max |altimeter_A - altimeter_B| over 4 s = {np.abs(Ya - Yb).max():.2e}")
print(f"  but after 4 s: y_A = {Xa[-1,0]:+.2f} m, y_B = {Xb[-1,0]:+.2f} m")

# ----------------------------------------------------------------------------
# Reconstruct x0 from n measurements
# ----------------------------------------------------------------------------
x0_true = np.array([0.5, -0.2, np.deg2rad(5), 0.3, 0.1, np.deg2rad(-3)])
est = {}
for name in ["GPS + altimeter (y,z)", "altimeter (z)"]:
    C = sensors[name]
    Wo = obsv(A, C)
    X, _ = simulate(A, B, x0_true, n - 1)
    Yvec = (X @ C.T).reshape(-1)                     # [y_0; y_1; ...; y_{n-1}]
    est[name] = np.linalg.pinv(Wo) @ Yvec
    print(f"\n  {name}: true x0 = {x0_true}\n{'':4s}estimated x0 = {est[name]}")

# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(15, 9.5))
gs = fig.add_gridspec(2, 3, hspace=0.4, wspace=0.32)

ax = fig.add_subplot(gs[0, 0])
cols = [C_RED, C_RED, C_ORANGE, C_GREEN]
bars = ax.barh(list(ranks.keys()), list(ranks.values()), color=cols)
ax.axvline(n, color=C_GREY, ls="--")
for b, v in zip(bars, ranks.values()):
    ax.text(v + 0.08, b.get_y() + b.get_height() / 2, f"{v}/6", va="center")
ax.set_xlim(0, 7)
ax.set_xlabel("rank W_o   (6 = observable)")
ax.set_title("Which sensors make the drone observable?")
ax.invert_yaxis()

ax = fig.add_subplot(gs[0, 1])
ax.plot(t, Ya[:, 0], color=C_BLUE, lw=4, alpha=0.5, label="drone A: altimeter reading")
ax.plot(t, Yb[:, 0], "--", color=C_RED, lw=2, label="drone B: altimeter reading")
ax.set_xlabel("time [s]")
ax.set_ylabel("measured z [m]")
ax.set_title("What the altimeter sees: identical")
ax.legend(fontsize=9)

ax = fig.add_subplot(gs[0, 2])
ax.plot(Xa[:, 0], Xa[:, 1], color=C_BLUE, lw=2, label="drone A (true path)")
ax.plot(Xb[:, 0], Xb[:, 1], color=C_RED, lw=2, label="drone B (true path)")
for k in (0, 20, 40):
    draw_quad(ax, Xa[k, 0], Xa[k, 1], Xa[k, 2], color=C_BLUE, arm=0.2, alpha=0.6)
    draw_quad(ax, Xb[k, 0], Xb[k, 1], Xb[k, 2], color=C_RED, arm=0.2, alpha=0.6)
ax.set_xlabel("y [m]")
ax.set_ylabel("z [m]")
ax.set_title("What really happens: totally different")
ax.legend(fontsize=9, loc="lower left")

ax = fig.add_subplot(gs[1, :2])
w = 0.27
xs = np.arange(n)
ax.bar(xs - w, x0_true, w, color=C_GREY, label="true x₀")
ax.bar(xs, est["GPS + altimeter (y,z)"], w, color=C_GREEN, label="estimate with GPS + altimeter")
ax.bar(xs + w, est["altimeter (z)"], w, color=C_RED, label="estimate with altimeter only")
ax.set_xticks(xs, ["y", "z", "θ", "v_y", "v_z", "ω"])
ax.axhline(0, color=C_GREY, lw=0.6)
ax.set_title("Reconstructing x₀ = pinv(W_o)·[y₀ … y₅]: exact when observable, "
             "only the visible part otherwise")
ax.legend(fontsize=9)

ax = fig.add_subplot(gs[1, 2])
ax.axis("off")
ax.text(0, 1, "Why?\n\n"
        "Altimeter only:  W_o = [C; CA; …]\n"
        "only ever involves z and v_z,\n"
        "so ker(W_o) contains every direction\n"
        "made of y, θ, v_y, ω.\n\n"
        "Drone B = drone A + a vector in ker(W_o)\n"
        "→ W_o x₀ᴮ = W_o x₀ᴬ\n"
        "→ identical measurements forever.\n\n"
        "Gyro only: rank 1 — the gyro sees the\n"
        "rotation rate, nothing else.\n\n"
        "GPS (y) alone is close (rank 4):\n"
        "y reveals v_y, θ and ω through the\n"
        "physics, but not the altitude.",
        va="top", fontsize=10, family="monospace")

fig.suptitle("Experiment 4 — observability: can we work out the state from the sensors?",
             fontsize=13, weight="bold")
save(fig, "exp4_observability.png")

# ----------------------------------------------------------------------------
# Animation: the sensor's view vs reality
# ----------------------------------------------------------------------------
fig, (axl, axr) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw=dict(width_ratios=[1, 2]))


def update(k):
    for a in (axl, axr):
        a.clear()
    axl.plot(t[:k + 1], Ya[:k + 1, 0], color=C_BLUE, lw=4, alpha=0.5)
    axl.plot(t[:k + 1], Yb[:k + 1, 0], "--", color=C_RED, lw=2)
    axl.set_xlim(0, t[-1])
    axl.set_ylim(0, 0.8)
    axl.set_title("altimeter readings (A and B)", fontsize=10)
    axl.set_xlabel("time [s]")
    axr.plot(Xa[:k + 1, 0], Xa[:k + 1, 1], ":", color=C_BLUE)
    axr.plot(Xb[:k + 1, 0], Xb[:k + 1, 1], ":", color=C_RED)
    draw_quad(axr, Xa[k, 0], Xa[k, 1], Xa[k, 2], color=C_BLUE, arm=0.2)
    draw_quad(axr, Xb[k, 0], Xb[k, 1], Xb[k, 2], color=C_RED, arm=0.2)
    axr.set_xlim(-2.8, 1.0)
    axr.set_ylim(-0.3, 1.2)
    axr.set_aspect("equal")
    axr.set_title(f"reality   t = {k*Ts:.1f} s", fontsize=10)
    axr.set_xlabel("y [m]")


save_anim(fig, update, range(0, N + 1), "exp4_same_sensor_different_reality.gif", fps=10)
