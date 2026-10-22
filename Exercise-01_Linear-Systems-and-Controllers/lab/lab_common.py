"""
lab_common.py — shared models and helpers for the Exercise 1 lab.

Two real systems are used throughout the lab:

1. PLANAR QUADCOPTER (a drone flying in a vertical plane), linearised about hover
   state  x = [y, z, theta, vy, vz, omega]   (6 states)
          y      horizontal position   [m]
          z      altitude              [m]   (deviation from hover height)
          theta  pitch / tilt angle    [rad]
          vy, vz velocities            [m/s]
          omega  pitch rate            [rad/s]
   input  u = [dF, tau]                (2 inputs)
          dF     total thrust minus hover thrust  [N]
          tau    pitching torque (difference between the two rotors) [N m]

   Linearised continuous-time physics (small tilt, sin(theta) ~ theta):
          y''     = -g * theta - c * y'      (tilting the thrust pushes sideways)
          z''     =  dF / m    - c * z'      (extra thrust lifts)
          theta'' =  tau / J                 (torque spins the drone)

2. ROBOT ARM JOINT (a single link hanging under gravity), linearised about hanging down
   state  x = [theta, omega]   input u = tau (motor torque)
          theta'' = -(g/l) theta - (b/(m l^2)) theta' + tau/(m l^2)

Both are converted to DISCRETE time x_{k+1} = A x_k + B u_k with a zero-order hold
(the input is held constant between samples), exactly the model class used in
the exercise.
"""
from pathlib import Path
import sys

import numpy as np
import matplotlib

if "--show" not in sys.argv:          # save to files unless the user asks for windows
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.linalg import expm
from scipy.signal import place_poles

np.set_printoptions(precision=4, suppress=True)

FIG_DIR = Path(__file__).resolve().parent / "figures"
FIG_DIR.mkdir(exist_ok=True)
SHOW = "--show" in sys.argv

# consistent colours across all figures
C_BLUE, C_ORANGE, C_GREEN, C_RED, C_PURPLE, C_GREY = (
    "#2a6fdb", "#e8833a", "#2e9e5b", "#d64545", "#7b52c7", "#888888")


# ----------------------------------------------------------------------------
# Models
# ----------------------------------------------------------------------------
def c2d(Ac, Bc, Ts):
    """Exact zero-order-hold discretisation:  A = e^{Ac Ts},  B = int_0^Ts e^{Ac s} ds Bc."""
    n, m = Bc.shape
    M = np.zeros((n + m, n + m))
    M[:n, :n], M[:n, n:] = Ac, Bc
    E = expm(M * Ts)
    return E[:n, :n], E[:n, n:]


def planar_quadcopter(Ts=0.1, m=0.5, J=0.004, drag=0.2, g=9.81):
    Ac = np.array([
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1],
        [0, 0, -g, -drag, 0, 0],
        [0, 0, 0, 0, -drag, 0],
        [0, 0, 0, 0, 0, 0]], dtype=float)
    Bc = np.array([
        [0, 0], [0, 0], [0, 0],
        [0, 0],
        [1 / m, 0],
        [0, 1 / J]], dtype=float)
    A, B = c2d(Ac, Bc, Ts)
    return dict(A=A, B=B, Ac=Ac, Bc=Bc, Ts=Ts,
                states=["y", "z", "theta", "vy", "vz", "omega"],
                inputs=["dF", "tau"])


def robot_arm(Ts=0.05, m=1.0, l=0.5, b=0.15, g=9.81):
    Ac = np.array([[0, 1], [-g / l, -b / (m * l**2)]])
    Bc = np.array([[0], [1 / (m * l**2)]])
    A, B = c2d(Ac, Bc, Ts)
    return dict(A=A, B=B, Ac=Ac, Bc=Bc, Ts=Ts, l=l, m=m, g=g,
                states=["theta", "omega"], inputs=["tau"])


def drone_pole_placement(A, B, radius):
    """Pole placement for the drone. It splits into two independent subsystems:
       altitude  (z, vz)                driven only by thrust dF
       lateral   (y, theta, vy, omega)  driven only by torque tau
    so we place 2 poles for altitude and 4 for lateral, then assemble the 2x6 gain K.
    All closed-loop eigenvalues end up between 0.85*radius and radius."""
    LAT, ALT = [0, 2, 3, 5], [1, 4]
    K = np.zeros((2, 6))
    p_lat = radius * np.array([1.0, 0.95, 0.9, 0.85])
    p_alt = radius * np.array([1.0, 0.9])
    K[1, LAT] = place_poles(A[np.ix_(LAT, LAT)], B[LAT, 1:2], p_lat).gain_matrix[0]
    K[0, ALT] = place_poles(A[np.ix_(ALT, ALT)], B[ALT, 0:1], p_alt).gain_matrix[0]
    return K


# ----------------------------------------------------------------------------
# Linear-algebra tools (the formulas from the Concepts file)
# ----------------------------------------------------------------------------
def ctrb(A, B, N=None):
    """W_c = [B, AB, ..., A^{N-1} B]   (N = n by default)."""
    N = A.shape[0] if N is None else N
    return np.hstack([np.linalg.matrix_power(A, i) @ B for i in range(N)])


def obsv(A, C, N=None):
    """W_o = [C; CA; ...; C A^{N-1}]."""
    N = A.shape[0] if N is None else N
    return np.vstack([C @ np.linalg.matrix_power(A, i) for i in range(N)])


def rank(M, tol=1e-9):
    return np.linalg.matrix_rank(M, tol=tol)


def distinct_eigs(A, decimals=6):
    """Distinct eigenvalues. Values are grouped by rounding, but the returned value is the
    group mean of the exact eigenvalues (rounding itself would shift lambda slightly and
    make lambda*I - A look non-singular in the PBH test)."""
    ev = np.linalg.eigvals(A)
    keys = np.round(ev, decimals)
    return np.array([ev[keys == k].mean() for k in np.unique(keys)])


def pbh_controllability(A, B):
    n = A.shape[0]
    return [(lam, rank(np.hstack([lam * np.eye(n) - A, B]))) for lam in distinct_eigs(A)]


def pbh_observability(A, C):
    n = A.shape[0]
    return [(lam, rank(np.vstack([lam * np.eye(n) - A, C]))) for lam in distinct_eigs(A)]


def simulate(A, B, x0, N, u=None, K=None, x_ref=None):
    """Simulate x_{k+1} = A x_k + B u_k for N steps.
    u      : array (N, m) of open-loop inputs, or
    K      : feedback gain, u_k = -K (x_k - x_ref)."""
    n, m = B.shape
    X = np.zeros((N + 1, n))
    U = np.zeros((N, m))
    X[0] = x0
    x_ref = np.zeros(n) if x_ref is None else x_ref
    for k in range(N):
        if K is not None:
            U[k] = -K @ (X[k] - x_ref)
        elif u is not None:
            U[k] = u[k]
        X[k + 1] = A @ X[k] + B @ U[k]
    return X, U


# ----------------------------------------------------------------------------
# Plot helpers
# ----------------------------------------------------------------------------
def unit_circle(ax, title="Eigenvalues in the complex plane"):
    t = np.linspace(0, 2 * np.pi, 400)
    ax.fill(np.cos(t), np.sin(t), color=C_GREEN, alpha=0.07)
    ax.plot(np.cos(t), np.sin(t), color=C_GREEN, lw=1.5)
    ax.axhline(0, color=C_GREY, lw=0.6)
    ax.axvline(0, color=C_GREY, lw=0.6)
    ax.set_aspect("equal")
    ax.set_xlabel("Real part")
    ax.set_ylabel("Imaginary part")
    ax.set_title(title)
    ax.text(0, -0.55, "stable region\n|λ| < 1", ha="center", va="center",
            color=C_GREEN, fontsize=9)


def plot_eigs(ax, eigs, color=C_BLUE, label=None, marker="x", size=90, jitter=0.035,
              annotate=True):
    """Scatter eigenvalues; repeated ones are fanned out slightly and annotated."""
    eigs = np.asarray(eigs)
    rounded = np.round(eigs, 4)
    vals, counts = np.unique(rounded, return_counts=True)
    first = True
    for v, c in zip(vals, counts):
        xs = np.full(c, v.real)
        ys = np.full(c, v.imag)
        if jitter and c > 1:
            ys = ys + np.linspace(-jitter, jitter, c)
        ax.scatter(xs, ys, marker=marker, s=size, color=color, lw=2,
                   label=label if first else None, zorder=5)
        first = False
        if c > 1 and annotate:
            ax.annotate(f"×{c}", (v.real, v.imag), xytext=(8, 8),
                        textcoords="offset points", color=color, fontsize=10, weight="bold")


def draw_quad(ax, y, z, theta, color=C_BLUE, alpha=1.0, arm=0.25, thrust=None):
    """Draw a simple 2-D drone: a bar with two rotors, rotated by theta (CCW)."""
    c, s = np.cos(theta), np.sin(theta)
    left = np.array([y - arm * c, z - arm * s])
    right = np.array([y + arm * c, z + arm * s])
    artists = []
    artists += ax.plot([left[0], right[0]], [left[1], right[1]], color=color, lw=4,
                       alpha=alpha, solid_capstyle="round")
    up = np.array([-s, c])  # body "up" direction
    for p in (left, right):
        mast = p + 0.06 * up
        artists += ax.plot([p[0], mast[0]], [p[1], mast[1]], color=color, lw=2, alpha=alpha)
        rot_l, rot_r = mast - 0.09 * np.array([c, s]), mast + 0.09 * np.array([c, s])
        artists += ax.plot([rot_l[0], rot_r[0]], [rot_l[1], rot_r[1]], color=C_GREY,
                           lw=2.5, alpha=alpha)
    artists.append(ax.scatter([y], [z], s=40, color=color, alpha=alpha, zorder=6))
    if thrust is not None:  # thrust arrow along body up
        tip = np.array([y, z]) + thrust * up
        artists.append(ax.annotate("", xy=tip, xytext=(y, z),
                                   arrowprops=dict(arrowstyle="->", color=C_ORANGE, lw=2)))
    return artists


def save(fig, name):
    path = FIG_DIR / name
    fig.savefig(path, dpi=110, bbox_inches="tight")
    print(f"  saved {path.relative_to(FIG_DIR.parent)}")
    if SHOW:
        plt.show()
    plt.close(fig)


def save_anim(fig, update, frames, name, fps=15):
    anim = FuncAnimation(fig, update, frames=frames, interval=1000 / fps, blit=False)
    path = FIG_DIR / name
    anim.save(path, writer=PillowWriter(fps=fps), dpi=75)
    print(f"  saved {path.relative_to(FIG_DIR.parent)}")
    if SHOW:
        plt.show()
    plt.close(fig)


def header(title):
    print("\n" + "=" * 72 + f"\n{title}\n" + "=" * 72)
