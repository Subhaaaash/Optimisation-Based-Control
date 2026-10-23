"""
Optimization-Based Control Methods - Lecture 2 examples, verified numerically.

1. Thermal pre-heating (N = 3, a = 0.5, r = 16, x0 = 0, T = 10):
   backward sweep p_k, optimal inputs u_k, trajectory x_k, and the
   "inputs double every step" property u_{k+1} = u_k / (1 - a).
   Cross-check: brute-force minimisation of the open-loop cost.
2. Deterministic scheduling of A, B, C, D (A before B, C before D).
3. Two-game chess match (timid vs bold), probability of winning the match.

Run:  python lecture_examples.py
"""
import itertools

import numpy as np
from scipy.optimize import minimize


# --------------------------------------------------------------------------
# 1. Thermal pre-heating
# --------------------------------------------------------------------------
def heating(N=3, a=0.5, r=16.0, T=10.0, x0=0.0):
    print("=== 1. Thermal pre-heating ===")
    bt = [a * (1 - a) ** (N - 1 - k) for k in range(N)]      # b~_k
    p = [0.0] * (N + 1)
    p[N] = r
    for k in range(N - 1, -1, -1):
        p[k] = p[k + 1] / (1 + bt[k] ** 2 * p[k + 1])
    print("b~_k =", bt)
    print("p_k  =", p, "(expected 2.56, 8/3, 3.2, 16)")

    x, us, xs = x0, [], [x0]
    for k in range(N):
        u = bt[k] * p[k] * (T - (1 - a) ** (N - k) * x)
        x = (1 - a) * x + a * u
        us.append(u)
        xs.append(x)
    print("u*_k =", np.round(us, 6), "(expected 3.2, 6.4, 12.8)")
    print("x_k  =", np.round(xs, 6), "(expected 0, 1.6, 4.0, 8.4)")
    print("ratio u_{k+1}/u_k =", [us[k + 1] / us[k] for k in range(N - 1)],
          "= 1/(1-a) =", 1 / (1 - a))

    def J(u):
        x = x0
        for uk in u:
            x = (1 - a) * x + a * uk
        return r * (x - T) ** 2 + np.sum(np.square(u))

    J_dp = p[0] * ((1 - a) ** N * x0 - T) ** 2
    res = minimize(J, np.zeros(N))
    print(f"J0(x0) from DP = {J_dp:.6f}, cost of DP inputs = {J(np.array(us)):.6f}, "
          f"brute force = {res.fun:.6f}")
    assert abs(J_dp - res.fun) < 1e-5 and abs(J_dp - J(np.array(us))) < 1e-9
    print()


# --------------------------------------------------------------------------
# 2. Deterministic scheduling (numbers from the lecture slide)
# --------------------------------------------------------------------------
def scheduling():
    print("=== 2. Deterministic scheduling ===")
    start = {"A": 5, "C": 3}                                  # S_A, S_C
    # setup transition costs C_mn (from operation m to operation n), read off the slide
    C = {"AB": 3, "AC": 2, "AD": 6, "BC": 4, "BD": 2,
         "CA": 1, "CB": 5, "CD": 4, "DA": 2, "DB": 3}

    def cost(seq):
        return start[seq[0]] + sum(C[seq[i] + seq[i + 1]] for i in range(3))

    feasible = [ "".join(s) for s in itertools.permutations("ABCD")
                 if s.index("A") < s.index("B") and s.index("C") < s.index("D")]
    costs = {s: cost(s) for s in feasible}
    for s, c in sorted(costs.items(), key=lambda t: t[1]):
        print(f"  {'->'.join(s)} : {c}")
    best = min(costs, key=costs.get)
    print("optimal schedule:", "->".join(best), "with cost", costs[best], "(expected C->A->B->D, 9)\n")
    assert best == "CABD" and costs[best] == 9


# --------------------------------------------------------------------------
# 3. Chess match (tie after two games -> sudden death played bold, win prob p_w)
# --------------------------------------------------------------------------
def chess(pw=0.45, pd=0.9):
    print(f"=== 3. Chess match, p_w = {pw}, p_d = {pd} ===")
    J2 = lambda x: 1.0 if x > 0 else (pw if x == 0 else 0.0)

    def q(x, J_next, action):
        if action == "timid":
            return pd * J_next(x) + (1 - pd) * J_next(x - 1)
        return pw * J_next(x + 1) + (1 - pw) * J_next(x - 1)

    J1, mu1 = {}, {}
    for x in (-1, 0, 1):
        vals = {a_: q(x, J2, a_) for a_ in ("timid", "bold")}
        mu1[x] = max(vals, key=vals.get)
        J1[x] = vals[mu1[x]]
        print(f"  k=1, x1={x:+d}: timid {vals['timid']:.4f}, bold {vals['bold']:.4f} -> {mu1[x]}")
    J1f = lambda x: J1[x]
    vals0 = {a_: q(0, J1f, a_) for a_ in ("timid", "bold")}
    mu0 = max(vals0, key=vals0.get)
    print(f"  k=0, x0= 0: timid {vals0['timid']:.4f}, bold {vals0['bold']:.4f} -> {mu0}")
    closed_form = pw * pd + pw ** 2 * (2 - pd - pw)
    print(f"P(win match) = {vals0[mu0]:.4f}  (closed form {closed_form:.4f}) > 1/2 although p_w < 1/2\n")
    assert mu1 == {-1: "bold", 0: "bold", 1: "timid"} and mu0 == "bold"


if __name__ == "__main__":
    heating()
    scheduling()
    chess()
