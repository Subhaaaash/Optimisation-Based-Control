"""
Optimization-Based Control Methods - Exercise 2, Task 2
Scalar linear-quadratic DP:  x_{k+1} = a x_k + b u_k (+ w_k),  N = 4
cost  q_N x_4^2 + sum_{k=0}^{3} (x_k^2 + u_k^2)

Checks
(a) symbolic: one backward DP step maps p_{k+1} x^2 to p_k x^2 with
        p_k = 1 + a^2 p_{k+1} / (1 + b^2 p_{k+1}),
        mu_k(x) = - a b p_{k+1} / (1 + b^2 p_{k+1}) * x
(b) q_N = 0  ->  mu_3 = 0
(c) stochastic: J_k = p_k x^2 + c_k with c_k = c_{k+1} + p_{k+1} sigma^2,
    gains unchanged; verified symbolically and by Monte Carlo.
Numerically: DP cost equals the brute-force open-loop minimum.

Run:  python task2_scalar_lq_dp.py
"""
import numpy as np
import sympy as sp
from scipy.optimize import minimize

# ---------------------------------------------------------------- symbolic
a, b, x, u, p, c, s2 = sp.symbols("a b x u p c sigma2", real=True)
w = sp.symbols("w", real=True)

print("=== (a) one backward step, deterministic ===")
Q = x**2 + u**2 + p * (a * x + b * u) ** 2          # x^2 + u^2 + J_{k+1}(ax+bu)
u_star = sp.solve(sp.diff(Q, u), u)[0]
print("argmin u      :", sp.simplify(u_star))
J = sp.simplify(Q.subs(u, u_star))
p_k = sp.simplify(J / x**2)
print("J_k(x) / x^2  :", p_k)
assert sp.simplify(p_k - (1 + a**2 * p / (1 + b**2 * p))) == 0
assert sp.simplify(sp.diff(Q, u, 2) - 2 * (1 + b**2 * p)) == 0
print("second derivative 2(1 + b^2 p) > 0 for p >= 0 -> minimum  OK")

print("\n=== (b) q_N = 0 ===")
print("mu_3(x) with p_4 = 0 :", u_star.subs(p, 0))
print("p_3 with p_4 = 0     :", p_k.subs(p, 0))

print("\n=== (c) one backward step with noise E[w]=0, E[w^2]=sigma^2 ===")
Qs = x**2 + u**2 + sp.expand(p * (a * x + b * u + w) ** 2 + c)
# take the expectation term by term: E[w] = 0, E[w^2] = sigma^2
Qs_E = Qs.subs(w**2, s2).subs(w, 0)
u_star_s = sp.solve(sp.diff(Qs_E, u), u)[0]
print("argmin u      :", sp.simplify(u_star_s), "(same as deterministic:",
      sp.simplify(u_star_s - u_star) == 0, ")")
Js = sp.expand(sp.simplify(Qs_E.subs(u, u_star_s)))
const = sp.simplify(Js.subs(x, 0))
print("constant term :", const, " -> c_k = c_{k+1} + p_{k+1} sigma^2")
assert sp.simplify(const - (c + p * s2)) == 0

# ---------------------------------------------------------------- numeric
def backward(a_, b_, qN, N=4):
    P = [0.0] * (N + 1)
    P[N] = qN
    L = [0.0] * N
    for k in range(N - 1, -1, -1):
        L[k] = a_ * b_ * P[k + 1] / (1 + b_**2 * P[k + 1])     # mu_k(x) = -L_k x
        P[k] = 1 + a_**2 * P[k + 1] / (1 + b_**2 * P[k + 1])
    return P, L


def cost(useq, a_, b_, qN, x0):
    xk, J_ = x0, 0.0
    for uk in useq:
        J_ += xk**2 + uk**2
        xk = a_ * xk + b_ * uk
    return J_ + qN * xk**2


print("\n=== DP vs brute force (random a, b, q_N, x0) ===")
rng = np.random.default_rng(0)
for _ in range(5):
    a_, b_ = rng.uniform(-2, 2), rng.uniform(0.3, 2) * rng.choice([-1, 1])
    qN, x0 = rng.choice([0, 0.5, 1, 10]), rng.uniform(-5, 5)
    P, L = backward(a_, b_, qN)
    res = minimize(cost, np.zeros(4), args=(a_, b_, qN, x0), tol=1e-12)
    print(f"a={a_:+.3f} b={b_:+.3f} qN={qN:>4} x0={x0:+.3f} | p0 x0^2 = {P[0]*x0**2:10.5f} "
          f"| brute force = {res.fun:10.5f}")
    assert abs(P[0] * x0**2 - res.fun) < 1e-5 * max(1, res.fun)

print("\n=== (c) Monte Carlo: same gains, extra cost sum p_{k+1} sigma^2 ===")
a_, b_, qN, x0, sigma = 1.2, 1.0, 1.0, 10.0, 0.5
P, L = backward(a_, b_, qN)
c0 = sum(P[k + 1] for k in range(4)) * sigma**2
M = 400_000
X = np.full(M, x0)
Jmc = np.zeros(M)
for k in range(4):
    U = -L[k] * X
    Jmc += X**2 + U**2
    X = a_ * X + b_ * U + sigma * rng.standard_normal(M)
Jmc += qN * X**2
print(f"predicted J0 = p0 x0^2 + c0 = {P[0]*x0**2:.4f} + {c0:.4f} = {P[0]*x0**2 + c0:.4f}")
print(f"Monte Carlo  = {Jmc.mean():.4f} +- {2*Jmc.std()/np.sqrt(M):.4f}")
