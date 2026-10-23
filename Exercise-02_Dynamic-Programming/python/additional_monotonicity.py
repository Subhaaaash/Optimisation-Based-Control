"""
Optimization-Based Control Methods - Exercise 2, Additional task (monotonicity of DP)

For a time-invariant finite MDP:  if J_{N-1} <= J_N pointwise, then J_k <= J_{k+1} for all k.
The script builds random problems that contain a free "stay" action (x' = x at zero cost),
which guarantees J_{N-1}(x) <= 0 + J_N(x) = J_N(x); it then runs the DP and checks the
chain J_0 <= J_1 <= ... <= J_N.
It also shows a case with J_{N-1} > J_N somewhere, where the reverse ordering appears.

Run:  python additional_monotonicity.py
"""
import numpy as np

rng = np.random.default_rng(3)


def dp(P, g, gN, N):
    J = [None] * (N + 1)
    J[N] = gN
    for k in range(N - 1, -1, -1):
        J[k] = (g + np.einsum("uxy,y->xu", P, J[k + 1])).min(axis=1)   # T applied to J_{k+1}
    return J


ok = 0
for trial in range(200):
    nS, nA, N = rng.integers(3, 8), rng.integers(2, 4), rng.integers(2, 10)
    P = rng.random((nA, nS, nS)); P /= P.sum(axis=2, keepdims=True)
    P[0] = np.eye(nS)                      # action 0 = "stay where you are" ...
    g = rng.uniform(-2, 5, (nS, nA))
    g[:, 0] = 0.0                          # ... at zero cost
    gN = rng.uniform(0, 40, nS)
    J = dp(P, g, gN, N)
    if np.all(J[N - 1] <= J[N]):           # hypothesis of the task
        assert all(np.all(J[k] <= J[k + 1] + 1e-12) for k in range(N))
        ok += 1
print(f"{ok} random problems satisfied J_(N-1) <= J_N; in all of them J_0 <= J_1 <= ... <= J_N held.")

# example with zero terminal cost: J_{N-1} >= J_N = 0, and the chain flips to non-increasing
nS, nA, N = 4, 2, 5
P = rng.random((nA, nS, nS)); P /= P.sum(axis=2, keepdims=True)
g = rng.uniform(0, 5, (nS, nA))
J = dp(P, g, np.zeros(nS), N)
print("\nwith g_N = 0 (so J_(N-1) >= J_N), stage-wise values at state 0:",
      np.round([J[k][0] for k in range(N + 1)], 3), "-> non-increasing in k")
