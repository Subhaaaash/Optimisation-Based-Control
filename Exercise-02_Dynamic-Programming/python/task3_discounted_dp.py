"""
Optimization-Based Control Methods - Exercise 2, Task 3
Finite-horizon discounted DP and the current-value function V_k = alpha^{-k} J_k.

On a random finite MDP (states, actions, random transitions and costs):
  * J_k from the standard recursion with stage cost alpha^k g and terminal alpha^N g_N
  * V_k from  V_N = g_N,  V_k = min_u E[ g + alpha V_{k+1} ]
  * check V_k == alpha^{-k} J_k and that both give the same optimal actions
  * show that V_k depends only on the number of steps to go (time invariance):
    the same operator T is applied N-k times, so V_k(N) == V_{k+1}(N+1).

Run:  python task3_discounted_dp.py
"""
import numpy as np

rng = np.random.default_rng(1)
nS, nA, N, alpha = 6, 3, 8, 0.8

P = rng.random((nA, nS, nS))
P /= P.sum(axis=2, keepdims=True)          # P[u, x, x'] = Prob(x' | x, u)
g = rng.uniform(0, 10, size=(nS, nA))       # expected stage cost g(x, u)
gN = rng.uniform(0, 10, size=nS)


def J_recursion(N):
    J = [None] * (N + 1)
    mu = [None] * N
    J[N] = alpha**N * gN
    for k in range(N - 1, -1, -1):
        Q = alpha**k * g + np.einsum("uxy,y->xu", P, J[k + 1])
        J[k], mu[k] = Q.min(axis=1), Q.argmin(axis=1)
    return J, mu


def T(V):
    """Bellman operator of the time-invariant discounted problem."""
    Q = g + alpha * np.einsum("uxy,y->xu", P, V)
    return Q.min(axis=1), Q.argmin(axis=1)


def V_recursion(N):
    V = [None] * (N + 1)
    mu = [None] * N
    V[N] = gN.copy()
    for k in range(N - 1, -1, -1):
        V[k], mu[k] = T(V[k + 1])
    return V, mu


J, muJ = J_recursion(N)
V, muV = V_recursion(N)
for k in range(N + 1):
    assert np.allclose(V[k], alpha ** (-k) * J[k])
print("V_k == alpha^{-k} J_k for every k:  OK")
print("same optimal actions mu_k at every stage:", all(np.array_equal(a, b) for a, b in zip(muJ, muV)))

V_long, _ = V_recursion(N + 1)
print("V_k depends only on steps-to-go (V_k for N equals V_{k+1} for N+1):",
      all(np.allclose(V[k], V_long[k + 1]) for k in range(N + 1)))

print("\nmagnitudes at k = 0 vs k = N:")
print(f"  max |J_N| = {np.abs(J[N]).max():.4e}   (shrinks like alpha^N)")
print(f"  max |V_N| = {np.abs(V[N]).max():.4e}   (stays on the scale of g)")
alpha_small, N_big = 0.5, 1200
print(f"  with alpha={alpha_small}, N={N_big}: alpha^N = {alpha_small**N_big:.1e} (underflows to 0.0 in float64)")

# value iteration: V_k converges as the number of steps to go grows
Vk = gN.copy()
diffs = []
for _ in range(60):
    Vn, _ = T(Vk)
    diffs.append(np.abs(Vn - Vk).max())
    Vk = Vn
print("\n||V_{j+1} - V_j||_inf for j = 0..5:", np.round(diffs[:6], 4))
print("ratios (should be <= alpha = 0.8):", np.round(np.array(diffs[1:6]) / np.array(diffs[:5]), 4))
