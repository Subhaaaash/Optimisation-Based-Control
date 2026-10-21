"""
Optimization-Based Control Methods - Exercise 1, Task 0
Reachability / controllability and observability of a discrete-time LTI system

    x_{k+1} = A x_k + B u_k
    y_k     = C x_k
"""
import numpy as np

np.set_printoptions(precision=4, suppress=True)

A = np.array([[1, 0, 0],
              [0, 2, 2],
              [0, 1, 3]], dtype=float)
B = np.array([[1], [0], [1]], dtype=float)
C = np.array([[1, 1, 1]], dtype=float)
n = A.shape[0]


def rank(M):
    return np.linalg.matrix_rank(M)


def basis_of_column_space(M, tol=1e-10):
    """Orthonormal basis of the column space (image) of M, via SVD."""
    U, s, _ = np.linalg.svd(M)
    r = int(np.sum(s > tol))
    return U[:, :r]


def basis_of_null_space(M, tol=1e-10):
    """Orthonormal basis of the null space (kernel) of M, via SVD."""
    _, s, Vt = np.linalg.svd(M)
    r = int(np.sum(s > tol))
    return Vt[r:].T


# ---------------------------------------------------------------
# 0. Eigenvalues of A (needed for the controllability question)
# ---------------------------------------------------------------
eigvals = np.linalg.eigvals(A)
print("Eigenvalues of A:", np.round(np.real_if_close(eigvals), 4))
print()

# ---------------------------------------------------------------
# 1. Reachability: R = [B, AB, A^2 B]
# ---------------------------------------------------------------
R = np.hstack([np.linalg.matrix_power(A, i) @ B for i in range(n)])
print("Reachability matrix R = [B, AB, A^2 B]:\n", R)
print(f"rank(R) = {rank(R)}  (n = {n}) ->",
      "reachable" if rank(R) == n else "NOT reachable")

# PBH (Hautus) test: eigenvalue lambda is uncontrollable if rank[lambda I - A, B] < n
print("\nPBH test for each distinct eigenvalue:")
for lam in np.unique(np.round(np.real(eigvals), 8)):
    r = rank(np.hstack([lam * np.eye(n) - A, B]))
    print(f"  lambda = {lam:g}: rank[lambda I - A, B] = {r}",
          "-> uncontrollable mode" if r < n else "-> controllable mode")
print("Discrete time: controllable (to the origin) <=> every uncontrollable "
      "eigenvalue is 0.")
print()

# ---------------------------------------------------------------
# 2. Reachable subspace and target states
# ---------------------------------------------------------------
print("Basis of reachable subspace (column space of R):\n",
      basis_of_column_space(R))
normal = basis_of_null_space(R.T)      # vector orthogonal to the reachable plane
print("Normal vector of the reachable plane:", normal.ravel())
print()

targets = {"x_f1": np.array([[4], [2], [3]], float),
           "x_f2": np.array([[6], [2], [8]], float)}

for name, xf in targets.items():
    reachable = rank(np.hstack([R, xf])) == rank(R)
    print(f"{name} = {xf.ravel()}: reachable? {reachable}")
    if not reachable:
        continue
    # Try N = 1, 2, ..., n steps. After N steps:
    # x_N = A^{N-1} B u_0 + A^{N-2} B u_1 + ... + B u_{N-1}
    for N in range(1, n + 1):
        R_N = np.hstack([np.linalg.matrix_power(A, N - 1 - i) @ B
                         for i in range(N)])
        u, *_ = np.linalg.lstsq(R_N, xf, rcond=None)
        if np.allclose(R_N @ u, xf):
            print(f"  minimum number of steps N = {N}, inputs u = {u.ravel()}")
            # verify by simulation
            x = np.zeros((n, 1))
            for k in range(N):
                x = A @ x + B * u[k, 0]
                print(f"    x_{k+1} = {x.ravel()}")
            break
print()

# ---------------------------------------------------------------
# 3. Observability: O = [C; CA; CA^2]
# ---------------------------------------------------------------
O = np.vstack([C @ np.linalg.matrix_power(A, i) for i in range(n)])
print("Observability matrix O = [C; CA; CA^2]:\n", O)
print(f"rank(O) = {rank(O)}  (n = {n}) ->",
      "observable" if rank(O) == n else "NOT observable")

print("\nPBH test for each distinct eigenvalue:")
for lam in np.unique(np.round(np.real(eigvals), 8)):
    r = rank(np.vstack([lam * np.eye(n) - A, C]))
    print(f"  lambda = {lam:g}: rank[lambda I - A; C] = {r}",
          "-> unobservable mode" if r < n else "-> observable mode")
print()

# ---------------------------------------------------------------
# 4. Unobservable / observable subspace and initial states
# ---------------------------------------------------------------
N_O = basis_of_null_space(O)
print("Basis of UNobservable subspace (kernel of O):", N_O.ravel())
print("Basis of observable subspace (row space of O):\n",
      basis_of_column_space(O.T))
print()

inits = {"x0_1": np.array([[1], [2], [3]], float),
         "x0_2": np.array([[-2], [4], [-2]], float)}

for name, x0 in inits.items():
    unobs_part = N_O @ (N_O.T @ x0)          # projection onto kernel of O
    print(f"{name} = {x0.ravel()}")
    print(f"  output sequence y_0..y_2 = {(O @ x0).ravel()}")
    print(f"  component in unobservable subspace = {unobs_part.ravel()}")
    if np.allclose(unobs_part, 0):
        print("  -> lies in the observable subspace: can be observed")
    elif np.allclose(O @ x0, 0):
        print("  -> lies entirely in the unobservable subspace: output is zero, "
              "indistinguishable from x0 = 0, cannot be observed")
    else:
        print("  -> only its observable part can be reconstructed")
