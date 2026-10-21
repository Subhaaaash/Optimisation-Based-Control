"""
Optimization-Based Control Methods - Exercise 1, Task 1
Tracking error of the state-feedback controller u_k = -K (x_k - x_ref)

Closed loop:  x_{k+1} = (A - BK) x_k + BK x_ref
Limit:        x_inf   = (I - (A - BK))^{-1} BK x_ref
Error:        x_inf - x_ref = (I - A + BK)^{-1} (A - I) x_ref
"""
import numpy as np

np.set_printoptions(precision=4, suppress=True)


def simulate(A, B, K, x_ref, x0, steps=200):
    x = x0.copy()
    for _ in range(steps):
        x = A @ x - B @ K @ (x - x_ref)
    return x


def analyse(name, A, B, K, x_ref):
    n = A.shape[0]
    I = np.eye(n)
    M = A - B @ K
    eig = np.linalg.eigvals(M)
    print(f"--- {name} ---")
    print("eig(A - BK) =", eig, "| stable:", bool(np.all(np.abs(eig) < 1)))
    x_sim = simulate(A, B, K, x_ref, np.zeros(n))
    x_formula = np.linalg.solve(I - M, B @ K @ x_ref)
    err_formula = np.linalg.solve(I - M, (A - I) @ x_ref)
    print("x_inf simulated :", x_sim)
    print("x_inf formula   :", x_formula)
    print("error simulated :", x_sim - x_ref)
    print("error formula   :", err_formula)
    print("A x_ref == x_ref ? ", np.allclose(A @ x_ref, x_ref), "\n")


# 1-D example: a = 0.5, b = 1, k = 0.5
analyse("1-D: a=0.5, b=1, k=0.5",
        np.array([[0.5]]), np.array([[1.0]]), np.array([[0.5]]), np.array([1.0]))

# 1-D sweep over k: the error never vanishes because a != 1
a, b = 0.5, 1.0
print("--- 1-D sweep (a=0.5, b=1) ---")
for k in [0.5, 1.0, 1.4, 1.49]:
    print(f"k={k:<5} pole={a-b*k:+.2f}  x_inf/x_ref={b*k/(1-a+b*k):.3f}  "
          f"error/x_ref={(a-1)/(1-a+b*k):+.3f}")
print()

# 1-D example with a = 1: zero error for any stabilizing k
analyse("1-D: a=1, b=1, k=0.5 (zero error expected)",
        np.array([[1.0]]), np.array([[1.0]]), np.array([[0.5]]), np.array([1.0]))

# System a) of the numerical experiment
A = np.array([[0.9, 0.5], [0.1, 0.3]])
B = np.array([[0, 1], [1, 0.0]])
K = np.array([[0, 0.2], [0.6, 0]])
analyse("System a), x_ref = [1, 0]", A, B, K, np.array([1.0, 0.0]))
