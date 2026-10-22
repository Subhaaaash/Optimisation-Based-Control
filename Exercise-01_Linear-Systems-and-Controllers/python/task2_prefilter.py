"""
Optimization-Based Control Methods - Exercise 1, Task 2
State feedback with reference prefilter:  u_k = -K (x_k - F x_ref)

Limit:   x_inf = (I - (A - BK))^{-1} BK F x_ref
Zero error  <=>  (A - BK + BKF) x_ref = x_ref
1-D:     f = (1 - a + bk) / (bk)
General (rank B = n, BK invertible):  F = (BK)^{-1} (I - A) + I
"""
import numpy as np

np.set_printoptions(precision=5, suppress=True)


def sim(A, B, K, F, x_ref, N=200):
    x = np.zeros(A.shape[0])
    for _ in range(N):
        x = A @ x - B @ K @ (x - F @ x_ref)
    return x


print("--- 1-D: a = 0.5, b = 1, k = 0.5 (same plant as Task 1) ---")
a, b, k = 0.5, 1.0, 0.5
f = (1 - a + b * k) / (b * k)
print(f"f = (1 - a + bk)/(bk) = {f}")
for ff, label in [(1.0, "no prefilter (Task 1)"), (f, "with f")]:
    x = sim(np.array([[a]]), np.array([[b]]), np.array([[k]]), np.array([[ff]]), np.array([1.0]))
    print(f"  {label:22s}: x_inf = {x[0]:.6f}   (formula bkf/(1-a+bk) = {b*k*ff/(1-a+b*k):.6f})")

print("\n--- Sensitivity: design f for a = 0.5, but the true plant has a = 0.55 ---")
x = sim(np.array([[0.55]]), np.array([[b]]), np.array([[k]]), np.array([[f]]), np.array([1.0]))
print(f"  x_inf = {x[0]:.6f}  -> error {x[0]-1:+.4f}  (the prefilter needs the exact model)")

print("\n--- Matrix case: systems a) and b) of the numerical experiment ---")
B = np.array([[0.0, 1.0], [1.0, 0.0]])
I = np.eye(2)
x_ref = np.array([1.0, 0.0])
for name, A in [("a", np.array([[0.9, 0.5], [0.1, 0.3]])),
                ("b", np.array([[1.0, 0.5], [0.1, 0.7]]))]:
    K = np.linalg.inv(B) @ (A - 0.5 * I)
    BK = B @ K
    F = np.linalg.inv(BK) @ (I - A) + I
    print(f"\n system {name}):  BK =\n{BK}\n  det(BK) = {np.linalg.det(BK):.4f}"
          f"\n  (BK)^-1 =\n{np.linalg.inv(BK)}\n  I - A =\n{I - A}\n  F =\n{F}")
    print(f"  check A - BK + BKF = I :\n{A - BK + BK @ F}")
    for xr in (x_ref, np.array([-2.0, 3.0])):
        print(f"  x_ref = {xr}:  simulated x_inf = {sim(A, B, K, F, xr)}")
