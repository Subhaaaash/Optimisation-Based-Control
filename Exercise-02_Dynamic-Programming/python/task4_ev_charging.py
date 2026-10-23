"""
Optimization-Based Control Methods - Exercise 2, Task 4
Tabular DP with hard constraints: EV battery charging, N = 2.

  x_k in {0,1,2},  u_k in {0,1},  x_{k+1} = x_k + u_k must stay in {0,1,2}
  g0 = 2u + u^2 (night),  g1 = 4u + u^2 (morning peak),  g2(x2) = 5 (2 - x2)^2

Prints the DP tables J2, J1, mu1, J0(0), mu0(0) and checks them against a
brute-force enumeration of all admissible schedules.

Run:  python task4_ev_charging.py
"""
import itertools

X, U = (0, 1, 2), (0, 1)
g = {0: lambda u: 2 * u + u**2, 1: lambda u: 4 * u + u**2}
g2 = lambda x: 5 * (2 - x) ** 2

J2 = {x: g2(x) for x in X}
print("J2:", J2)

J1, mu1 = {}, {}
for x in X:
    q = {u: g[1](u) + J2[x + u] for u in U if x + u in X}
    J1[x] = min(q.values())
    mu1[x] = [u for u in q if q[u] == J1[x]]
    print(f"x1={x}: " + ", ".join(f"u={u}: {g[1](u)} + {J2[x+u]} = {v}" for u, v in q.items())
          + f"  ->  J1={J1[x]}, mu1 in {mu1[x]}")

q0 = {u: g[0](u) + J1[0 + u] for u in U}
J0 = min(q0.values())
mu0 = [u for u in q0 if q0[u] == J0]
print("x0=0: " + ", ".join(f"u={u}: {g[0](u)} + {J1[u]} = {v}" for u, v in q0.items())
      + f"  ->  J0(0)={J0}, mu0 in {mu0}")

print("\nbrute force over all admissible (u0, u1):")
best = []
for u0, u1 in itertools.product(U, U):
    x1 = 0 + u0
    x2 = x1 + u1
    if x2 not in X:
        continue
    c = g[0](u0) + g[1](u1) + g2(x2)
    print(f"  (u0,u1)=({u0},{u1}) -> x2={x2}, cost {g[0](u0)} + {g[1](u1)} + {g2(x2)} = {c}")
    best.append((c, (u0, u1), x2))
cmin = min(b[0] for b in best)
print("optimal:", [(s, x2) for c, s, x2 in best if c == cmin], "cost", cmin)
assert cmin == J0 == 8
