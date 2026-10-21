---
title: "Exercise 1 — Step-by-Step Solutions"
subtitle: "Optimization-Based Control Methods · Linear Systems and Linear Controllers"
---

# Exercise 1 — Step-by-Step Solutions

**Course:** Optimization-Based Control Methods (TU Chemnitz, Automatic Control Group, Dr. Adrian Redder)
**Covers:** eigenvalue practice, Task 0 (all parts), Task 1 (all parts)
**Status:** Tasks 2, 3 and the numerical experiment will be added later

Every step is written out. Section references like "Concepts §5.4" point to [`01_Concepts.md`](01_Concepts.md), which explains the theory behind each step. The Python scripts in [`python/`](python/) verify every number.

---

## Contents

- [Part A — Eigenvalue practice](#part-a--eigenvalue-practice)
- [Part B — Task 0: Controllable and observable subspaces](#part-b--task-0-controllable-and-observable-subspaces)
- [Part C — Task 1: Tracking error for a state-feedback controller](#part-c--task-1-tracking-error-for-a-state-feedback-controller)

---

## Part A — Eigenvalue practice

**Method used throughout (Concepts §2.3):**
(1) form $\lambda I - A$; (2) take $\det(\lambda I - A)$; (3) find the roots; (4) check $\sum\lambda = \operatorname{tr}A$ and $\prod\lambda = \det A$.

**Stability rule (Concepts §4.1):** $x_{k+1} = Ax_k$ is asymptotically stable $\iff$ every eigenvalue has $|\lambda| < 1$.

### Problem 1

$$
A = \begin{bmatrix} 4 & 1 \\ 2 & 3 \end{bmatrix}
$$

**Step 1.** $\lambda I - A = \begin{bmatrix} \lambda - 4 & -1 \\ -2 & \lambda - 3 \end{bmatrix}$.

**Step 2.** $\det(\lambda I - A) = (\lambda - 4)(\lambda - 3) - (-1)(-2) = \lambda^2 - 7\lambda + 12 - 2 = \lambda^2 - 7\lambda + 10$.

*Shortcut check:* $\operatorname{tr}A = 7$ and $\det A = 4\cdot 3 - 1\cdot 2 = 10$ give the same polynomial $\lambda^2 - 7\lambda + 10$.

**Step 3.** Look for two numbers with product $10$ and sum $-7$: these are $-2$ and $-5$. So $(\lambda - 2)(\lambda - 5) = 0$.

$$
\boxed{\lambda = 2,\; 5}
$$

**Step 4.** $2 + 5 = 7 = \operatorname{tr}A$ ✓, $2\cdot 5 = 10 = \det A$ ✓.

**Stable?** No: $|2| > 1$ and $|5| > 1$.

### Problem 2

$$
A = \begin{bmatrix} 3 & 5 \\ 0 & -1 \end{bmatrix}
$$

$\lambda I - A = \begin{bmatrix} \lambda - 3 & -5 \\ 0 & \lambda + 1 \end{bmatrix}$, so
$\det = (\lambda - 3)(\lambda + 1) - (-5)(0) = (\lambda - 3)(\lambda + 1)$.

The $0$ below the diagonal removes the off-diagonal product. This is why **triangular matrices have their diagonal entries as eigenvalues**.

$$
\boxed{\lambda = 3,\; -1}
$$

Check: $3 + (-1) = 2 = \operatorname{tr}A$ ✓, $3\cdot(-1) = -3 = \det A$ ✓.
**Stable?** No: $|3| > 1$. ($|-1| = 1$ is also not strictly inside the unit circle.)

### Problem 3

$$
A = \begin{bmatrix} 0 & 1 \\ -2 & -3 \end{bmatrix}
$$

$\lambda I - A = \begin{bmatrix} \lambda & -1 \\ 2 & \lambda + 3 \end{bmatrix}$, so
$\det = \lambda(\lambda + 3) - (-1)(2) = \lambda^2 + 3\lambda + 2$.

Look for two numbers with product $2$ and sum $3$: these are $1$ and $2$. So $(\lambda + 1)(\lambda + 2) = 0$.

$$
\boxed{\lambda = -1,\; -2}
$$

Check: $-1 - 2 = -3 = \operatorname{tr}A$ ✓, $(-1)(-2) = 2 = \det A$ ✓.

*Note:* in this **controllable canonical form**, the bottom row $[-2,\,-3]$ holds the coefficients of $\lambda^2 + 3\lambda + 2$ with the signs flipped.

**Stable?** No: $|-2| = 2 > 1$ (and $|-1| = 1$ is on the boundary).

### Problem 4

$$
A = \begin{bmatrix} 1 & -2 \\ 2 & 1 \end{bmatrix}
$$

$\lambda I - A = \begin{bmatrix} \lambda - 1 & 2 \\ -2 & \lambda - 1 \end{bmatrix}$, so
$\det = (\lambda - 1)^2 - (2)(-2) = (\lambda - 1)^2 + 4$.

Set it to zero: $(\lambda - 1)^2 = -4 \Rightarrow \lambda - 1 = \pm\sqrt{-4} = \pm 2i$.

(Same result with the formula: $\operatorname{tr}A = 2$, $\det A = 1 + 4 = 5$, $\lambda = \frac{2 \pm \sqrt{4 - 20}}{2} = \frac{2 \pm 4i}{2}$.)

$$
\boxed{\lambda = 1 \pm 2i}
$$

Check: $(1+2i) + (1-2i) = 2$ ✓, $(1+2i)(1-2i) = 1 - 4i^2 = 1 + 4 = 5$ ✓.
Magnitude: $|\lambda| = \sqrt{1^2 + 2^2} = \sqrt 5 \approx 2.24$. Complex eigenvalues mean rotation, and $|\lambda| > 1$ means the state spirals outward.
**Stable?** No.

### Problem 5

$$
A = \begin{bmatrix} 0.5 & 0.2 \\ 0.1 & 0.4 \end{bmatrix}
$$

$\lambda I - A = \begin{bmatrix} \lambda - 0.5 & -0.2 \\ -0.1 & \lambda - 0.4 \end{bmatrix}$, so
$\det = (\lambda - 0.5)(\lambda - 0.4) - (0.2)(0.1) = \lambda^2 - 0.9\lambda + 0.20 - 0.02 = \lambda^2 - 0.9\lambda + 0.18$.

Quadratic formula: discriminant $= 0.9^2 - 4(0.18) = 0.81 - 0.72 = 0.09$, and $\sqrt{0.09} = 0.3$:

$$
\lambda = \frac{0.9 \pm 0.3}{2}
\quad\Rightarrow\quad
\boxed{\lambda = 0.6,\; 0.3}
$$

Check: $0.6 + 0.3 = 0.9$ ✓, $0.6 \cdot 0.3 = 0.18$ ✓.
**Stable?** **Yes**, since $|0.6| < 1$ and $|0.3| < 1$. This is the only stable matrix among the eight problems.

### Problem 6

$$
A = \begin{bmatrix} 1 & 0 & 0 \\ 4 & -2 & 0 \\ 7 & 5 & 0.5 \end{bmatrix}
$$

**Step 1.** Put $\lambda$ on the diagonal and flip every sign. The diagonal entries become $\lambda - 1$, $\lambda - (-2) = \lambda + 2$ and $\lambda - 0.5$:

$$
\lambda I - A = \begin{bmatrix} \lambda - 1 & 0 & 0 \\ -4 & \lambda + 2 & 0 \\ -7 & -5 & \lambda - 0.5 \end{bmatrix}.
$$

**Step 2.** Expand along row 1, $(\lambda - 1,\; 0,\; 0)$. The two zeros remove their terms:

$$
\det(\lambda I - A) = +(\lambda - 1)\det\begin{bmatrix} \lambda + 2 & 0 \\ -5 & \lambda - 0.5 \end{bmatrix} - 0\cdot(\dots) + 0\cdot(\dots).
$$

The $2\times 2$ determinant is $(\lambda + 2)(\lambda - 0.5) - (0)(-5) = (\lambda + 2)(\lambda - 0.5)$. Therefore

$$
\det(\lambda I - A) = (\lambda - 1)(\lambda + 2)(\lambda - 0.5).
$$

**Step 3.** Set each factor to zero:

$$
\boxed{\lambda = 1,\; -2,\; 0.5}
$$

**Step 4.** Sum: $1 - 2 + 0.5 = -0.5 = \operatorname{tr}A$ ✓. Product: $1\cdot(-2)\cdot 0.5 = -1 = \det A$ (the product of the diagonal for a triangular matrix) ✓.

**Why the shortcut works:** the entries below the diagonal ($4, 7, 5$) never contributed, because each time they would appear they are multiplied by a $0$ from above the diagonal.

**Stable?** No. $|-2| = 2 > 1$: that mode doubles in size and flips sign every step. $\lambda = 1$ never shrinks either. Only the $\lambda = 0.5$ mode dies out.

### Problem 7

$$
A = \begin{bmatrix} 2 & 0 & 0 \\ 0 & 1 & 4 \\ 0 & 2 & 3 \end{bmatrix}
$$

$$
\lambda I - A = \begin{bmatrix} \lambda - 2 & 0 & 0 \\ 0 & \lambda - 1 & -4 \\ 0 & -2 & \lambda - 3 \end{bmatrix}.
$$

Expand along row 1:

$$
\det = (\lambda - 2)\big[(\lambda - 1)(\lambda - 3) - (-4)(-2)\big] = (\lambda - 2)\big[\lambda^2 - 4\lambda + 3 - 8\big] = (\lambda - 2)(\lambda^2 - 4\lambda - 5).
$$

Factor: two numbers with product $-5$ and sum $-4$ are $-5$ and $+1$, so $\lambda^2 - 4\lambda - 5 = (\lambda - 5)(\lambda + 1)$.

$$
\boxed{\lambda = 2,\; 5,\; -1}
$$

Check: $2 + 5 - 1 = 6 = 2 + 1 + 3$ ✓. $2\cdot 5\cdot(-1) = -10$, and $\det A = 2\,(1\cdot 3 - 4\cdot 2) = 2(-5) = -10$ ✓.
**Stable?** No.

### Problem 8

$$
A = \begin{bmatrix} 2 & 0 & 1 \\ 0 & 3 & 0 \\ 1 & 0 & 2 \end{bmatrix}
$$

**Step 1.**

$$
\lambda I - A = \begin{bmatrix} \lambda - 2 & 0 & -1 \\ 0 & \lambda - 3 & 0 \\ -1 & 0 & \lambda - 2 \end{bmatrix}.
$$

**Step 2: choose the row.** The middle row $(0,\; \lambda - 3,\; 0)$ has two zeros, so only one term survives.

**Step 3: sign.** By the checkerboard $\begin{bmatrix} + & - & + \\ - & + & - \\ + & - & + \end{bmatrix}$, position $(2,2)$ has sign $+$.

**Step 4: expand.** Delete row 2 and column 2. What remains is $\begin{bmatrix} \lambda - 2 & -1 \\ -1 & \lambda - 2 \end{bmatrix}$:

$$
\det(\lambda I - A) = +(\lambda - 3)\det\begin{bmatrix} \lambda - 2 & -1 \\ -1 & \lambda - 2 \end{bmatrix}.
$$

**Step 5: the $2\times 2$ determinant.** $(\lambda - 2)(\lambda - 2) - (-1)(-1) = (\lambda - 2)^2 - 1$. Careful: $(-1)(-1) = +1$, and it is *subtracted*.

**Step 6: factor.** Use the difference of squares $a^2 - b^2 = (a - b)(a + b)$ with $a = \lambda - 2$, $b = 1$:

$$
(\lambda - 2)^2 - 1 = (\lambda - 3)(\lambda - 1).
$$

(Or expand: $\lambda^2 - 4\lambda + 4 - 1 = \lambda^2 - 4\lambda + 3 = (\lambda - 3)(\lambda - 1)$.)

$$
\det(\lambda I - A) = (\lambda - 3)^2(\lambda - 1)
\quad\Rightarrow\quad
\boxed{\lambda = 1,\; 3,\; 3}
$$

**Step 7: check.** $1 + 3 + 3 = 7 = 2 + 3 + 2$ ✓. $1\cdot 3\cdot 3 = 9$, and $\det A = 3\,(2\cdot 2 - 1\cdot 1) = 9$ (expanding along the middle row) ✓.

**Same result expanding along row 1** (signs $+,-,+$):

$$
\begin{aligned}
\det &= +(\lambda - 2)\det\begin{bmatrix} \lambda - 3 & 0 \\ 0 & \lambda - 2 \end{bmatrix} - 0 + (-1)\det\begin{bmatrix} 0 & \lambda - 3 \\ -1 & 0 \end{bmatrix} \\
&= (\lambda - 2)(\lambda - 3)(\lambda - 2) + (-1)\big[0\cdot 0 - (\lambda - 3)(-1)\big] \\
&= (\lambda - 3)(\lambda - 2)^2 - (\lambda - 3) = (\lambda - 3)\big[(\lambda - 2)^2 - 1\big].
\end{aligned}
$$

It is the same expression, but with more work. That is why choosing the row with the most zeros pays off.

**Stable?** No: $|3| > 1$, so that mode triples every step.

### Summary of Part A

| # | Eigenvalues | Stable? |
|---|---|---|
| 1 | $2,\;5$ | no |
| 2 | $3,\;-1$ | no |
| 3 | $-1,\;-2$ | no |
| 4 | $1 \pm 2i$ ($\lvert\lambda\rvert = \sqrt5$) | no |
| 5 | $0.6,\;0.3$ | **yes** |
| 6 | $1,\;-2,\;0.5$ | no |
| 7 | $2,\;5,\;-1$ | no |
| 8 | $1,\;3,\;3$ | no |

---

## Part B — Task 0: Controllable and observable subspaces

**System:**

$$
x_{k+1} = \underbrace{\begin{bmatrix} 1 & 0 & 0 \\ 0 & 2 & 2 \\ 0 & 1 & 3 \end{bmatrix}}_{A} x_k + \underbrace{\begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}}_{B} u_k,
\qquad
y_k = \underbrace{\begin{bmatrix} 1 & 1 & 1 \end{bmatrix}}_{C} x_k .
$$

Here $n = 3$ states, $m = 1$ input and $p = 1$ output.

### B.0 Preparation: eigenvalues and eigenvectors of $A$

From Concepts §2.7: $\det(\lambda I - A) = (\lambda - 1)^2(\lambda - 4)$, so $\lambda = 1, 1, 4$.

**Eigenvector for $\lambda = 4$:** solve $(4I - A)v = 0$:

$$
\begin{bmatrix} 3 & 0 & 0 \\ 0 & 2 & -2 \\ 0 & -1 & 1 \end{bmatrix} v = 0
\;\Rightarrow\; 3v_1 = 0,\; 2v_2 - 2v_3 = 0
\;\Rightarrow\; v = \begin{bmatrix} 0 \\ 1 \\ 1 \end{bmatrix}.
$$

Check: $A[0,1,1]^\top = [0,\,2+2,\,1+3]^\top = [0,4,4]^\top = 4\cdot[0,1,1]^\top$ ✓.

**Eigenvectors for $\lambda = 1$:** solve $(I - A)v = 0$:

$$
\begin{bmatrix} 0 & 0 & 0 \\ 0 & -1 & -2 \\ 0 & -1 & -2 \end{bmatrix} v = 0
\;\Rightarrow\; -v_2 - 2v_3 = 0
\;\Rightarrow\; v_2 = -2v_3,\; v_1 \text{ free}.
$$

This gives two independent eigenvectors, $[1,0,0]^\top$ and $[0,-2,1]^\top$, so the eigenspace is a plane. In particular, $[1,-2,1]^\top = [1,0,0]^\top + [0,-2,1]^\top$ is also an eigenvector for $\lambda = 1$. It will reappear in B.4.

---

### B.1 Task 0.1 — Is the system reachable? Is it controllable?

**Step 1: build $W_c = [\,B,\; AB,\; A^2B\,]$** (Concepts §5.1, $n = 3$, so we need up to $A^2B$).

$$
B = \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}.
$$

$$
AB = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 2 & 2 \\ 0 & 1 & 3 \end{bmatrix}\begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}
= \begin{bmatrix} 1\cdot1 + 0 + 0 \\ 0 + 2\cdot 0 + 2\cdot 1 \\ 0 + 1\cdot 0 + 3\cdot 1 \end{bmatrix}
= \begin{bmatrix} 1 \\ 2 \\ 3 \end{bmatrix}.
$$

$$
A^2B = A(AB) = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 2 & 2 \\ 0 & 1 & 3 \end{bmatrix}\begin{bmatrix} 1 \\ 2 \\ 3 \end{bmatrix}
= \begin{bmatrix} 1 \\ 2\cdot2 + 2\cdot3 \\ 1\cdot 2 + 3\cdot 3 \end{bmatrix}
= \begin{bmatrix} 1 \\ 10 \\ 11 \end{bmatrix}.
$$

$$
W_c = \begin{bmatrix} 1 & 1 & 1 \\ 0 & 2 & 10 \\ 1 & 3 & 11 \end{bmatrix}.
$$

**Step 2: rank by row reduction.**

$$
\begin{bmatrix} 1 & 1 & 1 \\ 0 & 2 & 10 \\ 1 & 3 & 11 \end{bmatrix}
\xrightarrow{R_3 - R_1}
\begin{bmatrix} 1 & 1 & 1 \\ 0 & 2 & 10 \\ 0 & 2 & 10 \end{bmatrix}
\xrightarrow{R_3 - R_2}
\begin{bmatrix} 1 & 1 & 1 \\ 0 & 2 & 10 \\ 0 & 0 & 0 \end{bmatrix}
$$

$$
\operatorname{rank}(W_c) = 2 < 3 = n .
$$

Cross-check with the determinant: $\det W_c = 1(2\cdot11 - 10\cdot3) - 1(0\cdot 11 - 10\cdot 1) + 1(0\cdot 3 - 2\cdot 1) = -8 + 10 - 2 = 0$ ✓.

The dependence is explicit: $A^2B = 5\,AB - 4\,B$, since $5[1,2,3] - 4[1,0,1] = [1,10,11]$. This is Cayley–Hamilton at work. With $p(\lambda) = \lambda^3 - 6\lambda^2 + 9\lambda - 4$, the minimal relation for this $B$ happens to come from the factor $(\lambda-1)(\lambda-4) = \lambda^2 - 5\lambda + 4$.

> **Conclusion:** the system is **not reachable** (not fully state controllable).

**Step 3: PBH test — which eigenvalue is the problem?** (Concepts §5.6)

*For $\lambda = 4$:*

$$
\begin{bmatrix} 4I - A & B \end{bmatrix} = \begin{bmatrix} 3 & 0 & 0 & 1 \\ 0 & 2 & -2 & 0 \\ 0 & -1 & 1 & 1 \end{bmatrix}.
$$

Row 1 is the only row with a nonzero first entry, so it is independent of the others. Rows 2 and 3 are not multiples of each other: $R_3 = -\tfrac12 R_2$ would need the last entry to be $0$, but it is $1$. So the rank is $3$, and **$\lambda = 4$ is controllable**.

*For $\lambda = 1$:*

$$
\begin{bmatrix} I - A & B \end{bmatrix} = \begin{bmatrix} 0 & 0 & 0 & 1 \\ 0 & -1 & -2 & 0 \\ 0 & -1 & -2 & 1 \end{bmatrix}.
$$

$R_3 = R_1 + R_2$, so the rank is $2 < 3$. **$\lambda = 1$ is uncontrollable.**

**Step 4: controllable in the sense of Def. 1.2?** Controllability (steering every $x_0$ to the origin) would require every uncontrollable mode to vanish by itself in finite time, which needs eigenvalue $0$. Here the uncontrollable eigenvalue is $\lambda = 1$. That mode neither grows nor decays, and the input cannot influence it.

> **Conclusion:** the system is **not controllable**. It is not even **stabilizable**: the uncontrollable eigenvalue $1$ is not strictly inside the unit circle, so no $K$ can make $A - BK$ asymptotically stable.

**Step 5: intuition — a conserved quantity** (Concepts §5.5). Take $w = [-1,-1,1]^\top$. Then
$w^\top A = -R_1 - R_2 + R_3 = [-1,\,-2+1,\,-2+3] = [-1,-1,1] = 1\cdot w^\top$, and
$w^\top B = -1 - 0 + 1 = 0$. So $z_k := w^\top x_k = x_{3,k} - x_{1,k} - x_{2,k}$ satisfies

$$
z_{k+1} = w^\top A x_k + w^\top B u_k = z_k + 0 .
$$

**The quantity $x_3 - x_1 - x_2$ never changes, whatever input is applied.** A system starting with $z_0 \neq 0$ can never reach the origin, where $z = 0$. That is exactly why it is not controllable.

**Python check** (`python/task0_reachability_observability.py`):

```text
Eigenvalues of A: [4. 1. 1.]
Reachability matrix R = [B, AB, A^2 B]:
 [[ 1.  1.  1.]
 [ 0.  2. 10.]
 [ 1.  3. 11.]]
rank(R) = 2  (n = 3) -> NOT reachable
PBH test for each distinct eigenvalue:
  lambda = 1: rank[lambda I - A, B] = 2 -> uncontrollable mode
  lambda = 4: rank[lambda I - A, B] = 3 -> controllable mode
```

---

### B.2 Task 0.2 — Reachable subspace and target states

**Step 1: reachable subspace** $= \operatorname{im}(W_c)$ (Concepts §5.4). Since $A^2B$ depends on the first two columns,

$$
\mathcal{R} = \operatorname{span}\left\{ \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}, \begin{bmatrix} 1 \\ 2 \\ 3 \end{bmatrix} \right\},
$$

which is a plane through the origin (dimension $2$ = rank).

**Step 2: describe the plane with one equation** using the cross product (Concepts §1.6):

$$
\begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix} \times \begin{bmatrix} 1 \\ 2 \\ 3 \end{bmatrix}
= \begin{bmatrix} 0\cdot 3 - 1\cdot 2 \\ 1\cdot 1 - 1\cdot 3 \\ 1\cdot 2 - 0\cdot 1 \end{bmatrix}
= \begin{bmatrix} -2 \\ -2 \\ 2 \end{bmatrix} \;\propto\; \begin{bmatrix} -1 \\ -1 \\ 1 \end{bmatrix}.
$$

So the plane is $-x_1 - x_2 + x_3 = 0$:

$$
\boxed{\;\mathcal{R} = \{\,x \in \mathbb{R}^3 : x_3 = x_1 + x_2\,\}\;}
$$

Check: $B$: $1 = 1 + 0$ ✓; $AB$: $3 = 1 + 2$ ✓; $A^2B$: $11 = 1 + 10$ ✓.
This normal vector is exactly the conserved-quantity vector $w$ from B.1 Step 5. Starting from $x_0 = 0$ we have $z = 0$ forever, so we stay in the plane $x_3 - x_1 - x_2 = 0$.

**Step 3: target $x_f^{(1)} = [4, 2, 3]^\top$.**
$x_1 + x_2 = 4 + 2 = 6 \neq 3 = x_3$. The target is not in $\mathcal{R}$.

> **$x_f^{(1)}$ is not reachable** from the origin, with any number of steps or any inputs.

(Rank test: $\operatorname{rank}[W_c \;\; x_f^{(1)}] = 3 \neq 2$.)

**Step 4: target $x_f^{(2)} = [6, 2, 8]^\top$.**
$x_1 + x_2 = 6 + 2 = 8 = x_3$ ✓. **Reachable.**

**Step 5: minimum number of steps.** From $x_0 = 0$, after $N$ steps $x_N = \sum_{j=0}^{N-1}A^{N-1-j}Bu_j$.

*Try $N = 1$:* $x_1 = Bu_0 = [u_0,\; 0,\; u_0]^\top$. The second component is always $0$, but we need $2$. **Impossible.**

*Try $N = 2$:* $x_2 = AB\,u_0 + B\,u_1$:

$$
u_0\begin{bmatrix} 1 \\ 2 \\ 3 \end{bmatrix} + u_1\begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix} = \begin{bmatrix} 6 \\ 2 \\ 8 \end{bmatrix}
\quad\Longleftrightarrow\quad
\begin{cases}
u_0 + u_1 = 6 \\
2u_0 = 2 \\
3u_0 + u_1 = 8
\end{cases}
$$

- Row 2: $u_0 = 1$.
- Row 1: $u_1 = 6 - 1 = 5$.
- Row 3 (consistency check): $3\cdot 1 + 5 = 8$ ✓.

> **$x_f^{(2)}$ is reachable in a minimum of $N = 2$ steps, with $u_0 = 1,\; u_1 = 5$.**
> (The solution is unique, because $B$ and $AB$ are independent.)

**Step 6: verify by simulation.**

$$
x_1 = A\cdot 0 + B\cdot 1 = \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix},
\qquad
x_2 = Ax_1 + B\cdot 5 = \begin{bmatrix} 1 \\ 2 \\ 3 \end{bmatrix} + \begin{bmatrix} 5 \\ 0 \\ 5 \end{bmatrix} = \begin{bmatrix} 6 \\ 2 \\ 8 \end{bmatrix} \checkmark
$$

**Python check:**

```text
x_f1 = [4. 2. 3.]: reachable? False
x_f2 = [6. 2. 8.]: reachable? True
  minimum number of steps N = 2, inputs u = [1. 5.]
    x_1 = [1. 0. 1.]
    x_2 = [6. 2. 8.]
```

---

### B.3 Task 0.3 — Is the system observable?

**Step 1: build $W_o = [\,C;\; CA;\; CA^2\,]$** (Concepts §6.1). A row vector times $A$ is a weighted sum of the rows of $A$.

$$
C = \begin{bmatrix} 1 & 1 & 1 \end{bmatrix}.
$$

$$
CA = 1\cdot[1,0,0] + 1\cdot[0,2,2] + 1\cdot[0,1,3] = \begin{bmatrix} 1 & 3 & 5 \end{bmatrix}.
$$

$$
CA^2 = (CA)A = 1\cdot[1,0,0] + 3\cdot[0,2,2] + 5\cdot[0,1,3] = [1,\; 6+5,\; 6+15] = \begin{bmatrix} 1 & 11 & 21 \end{bmatrix}.
$$

$$
W_o = \begin{bmatrix} 1 & 1 & 1 \\ 1 & 3 & 5 \\ 1 & 11 & 21 \end{bmatrix}.
$$

**Step 2: rank.**

$$
\xrightarrow{R_2 - R_1,\;R_3 - R_1}
\begin{bmatrix} 1 & 1 & 1 \\ 0 & 2 & 4 \\ 0 & 10 & 20 \end{bmatrix}
\xrightarrow{R_3 - 5R_2}
\begin{bmatrix} 1 & 1 & 1 \\ 0 & 2 & 4 \\ 0 & 0 & 0 \end{bmatrix}
$$

$$
\operatorname{rank}(W_o) = 2 < 3 .
$$

(Again $CA^2 = 5\,CA - 4\,C$: $5[1,3,5] - 4[1,1,1] = [1,11,21]$ ✓.)

> **Conclusion:** the system is **not observable**.

**Step 3: PBH test.**

*$\lambda = 4$:* $\begin{bmatrix} 4I - A \\ C \end{bmatrix} = \begin{bmatrix} 3 & 0 & 0 \\ 0 & 2 & -2 \\ 0 & -1 & 1 \\ 1 & 1 & 1 \end{bmatrix}$. Rows 1, 2 and 4 give $\det\begin{bmatrix} 3&0&0\\0&2&-2\\1&1&1\end{bmatrix} = 3(2\cdot1 - (-2)\cdot1) = 12 \neq 0$, so the rank is $3$. **Observable mode.**

*$\lambda = 1$:* $\begin{bmatrix} I - A \\ C \end{bmatrix} = \begin{bmatrix} 0 & 0 & 0 \\ 0 & -1 & -2 \\ 0 & -1 & -2 \\ 1 & 1 & 1 \end{bmatrix}$. Rows 2 and 3 are equal and row 1 is zero, so only $[0,-1,-2]$ and $[1,1,1]$ are independent. The rank is $2 < 3$. **Unobservable mode.**

Equivalently, the $\lambda = 1$ eigenvector $v = [1,-2,1]^\top$ has $Cv = 1 - 2 + 1 = 0$: the sensor cannot see it.
The system is also **not detectable**, because the unobservable eigenvalue $1$ is not inside the unit circle.

**Python check:**

```text
Observability matrix O = [C; CA; CA^2]:
 [[ 1.  1.  1.]
 [ 1.  3.  5.]
 [ 1. 11. 21.]]
rank(O) = 2  (n = 3) -> NOT observable
  lambda = 1: rank[lambda I - A; C] = 2 -> unobservable mode
  lambda = 4: rank[lambda I - A; C] = 3 -> observable mode
```

---

### B.4 Task 0.4 — Observable subspace and initial states

**Step 1: unobservable subspace $= \ker(W_o)$.** Solve $W_o v = 0$. Use the two independent rows:

$$
\begin{aligned}
v_1 + v_2 + v_3 &= 0 \qquad (1)\\
v_1 + 3v_2 + 5v_3 &= 0 \qquad (2)
\end{aligned}
$$

$(2) - (1)$: $2v_2 + 4v_3 = 0 \Rightarrow v_2 = -2v_3$.
Into $(1)$: $v_1 = -v_2 - v_3 = 2v_3 - v_3 = v_3$.
With $v_3 = t$: $v = t\,[1, -2, 1]^\top$.

Check with row 3: $1 - 22 + 21 = 0$ ✓.

$$
\boxed{\;\ker(W_o) = \operatorname{span}\{[1,\,-2,\,1]^\top\}\;}\quad\text{(a line, dimension } 3 - 2 = 1\text{)}
$$

**Step 2: observable subspace** $= \ker(W_o)^\perp$ = row space of $W_o$:

$$
\boxed{\;\mathcal{O} = \{\,x : x_1 - 2x_2 + x_3 = 0\,\} = \operatorname{span}\{[1,1,1]^\top,\;[1,3,5]^\top\}\;}
$$

Check: $[1,1,1]$: $1 - 2 + 1 = 0$ ✓; $[1,3,5]$: $1 - 6 + 5 = 0$ ✓.

**Step 3: initial state $x_0^{(1)} = [1, 2, 3]^\top$.**

- Test for membership in $\mathcal{O}$: $1 - 2\cdot2 + 3 = 0$ ✓. So $x_0^{(1)}$ lies in the observable subspace, and its unobservable part is zero.
- Its outputs: $\vec y = W_o x_0^{(1)} = [1+2+3,\; 1+6+15,\; 1+22+63]^\top = [6,\; 22,\; 86]^\top$, which is nonzero.

> **$x_0^{(1)}$ can be observed**: it is the unique state in the observable subspace that produces the outputs $6, 22, 86$.
> *Fine print:* since the system is not observable, any $x_0^{(1)} + c\,[1,-2,1]^\top$ would produce the same outputs. The measurements identify the observable component, and $x_0^{(1)}$ consists only of that component.

**Step 4: initial state $x_0^{(2)} = [-2, 4, -2]^\top$.**

- $x_0^{(2)} = -2\,[1, -2, 1]^\top$, so it lies entirely in $\ker(W_o)$.
- Outputs: $y_0 = -2 + 4 - 2 = 0$; $y_1 = -2 + 12 - 10 = 0$; $y_2 = -2 + 44 - 42 = 0$.
- Because $[1,-2,1]$ is an eigenvector for $\lambda = 1$, the state never moves ($x_k = x_0^{(2)}$ for all $k$), and $Cx_k = 0$ forever.

> **$x_0^{(2)}$ cannot be observed**: it produces $y_k = 0$ for all $k$, exactly like $x_0 = 0$.

**Python check:**

```text
Basis of UNobservable subspace (kernel of O): [-0.4082  0.8165 -0.4082]   (= multiple of [1,-2,1])
x0_1 = [1. 2. 3.]
  output sequence y_0..y_2 = [ 6. 22. 86.]
  -> lies in the observable subspace: can be observed
x0_2 = [-2.  4. -2.]
  output sequence y_0..y_2 = [0. 0. 0.]
  -> lies entirely in the unobservable subspace: cannot be observed
```

### Task 0 — summary

| Question | Answer |
|---|---|
| Eigenvalues of $A$ | $1, 1, 4$ |
| $\operatorname{rank} W_c$ | $2$ — **not reachable** |
| Uncontrollable mode (PBH) | $\lambda = 1$ — **not controllable**, not stabilizable |
| Reachable subspace | plane $x_3 = x_1 + x_2$ |
| $x_f^{(1)} = [4,2,3]$ | not reachable |
| $x_f^{(2)} = [6,2,8]$ | reachable, min. $N = 2$, $(u_0,u_1) = (1,5)$ |
| $\operatorname{rank} W_o$ | $2$ — **not observable** |
| Unobservable mode (PBH) | $\lambda = 1$ — not detectable |
| Unobservable subspace | $\operatorname{span}\{[1,-2,1]\}$ |
| Observable subspace | plane $x_1 - 2x_2 + x_3 = 0$ |
| $x_0^{(1)} = [1,2,3]$ | observable ($y = 6, 22, 86$) |
| $x_0^{(2)} = [-2,4,-2]$ | not observable ($y \equiv 0$) |

---

## Part C — Task 1: Tracking error for a state-feedback controller

**Setting:** $x_{k+1} = Ax_k + Bu_k$, controller $u_k = -K(x_k - x_{\text{ref}})$ with constant $x_{\text{ref}}$, and $K$ such that $A - BK$ is asymptotically stable.

### C.0 Setting up the closed loop

Substitute the controller into the system:

$$
\begin{aligned}
x_{k+1} &= Ax_k + B\big(-K(x_k - x_{\text{ref}})\big) \\
&= Ax_k - BKx_k + BKx_{\text{ref}} \\
&= \underbrace{(A - BK)}_{=:M}\,x_k + BKx_{\text{ref}} .
\end{aligned}
$$

So the closed loop is $x_{k+1} = Mx_k + c$ with constant $c = BKx_{\text{ref}}$, and by assumption all eigenvalues of $M$ satisfy $|\lambda| < 1$.

---

### C.1 Task 1.1 — $I - (A - BK)$ is invertible and $x_k \to x_\infty$

#### Part (a): invertibility of $I - M$

**Proof by contradiction.**

1. Suppose $I - M$ is **not** invertible.
2. Then there is a vector $v \neq 0$ with $(I - M)v = 0$ (Concepts §1.7).
3. Expand: $v - Mv = 0$, i.e. $Mv = v = 1\cdot v$.
4. So $v$ is an eigenvector of $M$ with eigenvalue $\lambda = 1$.
5. But $M$ is asymptotically stable, so every eigenvalue satisfies $|\lambda| < 1$. Since $|1| = 1$, this is a contradiction.
6. Therefore $I - M = I - (A - BK)$ **is invertible**. $\blacksquare$

> **Remember:** $I - M$ is singular $\iff$ $1$ is an eigenvalue of $M$.

#### Part (b): convergence to $x_\infty = (I - (A - BK))^{-1}BKx_{\text{ref}}$

**Step 1: find the equilibrium.** A limit must be an equilibrium (Concepts §4.3), so set $x_{k+1} = x_k = x_\infty$:

$$
x_\infty = Mx_\infty + BKx_{\text{ref}}
\;\Longrightarrow\;
x_\infty - Mx_\infty = BKx_{\text{ref}}
\;\Longrightarrow\;
(I - M)\,x_\infty = BKx_{\text{ref}} .
$$

By part (a), $I - M$ is invertible, so there is exactly one equilibrium:

$$
x_\infty = (I - M)^{-1}BKx_{\text{ref}} = \big(I - (A - BK)\big)^{-1}BKx_{\text{ref}} .
$$

**Step 2: show that the trajectory actually converges to it** (shift trick, Concepts §4.4). Define $d_k := x_k - x_\infty$:

$$
\begin{aligned}
d_{k+1} &= x_{k+1} - x_\infty \\
&= \big(Mx_k + BKx_{\text{ref}}\big) - \big(Mx_\infty + BKx_{\text{ref}}\big) &&\text{(both satisfy the same update)}\\
&= M(x_k - x_\infty) \\
&= Md_k .
\end{aligned}
$$

By induction, $d_k = M^k d_0$. Since $M$ is asymptotically stable, $M^k \to 0$, hence $d_k \to 0$ for **every** $x_0$:

$$
\lim_{k\to\infty} x_k = x_\infty = \big(I - (A - BK)\big)^{-1}BKx_{\text{ref}} . \qquad\blacksquare
$$

#### Alternative proof (via the lecture's Proposition 1.1)

The closed loop is a system with matrix $M$ and constant "input" $BKx_{\text{ref}}$. By the solution formula,

$$
x_k = M^kx_0 + \sum_{j=0}^{k-1}M^{k-1-j}BKx_{\text{ref}} = M^kx_0 + \Big(\sum_{i=0}^{k-1}M^i\Big)BKx_{\text{ref}} .
$$

As $k\to\infty$: $M^kx_0 \to 0$, and by the geometric series (Concepts §4.2) $\sum_{i=0}^{\infty}M^i = (I - M)^{-1}$. This gives the same $x_\infty$.

---

### C.2 Task 1.2 — Asymptotic tracking error and when it can be zero

#### Step 1: the one-dimensional case (as the hint suggests)

Let $A = a$, $B = b$, $K = k$, all scalars. The closed loop is $x_{k+1} = (a - bk)x_k + bk\,x_{\text{ref}}$.

- **Stability:** $|a - bk| < 1$.
- **Limit:** from C.1, $x_\infty = \dfrac{bk}{1 - (a - bk)}\,x_{\text{ref}} = \dfrac{bk}{1 - a + bk}\,x_{\text{ref}}$.
- **Error:** put both terms over a common denominator:

$$
\begin{aligned}
x_\infty - x_{\text{ref}}
&= \frac{bk}{1 - a + bk}x_{\text{ref}} - \frac{1 - a + bk}{1 - a + bk}x_{\text{ref}} \\
&= \frac{bk - 1 + a - bk}{1 - a + bk}\,x_{\text{ref}} \\
&= \boxed{\frac{a - 1}{1 - a + bk}\,x_{\text{ref}}}
\end{aligned}
$$

**Observation:** the numerator $a - 1$ **does not contain $k$**. Therefore:

- If $a = 1$: the error is $0$ for every stabilizing $k$ ✓.
- If $a \neq 1$ (and $x_{\text{ref}} \neq 0$): the error is **never** zero, whatever $k$ is.

**Can the error at least be made small?** A larger $bk$ makes the denominator larger, but stability restricts it: $|a - bk| < 1 \iff a - 1 < bk < a + 1$. So $1 - a + bk < 2$, and the error is bounded below in size by $\frac{|a - 1|}{2}|x_{\text{ref}}|$. It cannot be made arbitrarily small.

**Numerical example:** $a = 0.5$, $b = 1$. Stability requires $bk \in (-0.5,\, 1.5)$.

| $k$ | $a - bk$ (closed-loop pole) | $x_\infty / x_{\text{ref}} = \frac{bk}{1-a+bk}$ | error $/x_{\text{ref}}$ |
|---|---|---|---|
| $0.5$ | $0$ | $0.5/1.0 = 0.5$ | $-0.5$ |
| $1.0$ | $-0.5$ | $1.0/1.5 \approx 0.667$ | $\approx -0.333$ |
| $1.4$ | $-0.9$ | $1.4/1.9 \approx 0.737$ | $\approx -0.263$ |
| $\to 1.5$ | $\to -1$ (stability limit) | $\to 0.75$ | $\to -0.25$ |

You never get closer than 75 % of the way to the target.

#### Step 2: the general (matrix) case

Use the trick $x_{\text{ref}} = (I - M)^{-1}(I - M)\,x_{\text{ref}}$, which is valid because $I - M$ is invertible:

$$
\begin{aligned}
x_\infty - x_{\text{ref}}
&= (I - M)^{-1}BKx_{\text{ref}} - (I - M)^{-1}(I - M)x_{\text{ref}} \\
&= (I - M)^{-1}\big[BK - (I - M)\big]x_{\text{ref}} \\
&= (I - M)^{-1}\big[BK - I + A - BK\big]x_{\text{ref}} &&(M = A - BK)\\
&= (I - M)^{-1}(A - I)\,x_{\text{ref}} .
\end{aligned}
$$

$$
\boxed{\;x_\infty - x_{\text{ref}} = \big(I - A + BK\big)^{-1}(A - I)\,x_{\text{ref}}\;}
$$

For $n = 1$ this reduces to the scalar formula above.

#### Step 3: when is the error zero?

$(I - A + BK)^{-1}$ is invertible, so it maps only the zero vector to zero. Therefore

$$
x_\infty - x_{\text{ref}} = 0
\iff (A - I)\,x_{\text{ref}} = 0
\iff A\,x_{\text{ref}} = x_{\text{ref}} .
$$

> **Answer to Task 1.2.**
> The asymptotic tracking error is $x_\infty - x_{\text{ref}} = (I - A + BK)^{-1}(A - I)x_{\text{ref}}$.
> It is zero **if and only if $Ax_{\text{ref}} = x_{\text{ref}}$**, i.e. $x_{\text{ref}}$ is already an equilibrium of the *uncontrolled* system ($u = 0$), or equivalently $x_{\text{ref}} \in \ker(A - I)$ (an eigenvector of $A$ for eigenvalue $1$, or $x_{\text{ref}} = 0$).
> - If this holds, **any** stabilizing $K$ achieves $x_\infty = x_{\text{ref}}$.
> - If it does not hold, **no** choice of $K$ achieves zero error. The condition does not depend on $K$ at all.
> - Zero error for **every** reference requires $A = I$ (in 1-D: $a = 1$).

#### Step 4: interpretation (why this happens)

To hold the state at $x_{\text{ref}}$, the input must satisfy $x_{\text{ref}} = Ax_{\text{ref}} + B\bar u$, i.e. $B\bar u = (I - A)x_{\text{ref}}$. This input is nonzero unless $Ax_{\text{ref}} = x_{\text{ref}}$. But at $x = x_{\text{ref}}$ the P-controller outputs $u = -K\cdot 0 = 0$. The controller can only produce input by keeping an error, so the system settles where the error is just large enough to supply the needed input. This is the classic **steady-state error of a P-controller** (Concepts §7.3).

This motivates Task 2 (rescale the reference with a prefilter $F$) and Task 3 (add an integrator, i.e. a PI controller).

#### Step 5: numerical verification (`python/task1_tracking_error.py`)

| Case | Simulated $x_\infty$ | Formula |
|---|---|---|
| 1-D: $a = 0.5,\;b = 1,\;k = 0.5,\;x_{\text{ref}} = 1$ | $0.5$ | $\frac{0.5}{1.0} = 0.5$ ✓ |
| System a) of the numerical experiment, $K = \begin{bmatrix} 0 & 0.2 \\ 0.6 & 0 \end{bmatrix}$, $x_{\text{ref}} = [1,0]^\top$ | $[0.931,\;0.103]$ | $[0.931,\;0.103]$ ✓ |
| Error for the same case | $[-0.069,\;0.103]$ | $(I-A+BK)^{-1}(A-I)x_{\text{ref}} = [-0.069,\;0.103]$ ✓ |

For system a), the closed-loop eigenvalues are $0.445$ and $-0.045$ (stable), but $A[1,0]^\top = [0.9,0.1]^\top \neq [1,0]^\top$, so a steady-state error remains, as the theory predicts.
