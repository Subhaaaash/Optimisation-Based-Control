---
title: "Exercise 1 — Concepts"
subtitle: "Optimization-Based Control Methods · Linear Systems and Linear Controllers"
---

# Exercise 1 — Concepts

**Course:** Optimization-Based Control Methods (TU Chemnitz, Automatic Control Group, Dr. Adrian Redder)
**Topic:** Linear systems and linear controllers
**Covers:** everything needed for Tasks 0–3 and the numerical experiment

This document explains every idea from the beginning, assuming no previous control-systems course. Every formula is derived, and no steps are skipped. The worked solutions are in [`02_Solutions.md`](02_Solutions.md).

Notation follows the lecture script: $W_c$ is the controllability (reachability) matrix and $W_o$ is the observability matrix.

---

## Contents

1. [Linear algebra toolkit](#1-linear-algebra-toolkit)
2. [Eigenvalues and eigenvectors](#2-eigenvalues-and-eigenvectors)
3. [Discrete-time LTI systems](#3-discrete-time-lti-systems)
4. [Stability](#4-stability)
5. [Reachability, controllability, stabilizability](#5-reachability-controllability-stabilizability)
6. [Observability and detectability](#6-observability-and-detectability)
7. [State feedback, equilibria and tracking](#7-state-feedback-equilibria-and-tracking)
8. [Cheat sheet](#8-cheat-sheet)

---

## 1. Linear algebra toolkit

Control theory for linear systems is mostly linear algebra. This section collects every tool used later.

### 1.1 Vectors and matrices

A **vector** $x \in \mathbb{R}^n$ is a column of $n$ real numbers:

$$
x = \begin{bmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{bmatrix}.
$$

A **matrix** $A \in \mathbb{R}^{n \times m}$ has $n$ rows and $m$ columns. The entry in row $i$, column $j$ is $a_{ij}$.

The **transpose** $A^\top$ swaps rows and columns: $(A^\top)_{ij} = a_{ji}$. A column vector transposed becomes a row vector.

### 1.2 Matrix–vector and matrix–matrix products

**Rule:** entry $i$ of $Ax$ is the dot product of row $i$ of $A$ with $x$:

$$
(Ax)_i = \sum_{j=1}^{m} a_{ij}\, x_j .
$$

**Example.**

$$
\begin{bmatrix} 1 & 0 & 0 \\ 0 & 2 & 2 \\ 0 & 1 & 3 \end{bmatrix}
\begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}
=
\begin{bmatrix} 1\cdot 1 + 0\cdot 0 + 0\cdot 1 \\ 0\cdot 1 + 2\cdot 0 + 2\cdot 1 \\ 0\cdot 1 + 1\cdot 0 + 3\cdot 1 \end{bmatrix}
=
\begin{bmatrix} 1 \\ 2 \\ 3 \end{bmatrix}.
$$

**Equivalent view (very useful):** $Ax$ is a *weighted sum of the columns of $A$*, with the weights being the entries of $x$:

$$
Ax = x_1 \,(\text{column }1) + x_2\,(\text{column }2) + \dots + x_m\,(\text{column }m).
$$

**Row vector times matrix:** $c^\top A$ is a weighted sum of the *rows* of $A$. For example $[1,1,1]\,A$ is the sum of all rows of $A$.

**Matrix–matrix product:** $(AB)_{ij}$ is row $i$ of $A$ dotted with column $j$ of $B$. In general $AB \neq BA$.

**Matrix powers:** $A^0 = I$, $A^1 = A$, $A^2 = AA$, $A^k = A A^{k-1}$. To compute $A^2B$, it is easiest to compute $AB$ first and then $A(AB)$.

**Identity matrix** $I$: ones on the diagonal, zeros elsewhere, so $Ix = x$ for every $x$.

### 1.3 Linear combinations, span, linear independence

A **linear combination** of vectors $v_1,\dots,v_k$ is $c_1 v_1 + \dots + c_k v_k$ with real numbers $c_i$.

The **span** of $v_1,\dots,v_k$ is the set of *all* their linear combinations. In $\mathbb{R}^3$:

- the span of one nonzero vector is a **line** through the origin;
- the span of two non-parallel vectors is a **plane** through the origin;
- the span of three independent vectors is **all of $\mathbb{R}^3$**.

The vectors are **linearly independent** if the only way to get $c_1 v_1 + \dots + c_k v_k = 0$ is $c_1 = \dots = c_k = 0$. Otherwise, one of them can be written as a combination of the others, and it adds nothing new to the span.

### 1.4 Rank

The **rank** of a matrix is the number of linearly independent columns. It always equals the number of linearly independent rows. For $M \in \mathbb{R}^{n\times m}$: $\operatorname{rank}(M) \le \min(n,m)$. The matrix has **full rank** when equality holds.

**How to compute the rank by hand (row reduction / Gaussian elimination).**
Allowed operations, which never change the rank:

1. swap two rows;
2. multiply a row by a nonzero number;
3. add a multiple of one row to another row.

Use them to create zeros below the leading entries ("echelon form"). The rank is the number of nonzero rows left.

**Example.**

$$
\begin{bmatrix} 1 & 1 & 1 \\ 0 & 2 & 10 \\ 1 & 3 & 11 \end{bmatrix}
\xrightarrow{R_3 \leftarrow R_3 - R_1}
\begin{bmatrix} 1 & 1 & 1 \\ 0 & 2 & 10 \\ 0 & 2 & 10 \end{bmatrix}
\xrightarrow{R_3 \leftarrow R_3 - R_2}
\begin{bmatrix} 1 & 1 & 1 \\ 0 & 2 & 10 \\ 0 & 0 & 0 \end{bmatrix}
$$

There are two nonzero rows, so the rank is $2$.

**Another test for square matrices:** an $n\times n$ matrix has rank $n$ if and only if its determinant is nonzero (Section 1.7).

### 1.5 Column space (image) and null space (kernel)

For a matrix $M$:

- **Column space / image:** $\operatorname{im}(M) = \{ M v : v \text{ any vector}\}$, the span of the columns. Its dimension is $\operatorname{rank}(M)$.
- **Null space / kernel:** $\ker(M) = \{ v : Mv = 0 \}$, all vectors that $M$ sends to zero.

**Rank–nullity theorem.** For $M$ with $m$ columns:

$$
\operatorname{rank}(M) + \dim \ker(M) = m .
$$

**How to find a kernel:** solve $Mv = 0$ by row reduction, express the "pivot" variables in terms of the free variables, and read off a basis.

**Example.** $\ker\begin{bmatrix}1&1&1\\1&3&5\end{bmatrix}$:

$$
\begin{aligned}
v_1 + v_2 + v_3 &= 0 \\
v_1 + 3v_2 + 5v_3 &= 0
\end{aligned}
\quad\xrightarrow{\text{subtract}}\quad
2v_2 + 4v_3 = 0 \;\Rightarrow\; v_2 = -2v_3 .
$$

Then $v_1 = -v_2 - v_3 = 2v_3 - v_3 = v_3$. With $v_3 = t$: $v = t\,[1,\,-2,\,1]^\top$. The kernel is the line spanned by $[1,-2,1]^\top$.

**How to test if a vector $w$ lies in $\operatorname{im}(M)$:** $w \in \operatorname{im}(M)$ if and only if $\operatorname{rank}([M \;\; w]) = \operatorname{rank}(M)$. Appending $w$ must not add a new independent direction.

### 1.6 Orthogonality, orthogonal complement, cross product

Two vectors are **orthogonal** if their dot product is zero: $v^\top w = \sum_i v_i w_i = 0$.

The **orthogonal complement** $S^\perp$ of a subspace $S$ is the set of all vectors orthogonal to every vector in $S$. In $\mathbb{R}^3$, the complement of a line is a plane and vice versa.

A key fact: **the row space of $M$ is the orthogonal complement of $\ker(M)$**, because $Mv = 0$ says exactly that $v$ is orthogonal to every row of $M$.

**Describing a plane through the origin in $\mathbb{R}^3$ by one equation.**
If the plane is spanned by $p$ and $q$, a **normal vector** is the cross product

$$
p \times q = \begin{bmatrix} p_2 q_3 - p_3 q_2 \\ p_3 q_1 - p_1 q_3 \\ p_1 q_2 - p_2 q_1 \end{bmatrix},
$$

and the plane is $\{x : (p\times q)^\top x = 0\}$.

### 1.7 Determinants

The determinant is a number $\det(M)$ attached to a square matrix. The facts we need:

- $\det(M) \neq 0 \iff M$ is invertible $\iff \operatorname{rank}(M) = n \iff \ker(M) = \{0\}$.
- $\det(M) = 0 \iff$ there is a $v \neq 0$ with $Mv = 0$.

**$2\times 2$:**

$$
\det\begin{bmatrix} a & b \\ c & d \end{bmatrix} = ad - bc .
$$

**$3\times 3$ — cofactor (Laplace) expansion along a row.** Each position has a sign from the checkerboard

$$
\begin{bmatrix} + & - & + \\ - & + & - \\ + & - & + \end{bmatrix},
\qquad \text{sign of position } (i,j) = (-1)^{i+j}.
$$

Pick any row. For each entry of that row: take the entry, multiply by its checkerboard sign, multiply by the $2\times 2$ determinant that remains after deleting that entry's row and column (the **minor**). Add the three results. Along the first row:

$$
\det\begin{bmatrix} a & b & c \\ d & e & f \\ g & h & i \end{bmatrix}
= +a\det\begin{bmatrix} e & f \\ h & i \end{bmatrix}
 - b\det\begin{bmatrix} d & f \\ g & i \end{bmatrix}
 + c\det\begin{bmatrix} d & e \\ g & h \end{bmatrix}.
$$

**Rules of thumb.**

- Any row or column gives the same answer, so **expand along the one with the most zeros**: each zero kills a whole term.
- **Triangular matrix** (all zeros above *or* all zeros below the diagonal): the determinant is the product of the diagonal entries.
- $\det(AB) = \det(A)\det(B)$, $\det(A^\top) = \det(A)$.

### 1.8 Matrix inverse

$M^{-1}$ is the matrix with $M M^{-1} = M^{-1} M = I$. It exists exactly when $\det M \neq 0$. Then $Mx = c$ has the unique solution $x = M^{-1} c$.

**$2\times 2$ formula:**

$$
\begin{bmatrix} a & b \\ c & d \end{bmatrix}^{-1}
= \frac{1}{ad - bc}\begin{bmatrix} d & -b \\ -c & a \end{bmatrix}.
$$

For a scalar, $m^{-1} = 1/m$, which exists when $m \neq 0$.

---

## 2. Eigenvalues and eigenvectors

### 2.1 Definition and meaning

A nonzero vector $v$ is an **eigenvector** of a square matrix $A$, with **eigenvalue** $\lambda$, if

$$
A v = \lambda v .
$$

Most vectors change direction when multiplied by $A$. Eigenvectors don't: $A$ only stretches them by the factor $\lambda$. If $\lambda$ is negative, the vector also flips. If $\lambda$ is complex, the state rotates.

**Why control engineers care.** If $x_0 = v$ and $x_{k+1} = Ax_k$, then

$$
x_1 = Av = \lambda v,\quad x_2 = A(\lambda v) = \lambda^2 v,\quad \dots,\quad x_k = \lambda^k v .
$$

So along an eigenvector, the state is simply multiplied by $\lambda$ every step. Each eigenvalue is a **mode** of the system. The lecture also calls the eigenvalues of $A$ the **poles** of the system.

### 2.2 Deriving the characteristic equation

Start from $Av = \lambda v$ and rewrite it:

$$
Av = \lambda v
\;\iff\; Av = \lambda I v
\;\iff\; \lambda I v - A v = 0
\;\iff\; (\lambda I - A)\, v = 0 .
$$

(We wrote $\lambda v = \lambda I v$ so that we can factor out $v$. "$\lambda - A$" alone would not make sense, since a number minus a matrix is undefined.)

We need a **nonzero** $v$ solving $(\lambda I - A)v = 0$. By Section 1.7, such a $v$ exists exactly when $\lambda I - A$ is singular:

$$
\boxed{\;\det(\lambda I - A) = 0\;}
\qquad\text{(characteristic equation)}
$$

$p(\lambda) = \det(\lambda I - A)$ is a polynomial of degree $n$ in $\lambda$, called the **characteristic polynomial**. Its $n$ roots (counted with multiplicity, possibly complex) are the eigenvalues.

### 2.3 Procedure

1. **Form $\lambda I - A$:** put $\lambda$ on the diagonal and flip the sign of *every* entry of $A$. The diagonal becomes $\lambda - a_{ii}$ and an off-diagonal entry $a_{ij}$ becomes $-a_{ij}$.
2. **Take the determinant** (expand along the row or column with the most zeros).
3. **Find the roots** of the polynomial (factor it, or use the quadratic formula).
4. **Check** with trace and determinant (Section 2.5).
5. *(If needed)* **eigenvectors:** for each $\lambda$, solve $(\lambda I - A)v = 0$, which is a kernel computation as in Section 1.5.

### 2.4 The $2\times 2$ shortcut

For $A = \begin{bmatrix} a & b \\ c & d\end{bmatrix}$:

$$
\det(\lambda I - A)
= \det\begin{bmatrix} \lambda - a & -b \\ -c & \lambda - d \end{bmatrix}
= (\lambda - a)(\lambda - d) - (-b)(-c)
= \lambda^2 - (a+d)\lambda + (ad - bc).
$$

Since $a + d = \operatorname{tr}(A)$ (the trace, i.e. the sum of the diagonal entries) and $ad - bc = \det(A)$:

$$
\boxed{\;\lambda^2 - \operatorname{tr}(A)\,\lambda + \det(A) = 0\;}
\qquad\Rightarrow\qquad
\lambda_{1,2} = \frac{\operatorname{tr}(A) \pm \sqrt{\operatorname{tr}(A)^2 - 4\det(A)}}{2}.
$$

### 2.5 Checks that always hold (any size)

$$
\lambda_1 + \lambda_2 + \dots + \lambda_n = \operatorname{tr}(A),
\qquad
\lambda_1 \cdot \lambda_2 \cdots \lambda_n = \det(A).
$$

Use them every time to catch arithmetic mistakes.

### 2.6 Special structures (shortcuts)

- **Diagonal or triangular matrix:** the eigenvalues are the diagonal entries.
  *Why:* $\lambda I - A$ is also triangular, so its determinant is the product of its diagonal: $(\lambda - a_{11})(\lambda - a_{22})\cdots(\lambda - a_{nn})$. The entries off the diagonal are always multiplied by a zero during expansion.
- **Block structure:** if row $i$ and column $i$ are zero except for the diagonal entry $a_{ii}$, then $a_{ii}$ is an eigenvalue. The remaining eigenvalues are those of the smaller matrix left after deleting row $i$ and column $i$. More generally, for block-triangular matrices, the eigenvalues are the eigenvalues of the diagonal blocks.
- **Complex eigenvalues:** if $\operatorname{tr}(A)^2 - 4\det(A) < 0$, the roots are $\lambda = \alpha \pm i\beta$. Their magnitude is $|\lambda| = \sqrt{\alpha^2 + \beta^2}$. Complex eigenvalues mean rotation or oscillation.
- **Repeated eigenvalues:** a root that appears more than once, e.g. $(\lambda - 3)^2$. Its **algebraic multiplicity** is how often it is a root. Its **geometric multiplicity** is $\dim\ker(\lambda I - A)$, the number of independent eigenvectors.

### 2.7 Worked example

$$
A = \begin{bmatrix} 1 & 0 & 0 \\ 0 & 2 & 2 \\ 0 & 1 & 3 \end{bmatrix},
\qquad
\lambda I - A = \begin{bmatrix} \lambda - 1 & 0 & 0 \\ 0 & \lambda - 2 & -2 \\ 0 & -1 & \lambda - 3 \end{bmatrix}.
$$

Expand along row 1, which has two zeros:

$$
\det(\lambda I - A) = (\lambda - 1)\det\begin{bmatrix} \lambda - 2 & -2 \\ -1 & \lambda - 3 \end{bmatrix}
= (\lambda-1)\big[(\lambda-2)(\lambda-3) - (-2)(-1)\big].
$$

Expand the bracket: $(\lambda - 2)(\lambda - 3) = \lambda^2 - 5\lambda + 6$, and $(-2)(-1) = 2$, so the bracket is $\lambda^2 - 5\lambda + 4 = (\lambda - 1)(\lambda - 4)$. Therefore

$$
\det(\lambda I - A) = (\lambda - 1)^2(\lambda - 4)
\quad\Rightarrow\quad
\lambda = 1,\,1,\,4 .
$$

Check: $1 + 1 + 4 = 6 = 1 + 2 + 3 = \operatorname{tr}(A)$ ✓ and $1\cdot 1\cdot 4 = 4 = \det(A)$ ✓.

### 2.8 Cayley–Hamilton theorem

**Every square matrix satisfies its own characteristic polynomial:** if $p(\lambda) = \lambda^n + c_{n-1}\lambda^{n-1} + \dots + c_0$, then

$$
p(A) = A^n + c_{n-1}A^{n-1} + \dots + c_1 A + c_0 I = 0 .
$$

**Consequence:** $A^n = -c_{n-1}A^{n-1} - \dots - c_0 I$, so $A^n$, and therefore every higher power, is a linear combination of $I, A, \dots, A^{n-1}$. We use this in Section 5 to show that $n$ steps are always enough.

*Example:* for the $A$ above, $p(\lambda) = (\lambda-1)^2(\lambda-4) = \lambda^3 - 6\lambda^2 + 9\lambda - 4$, so $A^3 = 6A^2 - 9A + 4I$.

---

## 3. Discrete-time LTI systems

### 3.1 The model

$$
\boxed{\;
\begin{aligned}
x_{k+1} &= A x_k + B u_k \\
y_k &= C x_k \;(+\, D u_k)
\end{aligned}
\;}
$$

| Symbol | Name | Size | Meaning |
|---|---|---|---|
| $k$ | time index | $k = 0,1,2,\dots$ | sample number; real time is $t = k\,T_s$ with sampling time $T_s$ |
| $x_k$ | state | $n$ | everything you need to know about the system at time $k$ |
| $u_k$ | input | $m$ | what the controller commands (voltage, force, …) |
| $y_k$ | output | $p$ | what the sensors measure |
| $A$ | system matrix | $n\times n$ | how the state evolves by itself |
| $B$ | input matrix | $n\times m$ | how inputs enter the state |
| $C$ | output matrix | $p\times n$ | which combination of states the sensors see |
| $D$ | feedthrough | $p\times m$ | direct input-to-output path (often $0$) |

- **Linear:** the right-hand sides are linear in $x_k$ and $u_k$. So superposition holds: the response to a sum of inputs or initial states is the sum of the individual responses.
- **Time-invariant (LTI):** $A, B, C, D$ are constant and do not depend on $k$.
- **Single-input (SISO)** if $m = 1$, **multiple-input (MIMO)** if $m > 1$. The number of inputs is the number of columns of $B$.

### 3.2 Solution formula (derived by unrolling)

Apply the update equation repeatedly:

$$
\begin{aligned}
x_1 &= A x_0 + B u_0 \\
x_2 &= A x_1 + B u_1 = A(Ax_0 + Bu_0) + Bu_1 = A^2 x_0 + ABu_0 + Bu_1 \\
x_3 &= A x_2 + B u_2 = A^3 x_0 + A^2Bu_0 + ABu_1 + Bu_2 .
\end{aligned}
$$

The pattern: input $u_j$ is applied at time $j$ and then gets multiplied by $A$ once for every remaining step until time $k$, which is $k - 1 - j$ times. Hence

$$
\boxed{\; x_k = \underbrace{A^k x_0}_{\text{free response}} + \underbrace{\sum_{j=0}^{k-1} A^{k-1-j} B u_j}_{\text{forced response}} \;}
$$

(Formally: true for $k=0$; if true for $k$, then $x_{k+1} = A x_k + B u_k = A^{k+1}x_0 + \sum_{j=0}^{k-1}A^{k-j}Bu_j + Bu_k$, which is the formula for $k+1$.)

### 3.3 Modes: why eigenvalues decide everything

Suppose $A$ has $n$ independent eigenvectors $v_1,\dots,v_n$ with eigenvalues $\lambda_1,\dots,\lambda_n$ (this holds, for example, when all eigenvalues are distinct). Write the initial state as a combination $x_0 = c_1 v_1 + \dots + c_n v_n$. With $u = 0$, using $A^k v_i = \lambda_i^k v_i$:

$$
x_k = A^k x_0 = c_1 \lambda_1^k v_1 + \dots + c_n \lambda_n^k v_n .
$$

Each term behaves like $\lambda_i^k$:

| $\lambda$ | $\lambda^k$ as $k\to\infty$ |
|---|---|
| $\lvert\lambda\rvert < 1$ | $\to 0$ (decays) |
| $\lvert\lambda\rvert = 1$ | stays the same size ($\lambda = 1$: constant; $\lambda=-1$: flips sign; complex: rotates) |
| $\lvert\lambda\rvert > 1$ | grows without bound |
| $\lambda = 0$ | is exactly $0$ after one step |

> In discrete time, the dividing line is the **unit circle** $|\lambda| = 1$. In continuous time it would instead be the imaginary axis, $\operatorname{Re}\lambda = 0$. Don't mix the two.

---

## 4. Stability

### 4.1 Definitions (for $x_{k+1} = Mx_k$)

- **Asymptotically stable:** $x_k \to 0$ for every $x_0$.
  $\iff$ **all eigenvalues satisfy $|\lambda| < 1$** (strictly inside the unit circle).
  $\iff M^k \to 0$ as $k \to \infty$.
- **Lyapunov (marginally) stable:** $x_k$ stays bounded for every $x_0$, but need not go to $0$.
  This holds when all $|\lambda| \le 1$ and every eigenvalue on the unit circle has enough independent eigenvectors (geometric multiplicity = algebraic multiplicity). The lecture phrases this as: eigenvalues inside the unit circle, with at least one of absolute value 1.
- **Unstable:** some $x_0$ gives an unbounded $x_k$. This happens, for example, if any $|\lambda| > 1$.

In the lecture, "the matrix $M$ is asymptotically stable" is shorthand for "the trajectories of $x_{k+1} = Mx_k$ are asymptotically stable".

### 4.2 The geometric (Neumann) series and Proposition 1.1

**Lemma.** If all eigenvalues of $M$ satisfy $|\lambda|<1$, then $I - M$ is invertible and

$$
\sum_{j=0}^{\infty} M^j = (I - M)^{-1}.
$$

*Proof, step by step.*

1. **Invertibility.** If $I - M$ were singular, there would be $v \neq 0$ with $(I - M)v = 0$, i.e. $Mv = v$. Then $1$ would be an eigenvalue of $M$, which contradicts $|\lambda| < 1$. So $I - M$ is invertible.
2. **Telescoping sum.** Multiply out:
$$
\Big(\sum_{j=0}^{k} M^j\Big)(I - M) = (I + M + \dots + M^k) - (M + M^2 + \dots + M^{k+1}) = I - M^{k+1}.
$$
   Every middle term appears once with $+$ and once with $-$, so they cancel.
3. **Take the limit.** Since $M$ is asymptotically stable, $M^{k+1} \to 0$, so $\big(\sum_{j=0}^{\infty} M^j\big)(I - M) = I$.
4. Multiply both sides on the right by $(I - M)^{-1}$. $\blacksquare$

This is the matrix version of $1 + q + q^2 + \dots = \frac{1}{1-q}$ for $|q| < 1$.

**Proposition 1.1 (lecture).** If $M$ has all eigenvalues strictly inside the unit circle and the input is constant, $u_k = u$, then $x_{k+1} = Mx_k + Bu$ converges, for any $x_0$, to

$$
x_\infty = (I - M)^{-1} B u .
$$

*Proof.* By the solution formula with constant input, $x_k = M^k x_0 + \sum_{j=0}^{k-1} M^{k-1-j} B u = M^k x_0 + \big(\sum_{i=0}^{k-1} M^i\big) B u$, after re-indexing with $i = k-1-j$. As $k\to\infty$: $M^k x_0 \to 0$ and the sum tends to $(I - M)^{-1}$. $\blacksquare$

### 4.3 Equilibria (fixed points)

A state $\bar x$ is an **equilibrium** of $x_{k+1} = f(x_k)$ if $\bar x = f(\bar x)$: once there, you stay there.

**How to find it:** set $x_{k+1} = x_k = \bar x$ and solve. For $x_{k+1} = Mx_k + c$ with constant $c$:

$$
\bar x = M\bar x + c \iff (I - M)\bar x = c \iff \bar x = (I - M)^{-1} c \quad(\text{if } I - M \text{ is invertible}).
$$

**Important:** *if* a trajectory converges to some limit $x_\infty$, then $x_\infty$ must be an equilibrium. Take the limit on both sides of $x_{k+1} = Mx_k + c$ to get $x_\infty = Mx_\infty + c$.

### 4.4 The "shift to the equilibrium" trick

To prove that $x_{k+1} = Mx_k + c$ actually converges to $\bar x$, define the deviation $d_k := x_k - \bar x$. Then

$$
d_{k+1} = x_{k+1} - \bar x = (Mx_k + c) - (M\bar x + c) = M(x_k - \bar x) = M d_k .
$$

The constant $c$ cancels. So $d_k = M^k d_0 \to 0$ when $M$ is asymptotically stable, which means $x_k \to \bar x$. This trick is used in Tasks 1, 2 and 3.

---

## 5. Reachability, controllability, stabilizability

Question: **can the input $u$ steer the state where we want?**

### 5.1 Deriving the controllability matrix $W_c$

From the solution formula with $k = n$:

$$
x_n = A^n x_0 + A^{n-1}Bu_0 + A^{n-2}Bu_1 + \dots + ABu_{n-2} + Bu_{n-1}.
$$

Collect the coefficient matrices side by side and the inputs into one long vector:

$$
x_n - A^n x_0 =
\underbrace{\begin{bmatrix} B & AB & A^2B & \cdots & A^{n-1}B \end{bmatrix}}_{W_c \;\in\; \mathbb{R}^{n\times nm}}
\underbrace{\begin{bmatrix} u_{n-1} \\ u_{n-2} \\ \vdots \\ u_0 \end{bmatrix}}_{\vec u}.
$$

Note the **reversed order**: the *last* input $u_{n-1}$ multiplies $B$, and the *first* input $u_0$ multiplies $A^{n-1}B$, because it has had the most steps to be multiplied by $A$.

This is a linear system of equations $W_c \vec u = b$ with $b = x_n - A^n x_0$. It has a solution for **every** right-hand side $b$ exactly when the columns of $W_c$ span all of $\mathbb{R}^n$, i.e. when $\operatorname{rank}(W_c) = n$.

### 5.2 Definitions (as in the lecture)

> **Reachability (fully state controllable), Def. 1.1:** the system is reachable $\iff \operatorname{rank}(W_c) = n$. Then any $x_0$ can be steered to any $x^\ast$ in finite time.

> **Controllability, Def. 1.2:** the origin $x = 0$ can be reached from any initial state $x_0$ in finite time.

> **Stabilizability, Def. 1.7:** there exists $K$ such that $A - BK$ is asymptotically stable.

$$
\text{Reachable} \;\Rightarrow\; \text{Controllable} \;\Rightarrow\; \text{Stabilizable},
\qquad\text{(the reverse implications are false in general)}.
$$

Controllable but not reachable is possible when the "uncontrollable" part dies out on its own in finite time, i.e. it has eigenvalue $0$. Example: $A$ nilpotent ($A^k = 0$ for some $k$) and $B = 0$.

### 5.3 Why $n$ steps are always enough (Cayley–Hamilton)

From Section 2.8, $A^n B$ is a combination of $B, AB, \dots, A^{n-1}B$. The same holds for $A^{n+1}B$, etc. So **adding more steps never adds new directions**: the reachable set after $n$ steps is the reachable set after any number of steps.

This also explains a pattern you'll see in practice: once $A^{j}B$ depends on the previous columns, all later columns do too.

### 5.4 The reachable subspace

Starting from $x_0 = 0$, the states reachable in $N$ steps are

$$
x_N = \sum_{j=0}^{N-1} A^{N-1-j} B u_j \in \operatorname{span}\{B, AB, \dots, A^{N-1}B\}.
$$

**Reachable subspace** $= \operatorname{im}(W_c)$, the column space of $W_c$. Its dimension is $\operatorname{rank}(W_c)$.

- A target $x_f$ is reachable from $0$ $\iff$ $x_f \in \operatorname{im}(W_c)$ $\iff$ $\operatorname{rank}[W_c \;\; x_f] = \operatorname{rank}(W_c)$.
- **Minimum number of steps $N$:** the smallest $N$ with $x_f \in \operatorname{span}\{B, AB, \dots, A^{N-1}B\}$. Try $N = 1, 2, \dots$ in turn.
- **Input sequence:** solve $x_f = A^{N-1}Bu_0 + \dots + ABu_{N-2} + Bu_{N-1}$ for $u_0,\dots,u_{N-1}$.

### 5.5 Conserved quantities: an intuitive picture of an unreachable direction

If $\operatorname{rank}(W_c) < n$, there is a nonzero vector $w$ orthogonal to the reachable subspace, $w^\top W_c = 0$. In particular $w^\top B = 0$. Suppose additionally $w^\top A = \lambda\, w^\top$ ($w$ is a **left eigenvector**). Then the scalar $z_k := w^\top x_k$ evolves as

$$
z_{k+1} = w^\top A x_k + w^\top B u_k = \lambda\, w^\top x_k + 0 = \lambda z_k .
$$

**The input has no influence on $z_k$ at all.** If $\lambda = 1$, then $z_k$ is constant forever, a **conserved quantity**. That is exactly what makes a direction uncontrollable.

### 5.6 PBH (Popov–Belevitch–Hautus) eigenvalue test

$$
\text{Eigenvalue } \lambda \text{ of } A \text{ is uncontrollable} \iff \operatorname{rank}\begin{bmatrix} \lambda I - A & B \end{bmatrix} < n .
$$

*Why:* the rank drops exactly when some $w \neq 0$ satisfies $w^\top(\lambda I - A) = 0$ and $w^\top B = 0$. That is precisely the situation of Section 5.5: a left eigenvector that the input cannot touch.

You only need to test the eigenvalues of $A$. For any other $\lambda$, $\lambda I - A$ alone already has rank $n$.

**Reading the PBH result:**

| Uncontrollable eigenvalues | Conclusion |
|---|---|
| none | reachable |
| only $\lambda = 0$ | controllable (that part vanishes by itself in finite time), not reachable |
| only $\lvert\lambda\rvert < 1$ | stabilizable (decays by itself asymptotically), but not controllable if $\lambda \neq 0$ |
| some $\lvert\lambda\rvert \ge 1$ | not even stabilizable |

### 5.7 Pole placement (preview)

If the system is reachable, $K$ can be chosen so that $A - BK$ has **any** desired eigenvalues. If it is only stabilizable, the uncontrollable eigenvalues stay where they are whatever $K$ is, but the rest can be moved inside the unit circle.

---

## 6. Observability and detectability

Question: **can the internal state be determined from the measured outputs?**

### 6.1 Deriving the observability matrix $W_o$

Without loss of generality set $u_k = 0$. The known input contribution can always be subtracted. Then $x_k = A^k x_0$ and

$$
\begin{aligned}
y_0 &= Cx_0, \\
y_1 &= Cx_1 = CAx_0, \\
&\;\;\vdots \\
y_{n-1} &= CA^{n-1}x_0 .
\end{aligned}
\qquad\Longrightarrow\qquad
\begin{bmatrix} y_0 \\ y_1 \\ \vdots \\ y_{n-1} \end{bmatrix}
=
\underbrace{\begin{bmatrix} C \\ CA \\ \vdots \\ CA^{n-1} \end{bmatrix}}_{W_o \;\in\; \mathbb{R}^{np\times n}} x_0 .
$$

Recovering $x_0$ means solving $W_o x_0 = \vec y$. The solution is unique exactly when $\ker(W_o) = \{0\}$, i.e. $\operatorname{rank}(W_o) = n$. As before, Cayley–Hamilton says that more than $n$ measurements add no new information.

### 6.2 Definitions

> **Unobservable state, Def. 1.3:** $z \neq 0$ is unobservable if starting at $x_0 = z$ with zero input produces $y_k = 0$ for all $k$.

> **Observability, Def. 1.4/1.5:** the system is observable $\iff \operatorname{rank}(W_o) = n$.

> **Detectability, Def. 1.6:** every unobservable state decays to the origin by itself.

### 6.3 Unobservable and observable subspaces

- **Unobservable subspace** $= \ker(W_o)$: all $x_0$ with $W_o x_0 = 0$. These produce identically zero output, so they **cannot be distinguished from $x_0 = 0$**.
- **Observable subspace** $= \ker(W_o)^\perp$ = the row space of $W_o$.

**Why this matters.** Every $x_0$ splits uniquely as $x_0 = x_{\text{obs}} + x_{\text{unobs}}$, with $x_{\text{obs}}$ in the observable subspace and $x_{\text{unobs}}$ in $\ker(W_o)$. By linearity $\vec y = W_o x_{\text{obs}} + 0$, so the outputs determine only $x_{\text{obs}}$. If $x_0$ itself lies in the observable subspace ($x_{\text{unobs}} = 0$), it is the unique state in that subspace that explains the data. If it lies in the unobservable subspace, the output is zero and nothing can be learned.

**Test:** $x_0$ lies in the observable subspace $\iff$ $x_0$ is orthogonal to every basis vector of $\ker(W_o)$.

### 6.4 Duality and PBH test for observability

Observability is the mirror image of reachability:

$$
(A, C) \text{ observable} \iff (A^\top, C^\top) \text{ reachable},\qquad W_o(A,C) = W_c(A^\top, C^\top)^\top .
$$

**PBH test:**

$$
\text{Eigenvalue } \lambda \text{ is unobservable} \iff \operatorname{rank}\begin{bmatrix} \lambda I - A \\ C \end{bmatrix} < n
\iff \exists\, v\neq 0:\; Av = \lambda v,\; Cv = 0 .
$$

In words: an eigenvector that the sensor does not see. The system is **detectable** iff all unobservable eigenvalues satisfy $|\lambda| < 1$.

---

## 7. State feedback, equilibria and tracking

### 7.1 State feedback and the closed loop

A **linear state-feedback controller** chooses the input from the measured state:

$$
u_k = -K x_k, \qquad K \in \mathbb{R}^{m \times n}.
$$

Substituting it into the system:

$$
x_{k+1} = Ax_k + B(-Kx_k) = (A - BK)\,x_k .
$$

$A - BK$ is the **closed-loop system matrix**. Choosing $K$ changes its eigenvalues. $x_k \to 0$ for all $x_0$ $\iff$ $A - BK$ has all eigenvalues strictly inside the unit circle.

**Open loop vs closed loop:** an open-loop policy is a precomputed input sequence $u_0, u_1,\dots$ that ignores the actual state. A closed-loop policy $u_k = \mu(x_k)$ reacts to the measured state, so it can correct disturbances and model errors.

### 7.2 Tracking control

**Goal:** make $x_k \to x_{\text{ref}}$ for a given constant reference. The natural first attempt is a **MIMO P-controller**:

$$
u_k = -K(x_k - x_{\text{ref}}).
$$

Closed loop:

$$
x_{k+1} = Ax_k - BKx_k + BKx_{\text{ref}} = (A - BK)x_k + BK\,x_{\text{ref}} .
$$

This has the form $x_{k+1} = Mx_k + c$ with $M = A - BK$ and constant $c = BKx_{\text{ref}}$. Sections 4.3 and 4.4 apply directly.

- **Tracking error:** $e_k = x_k - x_{\text{ref}}$.
- **Asymptotic (steady-state) tracking error:** $e_\infty = x_\infty - x_{\text{ref}}$.

### 7.3 Why a P-controller generally leaves a steady-state error

At steady state the state must be an equilibrium. To **hold** the system at some state $\bar x$, the input $\bar u$ must satisfy

$$
\bar x = A\bar x + B\bar u \iff (I - A)\bar x = B\bar u .
$$

So **holding $x_{\text{ref}}$ requires an input $\bar u$ with $B\bar u = (I - A)x_{\text{ref}}$**, which is nonzero unless $Ax_{\text{ref}} = x_{\text{ref}}$. But at $x = x_{\text{ref}}$ the P-controller outputs $u = -K\cdot 0 = 0$. The controller can only produce a nonzero input when there is a nonzero error. The system therefore settles where the remaining error is just large enough to produce the input needed to hold it there.

**Everyday analogy:** a heater controlled by "power $\propto$ (target $-$ room temperature)" switches off exactly at the target. The room keeps losing heat, so it settles a little below the target.

The two standard fixes are the subject of Tasks 2 and 3:

1. **Reference prefilter $F$:** aim at a modified target $Fx_{\text{ref}}$, chosen so that the closed loop settles exactly at $x_{\text{ref}}$.
2. **Integral action:** add an integrator state $e^i_{k+1} = e^i_k + T_s e_k$ that keeps accumulating the error. It can only stop changing when $e_k = 0$. This is the discrete-time **MIMO PI controller**.

### 7.4 Which states can be held at all? (holdable equilibria)

Before asking whether a *particular* controller reaches $x_{\text{ref}}$, ask whether *any* controller can. If $x_k \to x_{\text{ref}}$ and $u_k \to \bar u$, then taking the limit of $x_{k+1} = Ax_k + Bu_k$ gives

$$
x_{\text{ref}} = Ax_{\text{ref}} + B\bar u
\quad\Longleftrightarrow\quad
\boxed{(I - A)\,x_{\text{ref}} = B\,\bar u}
$$

So $x_{\text{ref}}$ can be held **if and only if $(I - A)x_{\text{ref}} \in \operatorname{im}(B)$**. Test: $\operatorname{rank}[\,B \;\; (I-A)x_{\text{ref}}\,] = \operatorname{rank}B$.

Consequences:

- **$B$ square and invertible ($m = n$):** every state can be held, with $\bar u = B^{-1}(I - A)x_{\text{ref}}$.
- **Fewer inputs than states ($m < n$) and $I - A$ invertible:** the holdable states are exactly $\bar x = (I - A)^{-1}B\bar u$. They form a subspace of dimension at most $m$, so with one input it is a single line. Any reference off that line can **never** be held, by any controller: a P-controller, a prefilter and an integrator all fail. This is a property of the plant, not of the controller.
- If $Ax_{\text{ref}} = x_{\text{ref}}$, then $\bar u = 0$ works. That is exactly the Task 1 condition.

### 7.5 Reference prefilter (Task 2)

**Control law:** $u_k = -K(x_k - Fx_{\text{ref}})$. The controller is told to aim at the *modified* target $\tilde x_{\text{ref}} := Fx_{\text{ref}}$. Everything from Task 1 still applies with $\tilde x_{\text{ref}}$ in place of $x_{\text{ref}}$:

$$
x_\infty = (I - (A - BK))^{-1}BK\,F x_{\text{ref}} .
$$

**Zero-error condition.** Set $x_\infty = x_{\text{ref}}$ and multiply by $I - (A - BK)$:

$$
BKFx_{\text{ref}} = (I - A + BK)x_{\text{ref}}
\iff (A - BK + BKF)\,x_{\text{ref}} = x_{\text{ref}} .
$$

**1-D formula.** With $a, b, k, f$ scalars: $bkf = 1 - a + bk$, so

$$
f = \frac{1 - a + bk}{bk} = 1 + \frac{1 - a}{bk}\qquad(b \neq 0,\ k \neq 0).
$$

**Matrix formula** (for every $x_{\text{ref}}$ at once). We need $A - BK + BKF = I$, i.e. $BK(F - I) = I - A$. If $BK$ is invertible:

$$
\boxed{F = (BK)^{-1}(I - A) + I}
$$

$BK$ is $n\times n$ but has rank at most $\operatorname{rank}B \le m$. So invertibility needs $m \ge n$ and $\operatorname{rank}B = n$.

**Interpretation: a prefilter is feedforward.** Rewrite the law as

$$
u_k = -K(x_k - x_{\text{ref}}) + \underbrace{K(F - I)x_{\text{ref}}}_{u_{\text{ff}}} .
$$

With the formula for $F$: $Bu_{\text{ff}} = BK(F - I)x_{\text{ref}} = (I - A)x_{\text{ref}}$. That is exactly the **holding input** $\bar u$ of §7.4. The prefilter adds the input needed to hold the target, computed in advance from the model. At the target the feedback part is zero, and the feedforward part does the holding.

**Limitations.**

1. The feedforward is computed from $A$ and $B$. If the real plant differs, the computed holding input is wrong and an error comes back. Nothing measures or corrects this.
2. The general formula needs $BK$ invertible ($m \ge n$). With $m < n$, a prefilter can only fix a reference that is holdable (§7.4).

### 7.6 Block matrices

Task 3 stacks two vectors into one bigger state, so we need matrices made of blocks.

**Block multiplication** works like ordinary multiplication, with blocks in place of numbers, as long as the sizes fit. Keep the order of the factors, because blocks do not commute:

$$
\begin{bmatrix} P & Q \\ R & S \end{bmatrix}\begin{bmatrix} v \\ w \end{bmatrix}
= \begin{bmatrix} Pv + Qw \\ Rv + Sw \end{bmatrix}.
$$

*Example:* with $n = 2$, $\begin{bmatrix} x \\ e^i \end{bmatrix} \in \mathbb{R}^4$, and $\begin{bmatrix} T_s I & I \end{bmatrix}\begin{bmatrix} x \\ e^i \end{bmatrix} = T_s x + e^i$.

**Determinant of a block matrix with scalar-multiple-of-identity blocks.** If all four blocks are multiples of $I_n$:

$$
\det\begin{bmatrix} pI & qI \\ rI & sI \end{bmatrix} = (ps - qr)^n .
$$

*Why:* reorder the rows and columns so that the coordinates belonging together ($x_1$ with $e^i_1$, $x_2$ with $e^i_2$, …) sit next to each other. Reordering rows and columns in the same way does not change the determinant. The matrix becomes block-diagonal with $n$ copies of $\begin{bmatrix} p & q \\ r & s\end{bmatrix}$, and the determinant of a block-diagonal matrix is the product of the blocks' determinants. Consequently, the eigenvalues of such a matrix are those of the $2\times2$ matrix $\begin{bmatrix} p & q \\ r & s \end{bmatrix}$, each repeated $n$ times.

### 7.7 Integral action: the MIMO PI controller (Task 3)

$$
u_k = -K_p e_k - K_I e^i_k,
\qquad e_k = x_k - x_{\text{ref}},
\qquad e^i_{k+1} = e^i_k + T_s e_k,\quad e^i_0 = 0,
$$

with $K_p, K_I \in \mathbb{R}^{m\times n}$. Unrolling gives $e^i_k = T_s\sum_{j=0}^{k-1}e_j$, the running sum of all past errors. It is the discrete version of $\int e\,dt$.

**Why it removes the steady-state error.** If everything converges, the integrator equation in the limit reads $e^i_\infty = e^i_\infty + T_s e_\infty$. Hence $T_s e_\infty = 0$, and since $T_s > 0$, $e_\infty = 0$. The integrator can only stop changing when the error is exactly zero. Meanwhile $-K_I e^i_\infty$ supplies the holding input $\bar u$ of §7.4: the integrator *learns* it from the data instead of computing it from the model. That is why model errors don't matter.

**Augmented state-space model.** Substitute $u_k$ into the plant and stack $z_k = \begin{bmatrix} x_k \\ e^i_k\end{bmatrix}$:

$$
z_{k+1} =
\underbrace{\begin{bmatrix} A - BK_p & -BK_I \\ T_s I & I \end{bmatrix}}_{\tilde A} z_k
+ \underbrace{\begin{bmatrix} BK_p \\ -T_s I \end{bmatrix}}_{\tilde B} x_{\text{ref}} .
$$

This is the same kind of system as in Task 1, $z_{k+1} = \tilde A z_k + \tilde B x_{\text{ref}}$ with a constant "input" $x_{\text{ref}}$. The derivation is written out in full in the Solutions file (Part E).

**Proposition 1.3 (lecture).** If $\tilde A$ is asymptotically stable, then $x_\infty = x_{\text{ref}}$ for every constant $x_{\text{ref}}$.

**When can $\tilde A$ be made asymptotically stable?** Write $\tilde A = A_{\text{aug}} - B_{\text{aug}}[\,K_p\;K_I\,]$ with

$$
A_{\text{aug}} = \begin{bmatrix} A & 0 \\ T_s I & I \end{bmatrix},\qquad
B_{\text{aug}} = \begin{bmatrix} B \\ 0 \end{bmatrix}.
$$

This is ordinary state feedback on the augmented system, so we need $(A_{\text{aug}}, B_{\text{aug}})$ to be stabilizable. $A_{\text{aug}}$ is block-triangular, so its eigenvalues are those of $A$ plus $n$ eigenvalues equal to $1$ (from the integrators). PBH at $\lambda = 1$:

$$
\begin{bmatrix} I - A_{\text{aug}} & B_{\text{aug}} \end{bmatrix}
= \begin{bmatrix} I - A & 0 & B \\ -T_s I & 0 & 0 \end{bmatrix}.
$$

The middle block column is zero. The bottom block row $[-T_s I,\ 0,\ 0]$ already has rank $n$, and it uses up the first block column. The top block row then needs rank $n$ from the $B$ columns alone. So the rank is $n + \operatorname{rank}B$, and

$$
\lambda = 1 \text{ can be moved} \iff \operatorname{rank}B = n \;(\Rightarrow m \ge n).
$$

**With fewer inputs than states ($m < n$), $\lambda = 1$ is always an eigenvalue of $\tilde A$**, whatever $K_p, K_I$ are. $K_I$ is $m\times n$ with $m < n$, so there is a $w \neq 0$ with $K_I w = 0$. Then

$$
\tilde A\begin{bmatrix} 0 \\ w \end{bmatrix} = \begin{bmatrix} -BK_I w \\ w \end{bmatrix} = \begin{bmatrix} 0 \\ w \end{bmatrix}.
$$

That integrator direction has no effect on the input, so it just keeps summing up its error component forever.

**Practical fix when $m < n$:** integrate only $m$ chosen output errors, e.g. $e^i_{k+1} = e^i_k + T_s(Cx_k - r)$ with a $C$ of $m$ rows. That forces those outputs to their references. The remaining states settle wherever the holdable equilibrium puts them (§7.4).

**Merely Lyapunov stable $\tilde A$** (some $|\lambda| = 1$, none outside):

- $\lambda = e^{\pm i\varphi} \neq 1$: $\tilde A^k$ does not decay, so the transient never dies. The error oscillates forever.
- $\lambda = 1$: $I - \tilde A$ is singular, so there may be no equilibrium or infinitely many. The geometric-series argument of Task 1 breaks down. Typically one integrator component grows without bound while an error remains. This is **integrator windup**.

### 7.8 Choosing the gains (pole placement in practice)

The task only requires stability, so any stabilizing gains will do. Three hand-friendly methods are used in the Solutions:

1. **Prescribe the closed-loop matrix (when $B$ is invertible).** Choose a stable target matrix $D$ and set $A - BK = D$, i.e. $K = B^{-1}(A - D)$. With $D = 0.5I$, both eigenvalues are $0.5$.
2. **Coefficient matching (single input, $n = 2$).** With $K = [k_1, k_2]$:
   $$
   \operatorname{tr}(A - BK) = \operatorname{tr}A - b_1k_1 - b_2k_2,
   $$
   $$
   \det(A - BK) = \det A - k_1(b_1a_{22} - b_2a_{12}) - k_2(b_2a_{11} - b_1a_{21}).
   $$
   (The $k_1k_2$ terms cancel, so both are *linear* in $k$.) For desired eigenvalues $\lambda_1, \lambda_2$, set trace $= \lambda_1 + \lambda_2$ and determinant $= \lambda_1\lambda_2$, then solve two linear equations. `scipy.signal.place_poles` does the same automatically.
3. **PI with scalar blocks (when $B$ is invertible).** Choose $A - BK_p = \alpha I$ and $BK_I = \beta I$, i.e. $K_p = B^{-1}(A - \alpha I)$ and $K_I = \beta B^{-1}$. Then $\tilde A = \begin{bmatrix} \alpha I & -\beta I \\ T_s I & I\end{bmatrix}$, and by §7.6 its eigenvalues solve
   $$
   (\lambda - \alpha)(\lambda - 1) + \beta T_s = 0
   \iff \lambda^2 - (1 + \alpha)\lambda + (\alpha + \beta T_s) = 0 .
   $$
   For desired $\lambda_1, \lambda_2$: $\alpha = \lambda_1 + \lambda_2 - 1$ and $\beta = (\lambda_1\lambda_2 - \alpha)/T_s$.

**About $T_s$.** It only appears as the product $K_I T_s$, so a different sampling time just rescales $K_I$. Steady-state results never depend on it.

### 7.9 The three tracking controllers compared

| | P (Task 1) | P + prefilter (Task 2) | PI (Task 3) |
|---|---|---|---|
| law | $-K(x - x_{\text{ref}})$ | $-K(x - Fx_{\text{ref}})$ | $-K_pe - K_Ie^i$ |
| zero error if | $Ax_{\text{ref}} = x_{\text{ref}}$ | $(A - BK + BKF)x_{\text{ref}} = x_{\text{ref}}$ | $\tilde A$ asymptotically stable |
| needs | — | exact model; $BK$ invertible for all $x_{\text{ref}}$ ($m \ge n$) | $\operatorname{rank}B = n$ for a full-state integrator |
| holding input comes from | a remaining error | the model (feedforward) | the integrator (learned) |
| model error | error | error comes back | still zero |
| extra states | 0 | 0 | $n$ |
| typical drawback | steady-state error | sensitive to model mismatch | slower, overshoot, windup risk |

---

## 8. Cheat sheet

| Concept | Formula / test |
|---|---|
| Characteristic equation | $\det(\lambda I - A) = 0$ |
| $2\times2$ eigenvalues | $\lambda^2 - \operatorname{tr}(A)\lambda + \det(A) = 0$ |
| Eigenvalue checks | $\sum\lambda_i = \operatorname{tr}A$, $\prod\lambda_i = \det A$ |
| Solution | $x_k = A^kx_0 + \sum_{j=0}^{k-1}A^{k-1-j}Bu_j$ |
| Asymptotic stability (DT) | all $\lvert\lambda_i\rvert < 1$ |
| Geometric series | $\sum_{j\ge0}M^j = (I-M)^{-1}$ if $M$ stable |
| Constant input limit | $x_\infty = (I - A)^{-1}Bu$ |
| Equilibrium of $x^+ = Mx + c$ | $\bar x = (I - M)^{-1}c$ |
| $I - M$ singular | $\iff 1$ is an eigenvalue of $M$ |
| Controllability matrix | $W_c = [B,\;AB,\;\dots,\;A^{n-1}B]$ |
| Reachable | $\operatorname{rank} W_c = n$ |
| Reachable subspace | $\operatorname{im}(W_c)$ |
| $x_f$ reachable from 0 | $\operatorname{rank}[W_c\;x_f] = \operatorname{rank}W_c$ |
| PBH uncontrollable $\lambda$ | $\operatorname{rank}[\lambda I - A,\;B] < n$ |
| Controllable (DT) | all uncontrollable eigenvalues $= 0$ |
| Stabilizable | all uncontrollable eigenvalues $\lvert\lambda\rvert<1$ |
| Observability matrix | $W_o = [C;\;CA;\;\dots;\;CA^{n-1}]$ |
| Observable | $\operatorname{rank} W_o = n$ |
| Unobservable subspace | $\ker(W_o)$ |
| Observable subspace | $\ker(W_o)^\perp$ = row space of $W_o$ |
| PBH unobservable $\lambda$ | $\operatorname{rank}[\lambda I - A;\;C] < n$ |
| Detectable | all unobservable eigenvalues $\lvert\lambda\rvert<1$ |
| Closed loop (state feedback) | $x_{k+1} = (A - BK)x_k$ |
| P-tracking closed loop | $x_{k+1} = (A-BK)x_k + BKx_{\text{ref}}$ |
| P-tracking limit | $x_\infty = (I - A + BK)^{-1}BKx_{\text{ref}}$ |
| P-tracking error | $x_\infty - x_{\text{ref}} = (I - A + BK)^{-1}(A - I)x_{\text{ref}}$ |
| $x_{\text{ref}}$ holdable at all | $(I - A)x_{\text{ref}} \in \operatorname{im}B$, holding input $B\bar u = (I-A)x_{\text{ref}}$ |
| Prefilter limit | $x_\infty = (I - A + BK)^{-1}BKFx_{\text{ref}}$ |
| Prefilter, zero error | $(A - BK + BKF)x_{\text{ref}} = x_{\text{ref}}$ |
| Prefilter 1-D | $f = (1 - a + bk)/(bk)$ |
| Prefilter matrix | $F = (BK)^{-1}(I - A) + I$ ($BK$ invertible) |
| Integrator | $e^i_{k+1} = e^i_k + T_se_k$ |
| Augmented PI model | $\tilde A = \begin{bmatrix} A - BK_p & -BK_I \\ T_sI & I\end{bmatrix}$, $\tilde B = \begin{bmatrix} BK_p \\ -T_sI\end{bmatrix}$ |
| PI zero error | $\tilde A$ asymptotically stable $\Rightarrow x_\infty = x_{\text{ref}}$ |
| Full-state PI possible | $\operatorname{rank}B = n$; if $m < n$, $\lambda = 1 \in \operatorname{eig}(\tilde A)$ always |
| Block determinant | $\det\begin{bmatrix} pI & qI \\ rI & sI\end{bmatrix} = (ps - qr)^n$ |
