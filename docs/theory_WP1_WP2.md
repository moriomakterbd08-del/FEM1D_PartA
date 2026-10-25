# Mathematical Foundations & Proofs (WP1 & WP2)

## 1. WP1: Background Analysis in Sobolev Spaces

### 1.1 The Classical Boundary Value Problem
Consider the two-point Dirichlet problem:
$$
-u''(x) = f(x), \quad x \in (0, 1), \qquad u(0) = u(1) = 0 \tag{1}
$$
If $f \in C([0, 1])$, a classical solution is $u \in C^2([0, 1])$ satisfying (1) pointwise. However, for general source terms $f \in L^2(0, 1)$, we seek solutions in the generalized Sobolev space $H_0^1(0, 1)$.

### 1.2 Weak Derivatives and the Space $H_0^1(0, 1)$
- **$L^2(0, 1)$:** The Hilbert space of square-integrable functions with inner product $(u, v)_{L^2} = \int_0^1 u(x) v(x) \, dx$ and norm $\|u\|_{L^2} = \sqrt{(u, u)_{L^2}}$.
- **Weak Derivative:** A function $u \in L^2(0, 1)$ possesses a weak derivative $w \in L^2(0, 1)$ (denoted $u' = w$) if:
  $$
  \int_0^1 u(x) \phi'(x) \, dx = -\int_0^1 w(x) \phi(x) \, dx \quad \forall \phi \in C_0^\infty(0, 1)
  $$
- **Sobolev Space $H^1(0, 1)$:**
  $$
  H^1(0, 1) = \left\{ v \in L^2(0, 1) : v' \in L^2(0, 1) \right\}
  $$
  equipped with the norm $\|v\|_{H^1} = \sqrt{\|v\|_{L^2}^2 + \|v'\|_{L^2}^2}$.
- **Space $H_0^1(0, 1)$:** The closure of $C_0^\infty(0, 1)$ with respect to the $H^1$ norm, consisting of functions in $H^1(0, 1)$ that vanish on the boundary: $v(0) = v(1) = 0$.

---

### 1.3 Elementary Proof of the 1D Poincaré Inequality
**Theorem (1D Poincaré Inequality):**
For any $v \in H_0^1(0, 1)$,
$$
\|v\|_{L^2(0, 1)} \le \frac{1}{\pi} \|v'\|_{L^2(0, 1)} \tag{2}
$$

**Elementary Proof (using Fundamental Theorem of Calculus):**
Since $C_0^\infty(0, 1)$ is dense in $H_0^1(0, 1)$, let $v \in C_0^\infty(0, 1)$. Because $v(0) = 0$:
$$
v(x) = \int_0^x v'(t) \, dt
$$
Applying the Cauchy–Schwarz inequality:
$$
|v(x)|^2 = \left| \int_0^x 1 \cdot v'(t) \, dt \right|^2 \le \left( \int_0^x 1^2 \, dt \right) \left( \int_0^x |v'(t)|^2 \, dt \right) \le x \int_0^1 |v'(t)|^2 \, dt
$$
Integrating both sides over $x \in (0, 1)$:
$$
\int_0^1 |v(x)|^2 \, dx \le \left( \int_0^1 x \, dx \right) \left( \int_0^1 |v'(t)|^2 \, dt \right) = \frac{1}{2} \int_0^1 |v'(t)|^2 \, dt
$$
Taking square roots gives $\|v\|_{L^2} \le \frac{1}{\sqrt{2}} \|v'\|_{L^2}$.

To obtain the sharp constant $C_P = \frac{1}{\pi}$, expand $v$ in an orthonormal Fourier sine basis $\{ \sqrt{2} \sin(k\pi x) \}_{k=1}^\infty$:
$$
v(x) = \sum_{k=1}^\infty c_k \sqrt{2} \sin(k\pi x), \quad v'(x) = \sum_{k=1}^\infty c_k k\pi \sqrt{2} \cos(k\pi x)
$$
By Parseval's identity:
$$
\|v\|_{L^2}^2 = \sum_{k=1}^\infty c_k^2 \le \frac{1}{\pi^2} \sum_{k=1}^\infty c_k^2 (k\pi)^2 = \frac{1}{\pi^2} \|v'\|_{L^2}^2
$$
Taking the square root yields the sharp Poincaré inequality: $\|v\|_{L^2} \le \frac{1}{\pi} \|v'\|_{L^2}$. $\quad \blacksquare$

---

### 1.4 The Weak Formulation and Lax–Milgram Lemma
Multiply (1) by an arbitrary test function $v \in H_0^1(0, 1)$ and integrate over $(0, 1)$:
$$
-\int_0^1 u''(x) v(x) \, dx = \int_0^1 f(x) v(x) \, dx
$$
Integrating by parts on the left-hand side:
$$
-[u'(x)v(x)]_0^1 + \int_0^1 u'(x) v'(x) \, dx = \int_0^1 f(x) v(x) \, dx
$$
Since $v(0) = v(1) = 0$, the boundary term vanishes identically.

**Weak Problem:** Find $u \in H_0^1(0, 1)$ such that:
$$
a(u, v) = \ell(v) \quad \forall v \in H_0^1(0, 1) \tag{3}
$$
where:
$$
a(u, v) = \int_0^1 u'(x) v'(x) \, dx, \qquad \ell(v) = \int_0^1 f(x) v(x) \, dx
$$

**Theorem (Lax–Milgram Lemma):**
Let $V$ be a real Hilbert space with norm $\|\cdot\|_V$. Let $a: V \times V \to \mathbb{R}$ be a bilinear form and $\ell: V \to \mathbb{R}$ a linear functional satisfying:
1. **Continuity of $a$:** $|a(u, v)| \le M \|u\|_V \|v\|_V$ for all $u, v \in V$.
2. **Coercivity (V-ellipticity) of $a$:** $a(v, v) \ge \alpha \|v\|_V^2$ for some $\alpha > 0$ for all $v \in V$.
3. **Continuity of $\ell$:** $|\ell(v)| \le \|\ell\|_{V^*} \|v\|_V$ for all $v \in V$.

Then there exists a unique $u \in V$ such that $a(u, v) = \ell(v)$ for all $v \in V$, and $\|u\|_V \le \frac{1}{\alpha} \|\ell\|_{V^*}$.

**Application to 1D Poisson:**
Let $V = H_0^1(0, 1)$ with $\|v\|_V = \|v\|_{H^1} = \sqrt{\|v\|_{L^2}^2 + \|v'\|_{L^2}^2}$.
1. **Continuity:**
   $|a(u, v)| = \left| \int_0^1 u' v' dx \right| \le \|u'\|_{L^2} \|v'\|_{L^2} \le \|u\|_{H^1} \|v\|_{H^1} \implies M = 1$.
2. **Coercivity:**
   Using the Poincaré inequality $\|v\|_{L^2}^2 \le \frac{1}{\pi^2} \|v'\|_{L^2}^2$:
   $$
   \|v\|_{H^1}^2 = \|v\|_{L^2}^2 + \|v'\|_{L^2}^2 \le \left( 1 + \frac{1}{\pi^2} \right) \|v'\|_{L^2}^2 = \left( 1 + \frac{1}{\pi^2} \right) a(v, v)
   $$
   Therefore:
   $$
   a(v, v) \ge \alpha \|v\|_{H^1}^2 \quad \text{with } \alpha = \frac{1}{1 + \pi^{-2}} > 0
   $$
3. **Continuity of $\ell$:**
   $|\ell(v)| = \left| \int_0^1 f v dx \right| \le \|f\|_{L^2} \|v\|_{L^2} \le \|f\|_{L^2} \|v\|_{H^1}$.

By Lax–Milgram, the weak problem has a **unique solution** $u \in H_0^1(0, 1)$. $\quad \blacksquare$

---

## 2. WP2: FEM Theory for the Poisson Problem

### 2.1 The Galerkin Method & Galerkin Orthogonality
Let $V_h \subset V = H_0^1(0, 1)$ be a finite-dimensional subspace of continuous piecewise linear functions vanishing at $x = 0$ and $x = 1$.
The discrete Galerkin problem is:
$$
\text{Find } u_h \in V_h \text{ such that } a(u_h, v_h) = \ell(v_h) \quad \forall v_h \in V_h \tag{4}
$$
Subtracting (4) from the continuous equation (3) with $v = v_h \in V_h$:
$$
a(u - u_h, v_h) = 0 \quad \forall v_h \in V_h \tag{5}
$$
Equation (5) is the fundamental property of **Galerkin Orthogonality**: the discretization error $u - u_h$ is orthogonal to the discrete trial space $V_h$ with respect to the energy inner product $a(\cdot, \cdot)$.

---

### 2.2 Full Proof of Céa's Lemma
**Theorem (Céa's Lemma):**
Let $u \in V$ and $u_h \in V_h$ be the solutions of the continuous and Galerkin problems, respectively. Then:
$$
\|u - u_h\|_V \le \frac{M}{\alpha} \inf_{v_h \in V_h} \|u - v_h\|_V \tag{6}
$$

**Proof:**
For any arbitrary $v_h \in V_h$, write $u - u_h = (u - v_h) + (v_h - u_h)$.
Notice that $w_h \equiv u_h - v_h \in V_h$. Using coercivity:
$$
\alpha \|u - u_h\|_V^2 \le a(u - u_h, u - u_h)
$$
By Galerkin orthogonality, $a(u - u_h, w_h) = 0$. Therefore, we can add $a(u - u_h, v_h - u_h) = 0$:
$$
a(u - u_h, u - u_h) = a(u - u_h, u - u_h) + a(u - u_h, v_h - u_h) = a(u - u_h, u - v_h)
$$
Applying the continuity of $a(\cdot, \cdot)$:
$$
a(u - u_h, u - v_h) \le M \|u - u_h\|_V \|u - v_h\|_V
$$
Combining the inequalities:
$$
\alpha \|u - u_h\|_V^2 \le M \|u - u_h\|_V \|u - v_h\|_V
$$
Dividing by $\alpha \|u - u_h\|_V$ (if $u = u_h$, the lemma holds trivially):
$$
\|u - u_h\|_V \le \frac{M}{\alpha} \|u - v_h\|_V
$$
Since this holds for every $v_h \in V_h$, taking the infimum over all $v_h \in V_h$ yields:
$$
\|u - u_h\|_V \le \frac{M}{\alpha} \inf_{v_h \in V_h} \|u - v_h\|_V
$$
**Special Case for 1D Poisson:** Under the energy semi-norm $|v|_{H^1} = \sqrt{a(v, v)}$, $M = \alpha = 1$, giving:
$$
|u - u_h|_{H^1} = \inf_{v_h \in V_h} |u - v_h|_{H^1}
$$
meaning the FEM solution $u_h$ is the **best possible approximation** of $u$ in $V_h$ in the energy norm. $\quad \blacksquare$

---

### 2.3 1D Interpolation Error Estimates
Let $I_h u \in V_h$ be the piecewise linear nodal interpolant of $u$:
$$
I_h u(x_j) = u(x_j), \quad j = 0, 1, \dots, N
$$

**Theorem (Interpolation Error in $H^1$ and $L^2$):**
On an element $e = [x_{i-1}, x_i]$ of length $h_i$:
$$
\|(u - I_h u)'\|_{L^2(e)} \le \frac{h_i}{\pi} \|u''\|_{L^2(e)}, \qquad \|u - I_h u\|_{L^2(e)} \le \frac{h_i^2}{\pi^2} \|u''\|_{L^2(e)}
$$

**Proof:**
Set $\eta(x) = u(x) - I_h u(x)$. Since $\eta(x_{i-1}) = \eta(x_i) = 0$, by Rolle's theorem there exists $\xi \in (x_{i-1}, x_i)$ where $\eta'(\xi) = 0$.
For any $x \in e$:
$$
\eta'(x) = \int_\xi^x \eta''(t) \, dt = \int_\xi^x u''(t) \, dt
$$
because $(I_h u)'' = 0$ on the interior of $e$.
Applying Cauchy–Schwarz:
$$
|\eta'(x)|^2 \le |x - \xi| \left| \int_\xi^x |u''(t)|^2 \, dt \right| \le h_i \int_{x_{i-1}}^{x_i} |u''(t)|^2 \, dt
$$
Integrating over $e$:
$$
\int_{x_{i-1}}^{x_i} |\eta'(x)|^2 \, dx \le h_i^2 \int_{x_{i-1}}^{x_i} |u''(t)|^2 \, dt \implies \|\eta'\|_{L^2(e)} \le h_i \|u''\|_{L^2(e)}
$$
(With sharp Fourier analysis, the constant is $\frac{1}{\pi}$).
Furthermore, since $\eta \in H_0^1(e)$, applying Poincaré's inequality on element $e$ gives:
$$
\|\eta\|_{L^2(e)} \le \frac{h_i}{\pi} \|\eta'\|_{L^2(e)} \le \frac{h_i^2}{\pi^2} \|u''\|_{L^2(e)}
$$
Summing over all elements $e = 1, \dots, N$:
$$
|u - I_h u|_{H^1(0, 1)} \le C_1 h \|u''\|_{L^2(0, 1)}, \qquad \|u - I_h u\|_{L^2(0, 1)} \le C_2 h^2 \|u''\|_{L^2(0, 1)} \quad \blacksquare
$$

Combining with Céa's lemma yields the first main theoretical bound:
$$
|u - u_h|_{H^1(0, 1)} \le |u - I_h u|_{H^1(0, 1)} \le C h \|u''\|_{L^2(0, 1)} = \mathcal{O}(h)
$$

---

### 2.4 Aubin–Nitsche Duality Argument ($L^2$ Error Estimate)
Céa's lemma directly controls the energy norm $|u - u_h|_{H^1} = \mathcal{O}(h)$. To prove that the $L^2$ error converges at the optimal higher rate $\mathcal{O}(h^2)$, we invoke the **Aubin–Nitsche duality trick**.

**Theorem (Aubin–Nitsche):**
$$
\|u - u_h\|_{L^2(0, 1)} \le C h |u - u_h|_{H^1(0, 1)} \le C' h^2 \|u''\|_{L^2(0, 1)} = \mathcal{O}(h^2)
$$

**Proof:**
Consider the auxiliary (dual) adjoint problem: find $w \in H_0^1(0, 1)$ such that
$$
-w''(x) = u(x) - u_h(x), \quad x \in (0, 1), \qquad w(0) = w(1) = 0 \tag{7}
$$
In weak form, this means:
$$
a(v, w) = (v, u - u_h)_{L^2} \quad \forall v \in H_0^1(0, 1)
$$
Choosing the test function $v = u - u_h \in H_0^1(0, 1)$:
$$
\|u - u_h\|_{L^2}^2 = (u - u_h, u - u_h)_{L^2} = a(u - u_h, w)
$$
By Galerkin orthogonality, $a(u - u_h, v_h) = 0$ for all $v_h \in V_h$. Choosing $v_h = I_h w \in V_h$ (the nodal interpolant of $w$):
$$
\|u - u_h\|_{L^2}^2 = a(u - u_h, w - I_h w)
$$
Applying Cauchy–Schwarz to $a(\cdot, \cdot)$:
$$
\|u - u_h\|_{L^2}^2 \le |u - u_h|_{H^1} |w - I_h w|_{H^1}
$$
Using the interpolation estimate on $w$:
$$
|w - I_h w|_{H^1} \le C_1 h \|w''\|_{L^2}
$$
From the dual equation (7), $w'' = -(u - u_h)$, so $\|w''\|_{L^2} = \|u - u_h\|_{L^2}$. Substituting this back:
$$
\|u - u_h\|_{L^2}^2 \le |u - u_h|_{H^1} \left( C_1 h \|u - u_h\|_{L^2} \right)
$$
Dividing both sides by $\|u - u_h\|_{L^2}$:
$$
\|u - u_h\|_{L^2} \le C_1 h |u - u_h|_{H^1}
$$
Finally, substituting the Céa bound $|u - u_h|_{H^1} \le C h \|u''\|_{L^2}$:
$$
\|u - u_h\|_{L^2} \le C' h^2 \|u''\|_{L^2} = \mathcal{O}(h^2) \quad \blacksquare
$$

---

## 3. The Nodal Exactness Theorem (1D Green's Function)

**Theorem:**
Let $u$ be the solution of $-u'' = f$, $u(0)=u(1)=0$, and let $u_h \in V_h$ be the piecewise linear Galerkin FEM solution. If the load vector entries $F_j = \int_0^1 f(x) \phi_j(x) dx$ are integrated exactly, then:
$$
u_h(x_j) = u(x_j) \quad \forall \text{ nodes } x_j
$$

**Proof:**
The Green's function $G(x, y)$ for the two-point BVP satisfies $- \frac{\partial^2 G}{\partial y^2}(x, y) = \delta(x - y)$ with $G(x, 0) = G(x, 1) = 0$. Its closed-form expression is:
$$
G(x, y) = \begin{cases} y(1 - x), & 0 \le y \le x \\ x(1 - y), & x \le y \le 1 \end{cases}
$$
For any function $w \in H_0^1(0, 1)$, integration by parts gives:
$$
a(w, G(x, \cdot)) = \int_0^1 w'(y) \frac{\partial G}{\partial y}(x, y) \, dy = w(x)
$$
Setting $w = u - u_h$:
$$
u(x) - u_h(x) = a(u - u_h, G(x, \cdot))
$$
Now let $x = x_j$ be an interior mesh node. Observe that $y \mapsto G(x_j, y)$:
1. Is continuous on $[0, 1]$.
2. Vanishes at $y = 0$ and $y = 1$.
3. Has constant slope $1 - x_j$ on $[0, x_j]$ and constant slope $-x_j$ on $[x_j, 1]$.

Because $x_j$ is a node of the triangulation, the kink in $G(x_j, y)$ occurs exactly at a mesh boundary! Consequently, $G(x_j, \cdot)$ is continuous and piecewise linear on the mesh:
$$
G(x_j, \cdot) \in V_h
$$
By Galerkin orthogonality, $a(u - u_h, v_h) = 0$ for all $v_h \in V_h$. Choosing $v_h = G(x_j, \cdot) \in V_h$:
$$
u(x_j) - u_h(x_j) = a(u - u_h, G(x_j, \cdot)) = 0
$$
Thus, $u_h(x_j) = u(x_j)$ exactly at every mesh node. $\quad \blacksquare$
