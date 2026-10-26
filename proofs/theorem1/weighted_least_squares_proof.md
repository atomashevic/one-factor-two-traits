---
title: 'Online Resource 2 for “Fitting One Factor to Two Traits: Population Solutions and Fit Indices”'
subtitle: 'Global weighted least-squares fit for two exchangeable blocks'
---

## Scope

This Online Resource proves Proposition 1 of the article, the weighted
least-squares result (Section 5.1, with the proof outline in Appendix C). It establishes the
unrestricted global reduction used there for every positive within-block and
cross-block weight, including the zero-correlation endpoint and the degenerate
weight boundary. The argument uses elementary inequalities, compactness,
differentiation, and a quadratic minimization. Numerical optimization does not
enter the proof. The unweighted case is Theorem 2 of the article, proved in
Section 11 of Online Resource 1; Section 6 below relates the two results and
states the limits of the present statement, and Section 7 lists exact algebra
checks that supplement, but do not replace, the argument.

**Notation.** Symbols follow the article where they denote the same objects:
\(\omega=w_C/w_W\) is the weight ratio and \(\phi_\omega\) the threshold of
Proposition 1, and \(\varpi\) is a Karush–Kuhn–Tucker
multiplier. The constants \(C\), \(K\) and \(d\) in (W.1) and the bound
\(E_\omega\) in (W.11) are local to this resource. In the SymPy snippet of
Section 7, which `checks/exact_algebra.py` executes, `k` denotes \(\omega\)
and `mu` denotes \(\varpi\).

## Statement

Let each block contain an integer \(m\ge2\) items. Assume
\[
0<L<1,\qquad 0\le\phi<1,\qquad w_W>0,\qquad w_C>0.
\]
The population correlations equal \(L\) within blocks and \(c=L\phi\)
across blocks. Define
\[
\omega=\frac{w_C}{w_W},\qquad n=m-1,\qquad
C=\frac{m(m-1)}2,\qquad K=m(m-1),\qquad d=\omega m-m+1.
\tag{W.1}
\]
For \(l\in[-1,1]^{2m}\), minimize
\[
H(l)=\frac{G_w(l)}{w_W}
=\sum_{\mathrm{within}}(L-l_i l_j)^2
+\omega\sum_{\mathrm{across}}(c-l_i l_j)^2.
\tag{W.2}
\]
Each unordered within-block pair appears once. The cross-block sum contains
\(m^2\) pairs. The weights are constant within each pair class. They need not
equal one another.

**Theorem (Proposition 1 of the article: weighted least-squares solution).** The bounded problem attains its minimum. For \(0<\phi<1\), every
global minimizer has a common sign and is uniform within each block. Choose
the common sign to make its two block loadings \(\alpha,\beta>0\).

If \(d>0\), define
\[
\phi_\omega=\frac{d}{2\omega m}
=\frac12-\frac{m-1}{2\omega m}.
\tag{W.3}
\]
For \(0<\phi<\phi_\omega\), the two ordered solutions satisfy
\[
\alpha^2+\beta^2=L,\qquad
\alpha\beta=\frac{\omega mL\phi}{d}.
\tag{W.4}
\]
They differ by block interchange. Their minimum value is
\[
H_{\min}=CL^2-\frac{\omega m^2(m-1)}{d}L^2\phi^2.
\tag{W.5}
\]
For \(\phi\ge\phi_\omega\), including threshold equality, the unique
nonnegative solution is symmetric:
\[
\alpha^2=\beta^2=t_s
=L\frac{m-1+\omega m\phi}{m-1+\omega m},\qquad
H_{\min}=\frac{\omega m^2(m-1)}{m-1+\omega m}L^2(1-\phi)^2.
\tag{W.6}
\]
At \(\phi=\phi_\omega\), \(t_s=L/2\), and (W.5) and (W.6) agree.

If \(d\le0\), (W.6) gives the unique nonnegative global minimizer for
every \(\phi>0\). Thus the stated positive threshold occurs exactly when
\(\omega>(m-1)/m\). Multiplying each displayed minimum by \(w_W\) gives
the minimum of \(G_w\).

All positive-correlation minimizers satisfy \(0<l_i^2<L<1\). Uniqueness
above refers to the nonnegative loading vector. A simultaneous sign change
preserves the fitted covariance. Below the positive threshold, block
interchange gives two distinct fitted covariances. At and above it, the
fitted covariance is unique.

At \(\phi=0\), the complete classification is:

1. If \(d>0\), one block is fitted exactly and the other has zero loadings.
   For \(m\ge3\), the retained loadings equal \(\sqrt L\) up to their common
   sign. For \(m=2\), their magnitudes are
   \((r,L/r)\), \(L\le r\le1\). In all cases \(H_{\min}=CL^2\).
2. If \(d=0\), every pair of nonnegative block loadings satisfying
   \(\alpha^2+\beta^2=L\) is a minimizer, with \(H_{\min}=CL^2\).
   The two blocks can have independent signs. For \(m=2\), a zero block
   additionally permits the retained-pair family in item 1.
3. If \(d<0\), every minimizer is block-uniform with
   \(\alpha^2=\beta^2=L(m-1)/(m-1+\omega m)\).
   The block signs are independent, and (W.6) gives the minimum at zero.

The proof below establishes these endpoint qualifications separately.
It makes no claim about all nonglobal local minima.

## 1. Signs and the relaxed problem

The original domain is compact and \(H\) is continuous, so its minimum exists.
Replacing loadings by their absolute values cannot increase any residual:
\[
(r+z)^2-(r-z)^2=4rz\ge0\qquad(r,z\ge0).
\]
For \(\phi>0\), every target correlation and every weight is positive.
Thus opposite signs among nonzero minimizing loadings would strictly
increase the objective. The zero vector is not a minimum: assigning two
items in one block a small positive loading decreases their pair residual.
At a zero coordinate of a nonnegative, nonzero vector,
\[
\partial_iH=-2\sum_{j\ne i}w_{ij}R_{ij}l_j<0,
\]
where \(w_{ij}=1\) within blocks and \(w_{ij}=\omega\) across blocks.
Consequently every nonnegative minimizer is strictly positive when \(\phi>0\).

For nonnegative block loadings \(x_i,y_i\), define
\[
S_A=\sum_i x_i,\qquad Q_A=\sum_i x_i^2,\qquad
\xi=\frac{S_A^2-Q_A}{K},
\]
and define \(S_B,Q_B,\eta\) analogously. The quantities \(\xi,\eta\)
are the mean within-block pair products. Cauchy-Schwarz gives
\[
Q_A\ge m\xi,\qquad Q_B\ge m\eta,\qquad
S_A^2=Q_A+K\xi,\qquad S_B^2=Q_B+K\eta.
\tag{W.7}
\]
Jensen's inequality for each block's squared residuals gives
\[
H\ge f(\xi,\eta,Q_A,Q_B),
\]
\[
f=C[(L-\xi)^2+(L-\eta)^2]
+\omega\left[Q_AQ_B-2c\sqrt{(Q_A+K\xi)(Q_B+K\eta)}+m^2c^2\right].
\tag{W.8}
\]
The cross-block expression is exact. Enlarge the domain to
\[
\xi,\eta\ge0,\qquad Q_A\ge m\xi,\qquad Q_B\ge m\eta,
\tag{W.9}
\]
without upper bounds. When both norm inequalities are equalities, the
uniform loadings \(x_i=\sqrt\xi\), \(y_i=\sqrt\eta\) attain (W.8).
The loading bounds will be checked after minimization.

## 2. The relaxed minimum is attained for positive correlation

The constraints give \(S_A^2\le mQ_A\) and \(S_B^2\le mQ_B\). Hence
\[
f\ge C[(L-\xi)^2+(L-\eta)^2]
+\omega(\sqrt{Q_AQ_B}-mc)^2.
\tag{W.10}
\]
Thus \(\xi,\eta\), and \(T=Q_AQ_B\) remain bounded on every bounded
sublevel. An unbounded sequence in such a sublevel has a subsequence with
\(Q_A\to\infty\), \(Q_B\to0\), and \(\eta\to0\), after block
interchange if necessary. Along this subsequence,
\[
S_AS_B\le\sqrt{mT(1+K\xi/Q_A)}.
\]
Take further subsequences along which the remaining bounded quantities
converge. Dropping the nonnegative limiting term \(C(L-\xi)^2\) yields
\[
\begin{aligned}
\liminf f
&\ge CL^2+\omega\inf_{T\ge0}
\left[T-2c\sqrt{mT}+m^2c^2\right]\\
&=CL^2+\omega m(m-1)c^2=:E_\omega.
\end{aligned}
\tag{W.11}
\]

A finite feasible point has value strictly below \(E_\omega\) whenever
\(\phi>0\). This assertion uses only direct candidate evaluations.
For uniform block loadings, set \(u=\alpha^2+\beta^2\),
\(v=\alpha\beta\). Their objective is
\[
h(u,v)=C(u-L)^2+mdv^2-2\omega m^2cv+CL^2+\omega m^2c^2,
\qquad u\ge2v\ge0.
\tag{W.12}
\]

If \(d>0\) and \(0<\phi\le\phi_\omega\), the candidate
\((u,v)=(L,\omega mc/d)\) is feasible. Its value is (W.5), and
\[
E_\omega-h(L,\omega mc/d)
=\omega m(m-1)c^2\left(1+\frac m d\right)>0.
\tag{W.13}
\]
If \(d>0\) and \(\phi\ge\phi_\omega\), use the symmetric candidate
\((u,v)=(2t_s,t_s)\), whose value is (W.6).
At \(\phi_\omega\), this candidate equals the previous one.
As \(\phi\) increases from this threshold, \(E_\omega\) strictly
increases and the candidate value strictly decreases. Their gap remains
positive.

If \(d\le0\), use the same symmetric candidate. Its value at zero obeys
\[
h(2t_s,t_s)\big|_{\phi=0}
=\frac{\omega m^2(m-1)}{m-1+\omega m}L^2\le CL^2.
\]
Its value strictly decreases for \(0<\phi<1\), whereas
\(E_\omega>CL^2\). This also gives the required strict gap.

Therefore a minimizing sequence cannot escape to infinity. The relaxed
domain is closed, and a bounded minimizing sequence has a convergent
subsequence. Continuity of (W.8) proves attainment.

## 3. Both relaxed norm constraints are tight

Assume \(\phi>0\), and consider a relaxed minimizer. Both norms are
positive. If \(Q_A=0\) and \(Q_B>0\), then \(\xi=0\), and the
feasible perturbation \((\xi,Q_A)=(t^2,mt^2)\) changes the objective by
\[
-2\omega mc\,t\sqrt{Q_B+K\eta}+O(t^2)<0
\]
for small \(t>0\). If both norms vanish, small equal positive loadings
in both blocks improve the value. The argument is symmetric for \(Q_B=0\).
Also \(\xi,\eta>0\). At \(\xi=0\), \(Q_A>0\), a small increase in
\(\xi\) is feasible, and
\[
\frac1K\partial_\xi f=\xi-L-\omega c\frac{S_B}{S_A}<0.
\tag{W.14}
\]

Suppose \(Q_A>m\xi\). The variables \(Q_A,\xi\) are then locally
free, so their stationarity equations give
\[
Q_B=c\frac{S_B}{S_A},\qquad
\xi=L+\omega c\frac{S_B}{S_A}=L+\omega Q_B>L.
\tag{W.15}
\]

Both norm constraints cannot be slack. If \(Q_B>m\eta\), the corresponding
equations give \(\eta>L\) and \(Q_A=cS_A/S_B\). Consequently
\[
Q_AQ_B=c^2,
\qquad Q_AQ_B>m^2\xi\eta>m^2L^2,
\]
which contradicts \(c=L\phi<L\). Hence \(Q_B=m\eta\).

Vary \(\eta\) along this boundary, keeping \(Q_B=m\eta\).
This variation is two-sided because \(\eta>0\). Stationarity requires
\[
0=K(\eta-L)+\omega mQ_A-\omega m^2c\frac{S_A}{S_B}.
\tag{W.16}
\]
Since \(S_B^2=m^2\eta\), (W.15) gives
\[
S_A^2=\frac{c^2}{\eta},\qquad
c\frac{S_A}{S_B}=\frac{S_A^2}{m}.
\]
Substitute \(S_A^2=Q_A+K\xi\) into (W.16). The result is
\[
\eta=L+\omega m\xi>L.
\tag{W.17}
\]
But feasibility also gives \(S_A^2\ge m^2\xi\). Thus
\[
c^2=\eta S_A^2\ge m^2\xi\eta>m^2L^2,
\tag{W.18}
\]
again contradicting \(c<L\). This contradiction holds for every
\(\omega>0\). It does not require \(\omega m\ge1\).

Therefore \(Q_A=m\xi\), and the symmetric argument gives \(Q_B=m\eta\).
Every relaxed minimizer is attained by uniform loadings within both blocks.
In particular, the relaxed minimum equals the minimum of (W.12) on its cone.

## 4. Exact minimization on the cone

Nonnegative \(\alpha,\beta\) map onto \(u\ge2v\ge0\), because
\[
\alpha=\frac{\sqrt{u+2v}+\sqrt{u-2v}}2,\qquad
\beta=\frac{\sqrt{u+2v}-\sqrt{u-2v}}2
\]
recovers the solution with \(\alpha\ge\beta\).

For \(d>0\), (W.12) is strictly convex in \((u,v)\). Its unconstrained
minimum is \((L,\omega mc/d)\), feasible exactly when
\(\phi\le\phi_\omega\). Substitution proves (W.4) and (W.5).
For \(\phi\ge\phi_\omega\), take \((u,v)=(2t_s,t_s)\) on the edge
\(u=2v\). The constraint \(v\ge0\) is inactive. For the Lagrangian
\(h+\varpi(2v-u)\), stationarity holds with
\[
\varpi=2C(2t_s-L)
=\frac{m(m-1)L(2\omega m\phi-d)}{m-1+\omega m}.
\tag{W.19}
\]
This multiplier is nonnegative precisely when
\(\phi\ge\phi_\omega\). Convexity makes these conditions sufficient;
strict convexity gives uniqueness in \((u,v)\).
At equality the multiplier is zero, and both branches give
\((u,v)=(L,L/2)\). The loading magnitudes are uniquely
\(\alpha=\beta=\sqrt{L/2}\).

For \(d\le0\) and \(\phi>0\),
\[
\partial_v h=2mdv-2\omega m^2c<0.
\]
For each fixed \(u\), the minimum therefore occurs at \(v=u/2\).
Write \(u=2t\). On that edge,
\[
h(2t,t)=m(m-1)(L-t)^2+\omega m^2(c-t)^2.
\tag{W.20}
\]
This strictly convex quadratic has the unique minimum \(t=t_s\).
Its value is (W.6). The same edge calculation evaluates the value for
\(d>0\) above the threshold.

The asymmetric solutions have \(\alpha^2+\beta^2=L\), with both entries
positive for \(\phi>0\). The symmetric solution satisfies
\(0<t_s<L\) because \(\phi<1\). All resulting loadings therefore satisfy
the original bounds strictly.

For any original minimizer, its objective is at least the relaxed value,
which is attained by these bounded candidates. Equality forces its relaxed
coordinates to minimize (W.8). Section 3 then gives equality in both norm
constraints. Equality in Cauchy-Schwarz in (W.7) forces all loadings within
each block to coincide. This proves the unrestricted global reduction,
including exclusion of loading-boundary global minima for \(\phi>0\).

When \(d>0\), the positive asymmetric branch also satisfies
\[
(\alpha-\beta)^2
=L-\frac{2\omega mL\phi}{d}
=\frac{2\omega mL}{d}(\phi_\omega-\phi).
\tag{W.21}
\]
Thus the two branches merge continuously with the stated square-root form.

## 5. Zero correlation and all equality cases

Set \(\phi=c=0\). Absolute values still cannot increase the objective.
The within-block target is positive, so all nonzero loadings in a retained
block have a common sign. The cross-block target is zero, so the two block
signs are independent.

From (W.8) and (W.9),
\[
\begin{aligned}
H
&\ge C[(L-\xi)^2+(L-\eta)^2]+\omega Q_AQ_B\\
&\ge C[(L-\xi)^2+(L-\eta)^2]+\omega m^2\xi\eta\\
&=C(u-L)^2+mdv^2+CL^2,
\end{aligned}
\tag{W.22}
\]
where \(u=\xi+\eta\) and \(v=\sqrt{\xi\eta}\), so \(u\ge2v\ge0\).
Uniform block candidates attain the minima of this lower bound.

For \(d>0\), the last expression has its unique cone minimum at
\((u,v)=(L,0)\). Hence \((\xi,\eta)=(L,0)\) or \((0,L)\).
Choose \((L,0)\). Then \(Q_A>0\), so equality in (W.22) requires
\(Q_B=0\). Equality in Jensen's inequality requires every retained
within-block product to equal \(L\). For \(m\ge3\), positivity of these
products forces each retained magnitude to equal \(\sqrt L\).
For \(m=2\), only one product is fixed, giving exactly
\[
(|l_1|,|l_2|)=(r,L/r),\qquad L\le r\le1.
\tag{W.23}
\]
The two retained signs agree. This is the complete endpoint family, including
the two endpoint magnitude pairs with one loading magnitude equal to one.

For \(d=0\), the lower bound is minimized exactly when \(u=L\), with
arbitrary \(0\le v\le L/2\). If \(\xi,\eta>0\), equality requires
\(Q_AQ_B=m^2\xi\eta\). Since each norm separately satisfies (W.9), both
norm constraints must be equalities. Both blocks are therefore uniform,
with \(\alpha^2+\beta^2=L\). If one mean pair product is zero, equality
forces that entire block to vanish. The retained block then has exactly
the classification in (W.23), or equal magnitudes for \(m\ge3\).
Thus the zero-correlation weight boundary has a continuum of minimizing
covariances; it has no unique symmetric solution.

For \(d<0\), every fixed \(u>0\) has its minimum at \(v=u/2\), so
\(\xi=\eta\). Minimizing the resulting strictly convex quadratic gives
\[
\xi=\eta=\frac{L(m-1)}{m-1+\omega m}>0.
\]
Equality again forces both norm constraints to be tight. Thus both blocks
are uniform with these equal magnitudes, even for \(m=2\). Independent
block signs give two fitted covariances, distinguished by the sign of their
cross-block fitted correlations. Their objective values agree because
the cross-block target is zero.

These cases exhaust (W.22)'s equality conditions and prove every endpoint
assertion in the theorem.

## 6. Relation to the unweighted result and limits of the statement

At \(\omega=1\), \(d=1\) and \(\phi_\omega=1/(2m)\).
Equations (W.4) through (W.6) reduce to the ULS formulas of Theorem 2 of the
article, proved in Section 11 of Online Resource 1. The contradiction (W.18) supplies the
additional argument needed for arbitrary positive weights.

The theorem applies pointwise if the pair-class weights depend on the
population parameters. It does not establish how such weights vary, nor
does it establish the number of solutions to an implicit threshold equation.
It assumes identical weights within each pair class. General item-specific
weights and a complete classification of nonglobal local minima are outside
its scope. Zero weights, \(L=0\), \(L=1\), and \(\phi=1\) are also outside
the stated assumptions.

## 7. Reproducible exact algebra checks

The following SymPy snippet checks the displayed polynomial and rational
identities. In it, `k` denotes the weight ratio \(\omega\) and `mu` the
multiplier \(\varpi\). It supplements the analytic proof; it does not replace its
inequality, attainment, or equality-case arguments.

```python
import sympy as s

m, k, L, phi = s.symbols("m k L phi", positive=True)
a, b, u, v, r, xi, eta, QA = s.symbols("a b u v r xi eta QA")
n = m - 1
C, K, d, c = m*n/2, m*n, k*m-n, L*phi
h = C*(u-L)**2 + m*d*v**2 - 2*k*m**2*c*v + C*L**2 + k*m**2*c**2
direct = C*((L-a**2)**2 + (L-b**2)**2) + k*m**2*(c-a*b)**2
ts = L*(n+k*m*phi)/(n+k*m)
threshold = d/(2*k*m)
Ga = C*L**2 - k*m**2*n*c**2/d
Gs = k*m**2*n*L**2*(1-phi)**2/(n+k*m)
E = C*L**2 + k*m*n*c**2
mu = m*n*L*(2*k*m*phi-d)/(n+k*m)
edge = {u: 2*ts, v: ts}

checks = {
    "cone expansion": direct-h.subs({u: a**2+b**2, v: a*b}),
    "asymmetric value": h.subs({u: L, v: k*m*c/d})-Ga,
    "symmetric value": h.subs(edge)-Gs,
    "edge stationarity u": s.diff(h, u).subs(edge)-mu,
    "edge stationarity v": s.diff(h, v).subs(edge)+2*mu,
    "threshold loading": ts.subs(phi, threshold)-L/2,
    "threshold value": (Ga-Gs).subs(phi, threshold),
    "escape gap": E-Ga-k*m*n*c**2*(1+m/d),
    "escape square": r**2-2*c*s.sqrt(m)*r+m**2*c**2
                    -(r-s.sqrt(m)*c)**2-m*n*c**2,
    "boundary stationarity": K*(eta-L)+k*m*QA-k*m*(QA+K*xi)
                             -K*(eta-L-k*m*xi),
    "endpoint lower bound": C*((L-xi)**2+(L-eta)**2)+k*m**2*xi*eta
                            -C*(xi+eta-L)**2-m*d*xi*eta-C*L**2,
    "loading gap": L-2*k*m*c/d-(2*k*m*L/d)*(threshold-phi),
    "ULS specialization": h.subs(k, 1)
                          -(C*(u-L)**2+m*(v-m*c)**2+C*L**2-m**2*n*c**2),
}
for name, difference in checks.items():
    assert s.factor(difference) == 0, name
print(f"{len(checks)} exact identities passed")
```
