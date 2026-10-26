---
title: 'Online Resource 1 for “Fitting One Factor to Two Traits: Population Solutions and Fit Indices”'
subtitle: 'Symmetry breaking in a one-factor maximum likelihood fit to two correlated traits'
---

## Scope

This Online Resource gives the complete proof of Theorem 1 of the article, the maximum likelihood (ML) solution for two exchangeable blocks, whose main steps are summarized in Appendix A of the article. The statement of Theorem 1 in Section 1 also includes the behavior at the transition that the article states separately as Corollary 1 (continuity and \(C^1\)-but-not-\(C^2\) regularity of the minimum discrepancy, the loading-gap identity (1.9) and the square-root asymptotic (1.10)); these are proved in Lemma 7.2. The resource also proves the zero-uniqueness boundary criterion of Proposition 4 (Lemma 4 here; Appendix B of the article), the reflection and equal cross-loading reductions used for Proposition 3 (Section 10 here; Section 5.3 of the article), and the global unweighted least-squares (ULS) result of Theorem 2 (Section 11 here; Appendix C of the article). Proposition 1 of the article, the weighted least-squares result (Section 5.1), is proved in Online Resource 2, and blocks that differ in size or loading are treated in Online Resource 3. Section 12 reports numerical checks, which are separate from the proofs, and Section 13 lists the assumptions, exceptional cases and limits of the statements.

**Notation.** Symbols follow the article where they denote the same objects. Two exceptions are kept for traceability with the verification scripts `verify_symmetry_breaking.py` and `check_saddle_7_3_sympy.py`: \(N=2m-1=p-1\) is the article's \(\nu\) (the article uses \(N\) for sample size), and \(\phi_*\), \(q_*\) are the article's \(\varphi^*\), \(q^*\). The saddle script writes the block uniquenesses \(\psi_A,\psi_B\) of Sections 6 and 7 as \(x,y\). Section, lemma, corollary and equation numbers are local to this resource: Corollaries 8.1, C1 and C2 are distinct from Corollary 1 of the article, and equations (A.1) and (A.2) belong to the appendix at the end of this resource, not to Appendix A of the article.

## Status of the results

**Theorem 1 is proved for every integer \(m\ge2\), \(0<\lambda<1\), and \(0\le\phi<1\), with the stated two-item endpoint exception interpreted as a family of parameter representations.** No numerical calculation or computer-assisted certificate is used in the proof.

Beyond the value of the minimum discrepancy, the theorem gives closed forms for the asymmetric branch and its discrepancy, and for \(0<\phi<1\) every local minimum is global. The minimum discrepancy is \(C^1\), but not \(C^2\), at the threshold (Corollary 1 of the article).

The reflection to negative trait correlation (Corollary C1) is proved. The equal cross-loading reduction (Corollary C2) needs an explicit normalization: the formula \(\sigma_S^2/(1+\sigma_S^2)\) for the effective squared loading presupposes unstandardized residual variance one. The general statement is proved in Section 10. For ULS, both the local threshold and a full global theorem are proved. A complete classification of nonglobal ULS local minima is not claimed.

All factor-analysis identities used below are rederived. The only background tools are elementary finite-dimensional linear algebra, Perron's theorem for positive matrices (Horn & Johnson, 2013, Ch. 8), differentiation, compactness, Cauchy–Schwarz, and convexity of a squared residual. This resource makes no literature or novelty claim.

## 1. Notation and statement of Theorem 1

Throughout, write
\[
L=\lambda^2,\quad a=1-L,\quad n=m-1,\quad N=2m-1=p-1,
\]
\[
b=a+mL=1+nL,\qquad h=mL\phi,\qquad
q=\frac hb=\frac{mL\phi}{b},\qquad q_* =\frac{L}{2-L}.
\tag{1.1}
\]
Thus \(a>0\), \(b>a\), \(0\le h<b\), and \(0\le q<1\). Set
\[
\phi_* =\frac{b}{m(2-L)},\qquad
F_0=-n\log a-\log b.
\tag{1.2}
\]
Let \(R\) have unit diagonal, within-block off-diagonal entries \(L\), and between-block entries \(L\phi\). We minimize
\[
F(\Sigma)=\operatorname{tr}(\Sigma^{-1}R)
 +\log\det\Sigma-\log\det R-p
\tag{1.3}
\]
over \(\Sigma=ll^T+\Psi\succ0\), with \(\Psi\) diagonal and nonnegative.

### Theorem 1 — Exact ML minimum and minimizers. **Proved.**

The infimum is attained and
\[
\boxed{
F_{\min}(\phi)=
\begin{cases}
F_0+n\log(1-q^2),&0\le\phi\le\phi_*,\\[3pt]
N\log\!\left(1+\dfrac{mL(1-\phi)}{aN}\right)
-\log\!\left(1+\dfrac{mL(1-\phi)}a\right),
&\phi_*\le\phi<1.
\end{cases}}
\tag{1.4}
\]
The two expressions coincide at the threshold.

For \(\phi\ge\phi_*\), the minimizing covariance is unique. Its parameters, up to a simultaneous sign change of all loadings, are
\[
l_i=l_s,\quad
l_s^2=\frac{L(n+m\phi)}N,\quad
\psi_i=\psi_s=a+\frac{mL(1-\phi)}N=1-l_s^2.
\tag{1.5}
\]

For \(0<\phi<\phi_*\), there are exactly two minimizing covariances, interchanged by swapping the blocks. Choose nonnegative loadings and the block order so that \(\alpha>\beta\). Then
\[
l_i=\alpha\ (i\in A),\qquad l_i=\beta\ (i\in B),
\qquad \psi_i=1-l_i^2,
\]
\[
\boxed{\alpha\beta=q,\qquad
\alpha^2+\beta^2=L+(2-L)q^2.}
\tag{1.6}
\]
Equivalently, if \(S=L+(2-L)q^2\),
\[
\boxed{
\alpha=\frac{\sqrt{S+2q}+\sqrt{S-2q}}2,
\qquad
\beta=\frac{\sqrt{S+2q}-\sqrt{S-2q}}2.}
\tag{1.7}
\]
Both uniquenesses are strictly positive. The only parameter symmetries are simultaneous sign change and block interchange.

At \(\phi=0\), exactly one block is fitted exactly and the other is the identity. For \(m\ge3\), this means retained loadings \(\lambda\), retained uniquenesses \(a\), and zero loadings with uniqueness one in the other block, up to the same symmetries. For \(m=2\), the retained block has the complete family
\[
\boxed{
(l_1,l_2)=(t,L/t),\quad L\le t\le1,\qquad
(\psi_1,\psi_2)=(1-t^2,1-L^2/t^2),}
\tag{1.8}
\]
with the other block's loadings zero and uniquenesses one, and with either block retained. Include simultaneous negatives. The endpoints in (1.8) have one zero uniqueness. All values of \(t\) give the same covariance for a chosen retained block: there are two minimizing covariances, although infinitely many parameter representations.

For \(0<\phi<1\), these global minimizers exhaust all local minima. The symmetric stationary solution is a saddle for every \(0\le\phi<\phi_*\). At \(\phi_*\) it is the unique minimizing covariance, with a quartic degeneracy in the symmetry-breaking direction.

Finally,
\[
\boxed{|\alpha-\beta|^2=(2-L)(q_*-q)(1-q),}
\tag{1.9}
\]
so, as \(\phi\uparrow\phi_*\),
\[
|\alpha-\beta|\sim C(\phi_*-\phi)^{1/2},
\qquad
\boxed{C=\lambda\sqrt{\frac{2m(1-\lambda^2)}{1+(m-1)\lambda^2}}.}
\tag{1.10}
\]
The function \(F_{\min}\) is continuous on \([0,1)\), real-analytic away from \(\phi_*\), and \(C^1\) but not \(C^2\) at \(\phi_*\).

## 2. Spectrum, coercivity, and existence

### Lemma 2.1 — Spectrum of the population matrix. **Proved.**

On the subspace of vectors whose coordinates sum to zero separately within each block, \(R\) acts as \(aI\). This subspace has dimension \(2n\). In the orthonormal block-mean basis \((1_A/\sqrt m,1_B/\sqrt m)\), the remaining matrix is
\[
\begin{pmatrix}b&h\\h&b\end{pmatrix}.
\]
Consequently the eigenvalues of \(R\) are \(a\), with multiplicity \(2n\), and \(b+h,b-h\). They are strictly positive, since \(b-h=a+mL(1-\phi)>0\), and
\[
\boxed{\det R=a^{2n}(b^2-h^2)=a^{2n}b^2(1-q^2).}
\tag{2.1}
\]

### Lemma 2.2 — Existence, including boundary representations. **Proved.**

Let \(g(t)=t-\log t-1\) for \(t>0\). Its derivative is \(1-1/t\), so \(g\ge0\), and \(g(t)\to\infty\) at both zero and infinity. If \(\mu_i\) are the eigenvalues of \(R^{1/2}\Sigma^{-1}R^{1/2}\), then
\[
F(\Sigma)=\sum_i g(\mu_i).
\]
On any finite sublevel \(F\le M\), every \(\mu_i\) lies in a fixed compact subinterval of \((0,\infty)\). Therefore \(\Sigma\) is uniformly bounded above and uniformly positive definite on that sublevel.

The model is closed relative to the positive-definite cone. Indeed, if \(\Sigma_j=l_jl_j^T+\Psi_j\to\Sigma\succ0\), then
\[
(l_j)_i^2\le(\Sigma_j)_{ii},\qquad
0\le(\Psi_j)_{ii}\le(\Sigma_j)_{ii}.
\]
Thus the parameters have a convergent subsequence, whose limit represents \(\Sigma\) with nonnegative uniquenesses. A nonempty finite sublevel is consequently compact. For example, \(I\) belongs to the model and has finite discrepancy. Taking a minimizing sequence in such a sublevel proves attainment. \(\square\)

## 3. Concentration and interior stationarity

### Lemma 3.1 — Concentration for a fixed positive uniqueness matrix. **Proved.**

Let \(\Psi\succ0\) be diagonal and put
\[
B=\Psi^{-1/2}R\Psi^{-1/2},\qquad v=\Psi^{-1/2}l,
\qquad r=v^Tv.
\]
The rank-one inverse and determinant identities give
\[
F=\operatorname{tr}B-\log\det B-p
-\frac{v^TBv}{1+r}+\log(1+r).
\tag{3.1}
\]
At fixed \(r\), minimize by choosing a leading eigenvector of \(B\). If its eigenvalues are \(\zeta_1\ge\cdots\ge\zeta_p\), so that the largest is \(\zeta_1\), the remaining derivative is
\[
\frac{d}{dr}\left[-\frac{\zeta_1r}{1+r}+\log(1+r)\right]
=\frac{1+r-\zeta_1}{(1+r)^2}.
\]
Thus \(r=\max(\zeta_1-1,0)\), and
\[
\min_l F=
\begin{cases}
\sum_{j=2}^p g(\zeta_j),&\zeta_1>1,\\
\sum_{j=1}^p g(\zeta_j),&\zeta_1\le1.
\end{cases}
\tag{3.2}
\]
At \(\zeta_1=1\), the two expressions agree and the minimizing loading is zero. Every nonzero local minimum over \(l\) also selects a leading eigendirection: a small rotation toward an eigendirection with a larger eigenvalue would otherwise decrease (3.1).

### Lemma 3.2 — Interior score equations and diagonal reproduction. **Proved.**

Write \(K=\Sigma^{-1}\) and
\[
H=K-KRK=K(\Sigma-R)K.
\]
At an interior stationary point, differentiation with respect to \(l\) and each uniqueness gives
\[
Hl=0,\qquad H_{ii}=0.
\tag{3.3}
\]
Put \(u=\Psi^{-1}l\) and \(\gamma=l^T\Psi^{-1}l\). Since \(Kl=u/(1+\gamma)\), multiplying \(Hl=0\) by \(\Sigma\) yields
\[
\boxed{Ru=(1+\gamma)l.}
\tag{3.4}
\]
Also \(\Sigma u=(1+\gamma)l\), so \((\Sigma-R)u=0\). Substituting \(K=\Psi^{-1}-uu^T/(1+\gamma)\) into \(H\) now gives
\[
H=\Psi^{-1}(\Sigma-R)\Psi^{-1}.
\]
The diagonal equations in (3.3) therefore imply
\[
\boxed{\Sigma_{ii}=1,\qquad\psi_i=1-l_i^2.}
\tag{3.5}
\]
Conversely, (3.4) and (3.5) imply (3.3), so they are sufficient for interior stationarity as well.

The equivalent row form, used below, is
\[
\boxed{\sum_{j\ne i}u_j(R_{ij}-l_i l_j)=0.}
\tag{3.6}
\]
To verify it, subtract the diagonal contribution from (3.4) and use \(u_i(1-l_i^2)=l_i\).

### Lemma 3.3 — Positivity at an interior local minimum when \(\phi>0\). **Proved.**

If \(l=0\), uniqueness stationarity gives \(\Psi=I\). But \(\lambda_{\max}(R)=b+h>1\), so (3.1) decreases for a sufficiently small nonzero loading in that leading direction. Hence \(l=0\) is not a local minimum.

For \(l\ne0\), Lemma 3.1 requires \(v=\Psi^{-1/2}l\) to be a leading eigenvector of \(B=\Psi^{-1/2}R\Psi^{-1/2}\). When \(\phi>0\), every entry of \(B\) is strictly positive. By Perron's theorem for positive matrices (Horn & Johnson, 2013, Ch. 8), the spectral radius \(\rho(B)\) is then an algebraically simple eigenvalue, with an eigenvector whose entries are all strictly positive. Since \(B\) is symmetric, its eigenvalues are real and bounded in modulus by \(\rho(B)\). So \(\rho(B)=\zeta_1\) is the leading eigenvalue, and its one-dimensional eigenspace is spanned by that positive vector. Hence \(v\), and therefore \(l=\Psi^{1/2}v\), is strictly positive up to a global sign change. Equation (3.5) then gives \(0<l_i<1\). \(\square\)

## 4. Heywood boundary points

### Lemma 4: Boundary local minima are excluded except at \(m=2,\phi=0\). **Proved.**

**Step 0: at most one zero uniqueness.** A positive-definite one-factor covariance can have at most one zero uniqueness. If \(\psi_j=\psi_k=0\) with \(j\ne k\), the principal submatrix of \(\Sigma\) on \(\{j,k\}\) equals that of the rank-one matrix \(ll^T\), which is singular.

**Step 1: coordinates on the stratum.** Fix \(j\) and consider the stratum \(\psi_j=0\). Positive definiteness gives \(\Sigma_{jj}=l_j^2>0\). Fixing the sign of \(l_j\), the map
\[
(l,\ \psi_{-j})\ \longmapsto\ (t,\delta,\psi_{-j}),\qquad
t=l_j^2,\quad \delta_i=l_i/l_j\ (i\ne j),
\]
is a diffeomorphism of the stratum onto \((0,\infty)\times\mathbb R^{p-1}\times(0,\infty)^{p-1}\), with the other uniquenesses positive by Step 0.

**Step 2: factorization of the stratum covariance.** Let \(D=\operatorname{diag}(d_1,\dots,d_p)\) with \(d_j=t\) and \(d_i=\psi_i\) for \(i\ne j\), and let
\[
T=I+\sum_{i\ne j}\delta_i\,e_ie_j^T .
\]
If \(\xi\sim N(0,D)\) and \(X=T\xi\), then \(X_j=\xi_j\) and \(X_i=\delta_iX_j+\xi_i\). Hence
\[
\Sigma=TDT^T,
\]
because \(\operatorname{Cov}(X_j)=t=l_j^2\), \(\operatorname{Cov}(X_i,X_j)=\delta_it=l_il_j\), and \(\operatorname{Cov}(X_i,X_k)=\delta_i\delta_kt+\mathbf 1\{i=k\}\psi_i=l_il_k+\mathbf 1\{i=k\}\psi_i\). Since \((e_ie_j^T)(e_ke_j^T)=0\) for \(i,k\ne j\), we have \(T^{-1}=I-\sum_{i\ne j}\delta_ie_ie_j^T\) and \(\det T=1\). Therefore
\[
\log\det\Sigma=\log t+\sum_{i\ne j}\log\psi_i,
\qquad
x^T\Sigma^{-1}x=\frac{x_j^2}{t}+\sum_{i\ne j}\frac{(x_i-\delta_ix_j)^2}{\psi_i}.
\tag{4.0}
\]

**Step 3: the stratum objective and its unique minimizer.** Taking the expectation of (4.0) under \(x\sim N(0,R)\) gives \(\operatorname{tr}(\Sigma^{-1}R)\). Hence, up to the constant \(-\log\det R-p\),
\[
F=\log t+\frac1t+\sum_{i\ne j}\left[\log\psi_i+\frac{1-2\delta_iR_{ij}+\delta_i^2}{\psi_i}\right].
\tag{4.1}
\]
The objective is separable. In \(t\), the derivative \(1/t-1/t^2\) vanishes only at \(t=1\). For each \(i\), the \(\delta_i\)-derivative \(2(\delta_i-R_{ij})/\psi_i\) vanishes only at \(\delta_i=R_{ij}\), and the \(\psi_i\)-derivative \(1/\psi_i-c_i/\psi_i^2\), with \(c_i=1-2\delta_iR_{ij}+\delta_i^2\ge 1-R_{ij}^2>0\), vanishes only at \(\psi_i=c_i\). So the stratum has a unique stationary point:
\[
t=1,\quad\delta_i=R_{ij},\quad\psi_i=1-R_{ij}^2 .
\tag{4.2}
\]
Each separable term is coercive and has a single stationary point, so (4.2) is also the stratum's global minimizer. A boundary local minimum of the full problem is in particular a local minimum on its stratum, so it must be the **star covariance** \(\Sigma_\star\): choose \(l_j=1\), so \(l_i=R_{ij}\) and \(\psi_i=1-R_{ij}^2\) for \(i\ne j\), with \(\psi_j=0\).

**Step 4: the inward derivative.** Increasing \(\psi_j\) from zero is a feasible move, \(\Sigma_\star\mapsto\Sigma_\star+\epsilon e_je_j^T\) with \(\epsilon\ge0\). Since \(\Sigma_\star\succ0\), \(F\) is differentiable in a neighbourhood of \(\Sigma_\star\) in the space of symmetric matrices, and with \(K=\Sigma_\star^{-1}\),
\[
\frac{d}{d\epsilon}F(\Sigma_\star+\epsilon e_je_j^T)\Big|_{\epsilon=0^+}
=\operatorname{tr}\big[(K-KRK)e_je_j^T\big]=K_{jj}-(KRK)_{jj}.
\]
Write
\[
w_i=\frac{R_{ij}}{1-R_{ij}^2}\ (i\ne j),\qquad
s=\sum_{i\ne j}w_iR_{ij}=\sum_{i\ne j}\frac{R_{ij}^2}{1-R_{ij}^2}.
\]
Reading the quadratic form (4.0) at \(t=1\), \(\delta_i=R_{ij}\), the \(j\)-th column of \(K\) is
\[
Ke_j=(1+s)\,e_j-\sum_{i\ne j}w_ie_i,
\]
so \(K_{jj}=1+s\). Using \(R_{jj}=R_{ii}=1\),
\[
(KRK)_{jj}=(Ke_j)^TR\,(Ke_j)
=(1+s)^2-2(1+s)s+\sum_{i\ne j}w_i^2+\sum_{\substack{i\ne k\\ i,k\ne j}}w_iw_kR_{ik}.
\]
Hence
\[
K_{jj}-(KRK)_{jj}
=(1+s)\big[1-(1+s)+2s\big]-\sum_{i\ne j}w_i^2-\sum_{i\ne k}w_iw_kR_{ik}
=s+s^2-\sum_{i\ne j}w_i^2-\sum_{i\ne k}w_iw_kR_{ik}.
\]
Expand \(s^2=\sum_{i\ne j}w_i^2R_{ij}^2+\sum_{i\ne k}w_iw_kR_{ij}R_{kj}\). The single-index terms cancel exactly:
\[
s+\sum_{i\ne j}w_i^2R_{ij}^2-\sum_{i\ne j}w_i^2
=\sum_{i\ne j}\frac{R_{ij}^2}{(1-R_{ij}^2)^2}\Big[(1-R_{ij}^2)+R_{ij}^2-1\Big]=0.
\]
What remains is
\[
\boxed{
\left.\frac{d}{d\epsilon}F(\Sigma_\star+\epsilon e_je_j^T)\right|_{0^+}
=-\sum_{\substack{i\ne k\\i,k\ne j}}
\frac{R_{ij}R_{kj}\,(R_{ik}-R_{ij}R_{kj})}
{(1-R_{ij}^2)(1-R_{kj}^2)}.}
\tag{4.3}
\]
This holds for any positive-definite correlation matrix \(R\). Positive definiteness implies \(|R_{ij}|<1\) for \(i\ne j\) and makes the stated discrepancy finite; the two-block structure has not been used yet.

**Interpretation.** Let \(r_i=X_i-R_{ij}X_j\) be the residual of item \(i\) after regression on item \(j\). Under \(R\), \(\operatorname{Var}(r_i)=1-R_{ij}^2\) and \(\operatorname{Cov}(r_i,r_k)=R_{ik}-R_{ij}R_{kj}\). Equation (4.3) is therefore
\[
\left.\frac{d}{d\epsilon}F\right|_{0^+}
=\sum_{i\ne j}w_i^2\operatorname{Var}(r_i)-\operatorname{Var}\Big(\sum_{i\ne j}w_ir_i\Big).
\tag{4.4}
\]
The star covariance treats the residuals as independent. Moving item \(j\) off the boundary pays off exactly when the true residuals, weighted by \(w\), co-vary positively: conditioning on item \(j\) leaves shared variance that only a common factor distinct from \(X_j\) can absorb.

**Step 5: strict sign in the two-block population.** Take \(j\in A\) without loss of generality. Every weight product \(w_iw_k\) is nonnegative, and the residual covariances are
\[
L-L^2=aL\quad(i,k\in A),\qquad
L-L^2\phi^2>0\quad(i,k\in B),\qquad
L\phi-L^2\phi=aL\phi\quad(i\in A,\ k\in B),
\]
all nonnegative. Hence (4.3) is at most zero. It is strictly negative as soon as one summand is positive:
- if \(\phi>0\), any pair \(i\ne k\) in \(B\) has \(w_iw_k=L^2\phi^2/(1-L^2\phi^2)^2>0\) and positive residual covariance;
- if \(\phi=0\) and \(m\ge3\), any pair \(i\ne k\) in \(A\mathbin{\backslash}\{j\}\) has \(w_iw_k=L^2/(1-L^2)^2>0\) and residual covariance \(aL>0\).

In both cases increasing \(\psi_j\) strictly decreases \(F\), so \(\Sigma_\star\) is not a local minimum. Since \(\Sigma_\star\) was the only candidate on the stratum (Step 3), no boundary point is a local minimum.

Only \(m=2,\phi=0\) remains. Then \(A\mathbin{\backslash}\{j\}\) has one item, and every weight from \(B\) is zero, so (4.3) vanishes. The star covariance fits the two-item retained block exactly and sets the other block to the identity. It is an endpoint representation of (1.8); its global optimality and completeness are proved in Section 8. \(\square\)

*Numerical check (not part of the proof).* On 400 random correlation matrices with \(p=3,\dots,8\) and random \(j\), the closed-form precision (4.0), the intermediate identity, (4.3) and the variance form (4.4) matched direct matrix computation to \(4\times10^{-14}\). One-sided finite differences agreed to within the expected \(O(\epsilon)\) error.

## 5. Unrestricted within-block reduction

### Lemma 5.1 — Two distinct stationary loadings in a block must overfit their covariance. **Proved.**

Consider a positive interior stationary point with \(0<\phi<1\), and set \(H_0=1+\gamma\). For \(i\in A\), (3.4) becomes
\[
H_0l_i-a\frac{l_i}{1-l_i^2}=C_A,
\qquad C_A=L(U_A+\phi U_B)>0,
\quad U_A=\sum_{i\in A}u_i.
\tag{5.1}
\]
The same holds in \(B\). The function
\[
f(t)=H_0t-\frac{at}{1-t^2},\quad0<t<1,
\qquad
f''(t)=-\frac{2at(3+t^2)}{(1-t^2)^3}<0
\]
is strictly concave. Thus a horizontal line meets its graph at most twice, so a block has at most two distinct positive loading values.

Suppose two distinct values \(x<y\) occur in a block. Put
\[
z=xy,\quad P=(1-x^2)(1-y^2),\quad
\gamma_{\mathrm{rest}}=\sum_{k\ne i,j}\frac{l_k^2}{1-l_k^2}.
\]
The identity
\[
\frac{y/(1-y^2)-x/(1-x^2)}{y-x}
=\frac{1+xy}{(1-x^2)(1-y^2)}
\]
and \(f(x)=f(y)\) give \(H_0P=a(1+z)\). On the other hand,
\[
H_0P=1-z^2+\gamma_{\mathrm{rest}}P.
\]
Consequently
\[
\boxed{\gamma_{\mathrm{rest}}P=(1+z)(z-L).}
\tag{5.2}
\]
All loadings are positive and there are at least two other items, so \(\gamma_{\mathrm{rest}}>0\). It follows that
\[
\boxed{xy>L.}
\tag{5.3}
\]

### Lemma 5.2 — Every positive interior stationary point is block-uniform. **Proved.**

Suppose block \(A\) has its two values \(x<y\). By (5.3), \(xy>L\), hence \(y>\sqrt L\). In the score row (3.6) of an item loading \(y\), every within-block term is strictly negative: every fitted product is at least \(xy>L\). A cross-block residual must therefore be positive. If \(z=\min_{j\in B}l_j\), this implies
\[
yz<L\phi.
\tag{5.4}
\]
In particular \(z<\sqrt L\). In the row of an item loading \(z\), every cross-block residual is positive, because every loading in \(A\) is at most \(y\). If \(B\) were uniform, its within-block residuals \(L-z^2\) would also all be positive, contradicting (3.6). Hence \(B\) also has two values, say \(z<w\).

Applying (5.3) in \(B\) gives \(zw>L\). At a row loading \(w\), all within-block terms are negative, so some cross-block residual must be positive. Since \(x\) is the smallest loading in \(A\),
\[
wx<L\phi.
\tag{5.5}
\]
Multiplication of (5.4) and (5.5) gives
\[
(xy)(zw)<L^2\phi^2<L^2,
\]
contradicting \((xy)(zw)>L^2\). Neither block can be mixed. \(\square\)

**Consequence.** By Lemmas 3.3 and 4, every local minimum, and in particular every global minimum, for \(0<\phi<1\) is an interior block-uniform solution with strictly positive loadings after a global sign change. This proves the unrestricted reduction; it is not an assumption based on permutation symmetry or numerical evidence.

## 6. Exact global analysis after reduction

### Lemma 6.1 — Eigenvalue ordering and scale profiling. **Proved.**

Let \(\Psi=\operatorname{diag}(\psi_AI_m,\psi_BI_m)\), with \(\psi_A,\psi_B>0\). The whitened population matrix has contrast eigenvalues \(a/\psi_A,a/\psi_B\), each of multiplicity \(n\), and the two eigenvalues of
\[
\begin{pmatrix}b/\psi_A&h/\sqrt{\psi_A\psi_B}\\h/\sqrt{\psi_A\psi_B}&b/\psi_B\end{pmatrix}.
\tag{6.1}
\]
Its larger eigenvalue is at least \(\max(b/\psi_A,b/\psi_B)\), strictly larger than both contrast eigenvalues since \(b>a\). Therefore the largest eigenvalue always comes from (6.1). The only top-eigenvalue degeneracy is \(h=0,\psi_A=\psi_B\).

Write \(\psi_A=e^{s+d},\psi_B=e^{s-d}\), and define
\[
T_\pm(d)=b\cosh d\pm\sqrt{b^2\sinh^2d+h^2},
\qquad Q(d)=2na\cosh d+T_-(d).
\tag{6.2}
\]
Before multiplication by \(e^{-s}\), the leading eigenvalue is \(T_+\) and the sum of the other \(N\) eigenvalues is \(Q\). In the active region \(e^{-s}T_+>1\), Lemma 3.1 gives
\[
F_{\mathrm{conc}}=e^{-s}Q+Ns-2n\log a-\log T_- -N.
\tag{6.3}
\]
It is strictly convex in \(s\), and its unique optimum is
\[
\boxed{e^{s_*(d)}=Q(d)/N.}
\tag{6.4}
\]
This optimum is active: every nonleading eigenvalue is at most \(T_+\), and the \(2n\) contrast eigenvalues are strictly smaller, so \(Q<NT_+\).

No point in the inactive region can minimize the scale profile. There \(s\ge\log T_+\), and the diagonal-model objective has derivative
\[
p-e^{-s}\operatorname{tr}B_0>0,
\]
where \(B_0\) is the matrix before multiplication by \(e^{-s}\). Indeed \(\operatorname{tr}B_0<pT_+\). Decreasing \(s\) improves the objective. Thus (6.4) is the global scale optimum over all real \(s\), including both active and inactive regions.

The remaining discrepancy is
\[
\mathcal H(d)=N\log\frac{Q(d)}N-2n\log a-\log T_-(d).
\tag{6.5}
\]

### Lemma 6.2 — An exact monotonicity formula. **Proved.**

Put \(\Delta=b^2-h^2>0\), \(k=2na/b>0\), and
\[
w=\frac{\sqrt{b^2\sinh^2d+h^2}}{b\cosh d}
=\sqrt{1-\frac{1-q^2}{\cosh^2d}}.
\tag{6.6}
\]
As \(|d|\) increases from zero to infinity, \(w\) increases strictly from \(q\) to one. Moreover
\[
b\cosh d=\frac{\sqrt\Delta}{\sqrt{1-w^2}},\quad
T_-=b\cosh d(1-w),\quad
Q=b\cosh d(k+1-w).
\]
Substitution in (6.5) gives
\[
\mathcal H(w)=n\log\Delta-2n\log a-N\log N
+N\log(k+1-w)-m\log(1-w)-n\log(1+w).
\tag{6.7}
\]
Differentiating and collecting terms yields
\[
\frac{d\mathcal H}{dw}
=\frac{(k-2n)+(Nk+2n)w}{(k+1-w)(1-w^2)}
=\boxed{\frac{(Nk+2n)(w-q_*)}{(k+1-w)(1-w^2)}}.
\tag{6.8}
\]
Here
\[
\frac{2n-k}{Nk+2n}=\frac{b-a}{b+Na}=\frac{L}{2-L}=q_*.
\]
The denominator is positive. The profile decreases strictly up to \(q_*\), increases strictly afterwards, and diverges to infinity as \(w\uparrow1\). Its unique minimizing value is therefore
\[
\boxed{w_{\min}=\max(q,q_*).}
\tag{6.9}
\]
This proves the global comparison and excludes extra local minima in the reduced problem. It does not rely on curvature alone.

### Lemma 6.3 — Symmetric and asymmetric solutions. **Proved.**

The condition \(q\ge q_*\) is exactly \(\phi\ge\phi_*\). In this case (6.9) requires \(d=0\). Equations (6.2) and (6.4) give
\[
\psi_A=\psi_B=\frac{2na+b-h}N=a+\frac{mL(1-\phi)}N=\psi_s.
\]
The loading obtained from the leading eigenvector is uniform and has square \(1-\psi_s\). Thus (1.5) follows. Substitution in (6.5) gives the symmetric expression in (1.4).

For \(0\le q<q_*\), the two choices of \(d\) satisfy
\[
\cosh^2d=\frac{1-q^2}{1-q_*^2}
=\frac{(2-L)^2(1-q^2)}{4a}.
\tag{6.10}
\]
Using \(k+1-q_*=2aN/[b(2-L)]\) in (6.4) gives
\[
\boxed{\psi_A+\psi_B=(2-L)(1-q^2),\qquad \psi_A\psi_B=a(1-q^2).}
\tag{6.11}
\]
At \(q=0\), the roots are \(\{a,1\}\). For \(0<q<q_*\), they are two distinct numbers in \((0,1)\). For completeness, their sum and product are positive and their discriminant is
\[
(1-q^2)\{L^2-(2-L)^2q^2\}>0.
\]
Also \((1-\psi_A)(1-\psi_B)=q^2>0\) and \(\psi_A+\psi_B<2\); hence both are below one. At \(q=q_*\) they coincide.

Define \(\alpha=\sqrt{1-\psi_A}\), \(\beta=\sqrt{1-\psi_B}\), with the order chosen so \(\alpha\ge\beta\). Equation (6.11) gives exactly (1.6) and (1.7). These loadings really are the concentrated optimum, rather than merely reproducing the diagonal: from (6.11),
\[
1/\psi_A+1/\psi_B=(2-L)/a,\qquad
\gamma=m(1/\psi_A+1/\psi_B-2)=mL/a,
\]
and, using \(h=bq=b\alpha\beta\),
\[
b\alpha/\psi_A+h\beta/\psi_B
=b\alpha(1/\psi_A+1/\psi_B-1)=b\alpha/a.
\tag{6.12}
\]
The analogous equation holds for \(\beta\). The two mean-space eigenvalues are therefore \(b/a=1+\gamma\) and \(b\). The former is the leading eigenvalue and the whitened loading norm is \(\gamma\), exactly as required by Lemma 3.1. At \(q=0\) the leading eigenvalue at these optimal uniquenesses remains simple.

The nonleading eigenvalues sum to
\[
na(1/\psi_A+1/\psi_B)+b=n(2-L)+b=N,
\]
and their product is
\[
\left(\frac{a^2}{\psi_A\psi_B}\right)^n b
=\left(\frac{a}{1-q^2}\right)^n b.
\]
Equation (3.2) thus gives
\[
\boxed{F_a=-n\log a-\log b+n\log(1-q^2).}
\tag{6.13}
\]

### Completion of the global and local statements for \(\phi>0\). **Proved.**

Existence is established in Lemma 2.2. Every local minimum is positive, interior, and block-uniform by Sections 3–5. At such a point its loading is the unique fixed-\(\Psi\) minimizer up to sign. This concentrated loading varies continuously with nearby block-uniform uniquenesses. Indeed, the leading eigenvalue is continuous by the Rayleigh variational formula; any subsequential limit of normalized positive leading eigenvectors is the unique normalized positive leading eigenvector at the limiting matrix, so the eigenvectors converge as well. Their optimal loading norms also converge. A local minimum of the full problem must therefore be a local minimum of the profile in Section 6. Formula (6.8) and the unique scale optimum leave exactly the minima in Theorem 1.

Conversely, at least one of these candidates is a full global minimum by existence and reduction. The only candidates below threshold are block interchanges with equal objective, so both are global. Above threshold, including equality, there is only the symmetric candidate, hence it is global and unique as a covariance. Every listed covariance has nonzero off-diagonal products when \(\phi>0\), which identify its loadings up to global sign: for distinct \(i,j,k\), \(l_i^2=\Sigma_{ij}\Sigma_{ik}/\Sigma_{jk}\). This proves parameter uniqueness too.

There are consequently no asymmetric local minima above threshold, and no interval of coexistence with a globally optimal symmetric fit.

## 7. The symmetric saddle and the transition

### Lemma 7.1 — A smooth descending direction below threshold, including \(\phi=0\). **Proved.**

The symmetric point (1.5) satisfies both equations of Lemma 3.2, so it is stationary for every \(0\le\phi<1\). For \(\phi>0\), (6.8) already proves negative symmetry-breaking curvature below threshold. At \(\phi=0\), however, the concentrated eigenvalue has a degeneracy at equal uniquenesses; a cusp of a concentrated objective alone would not establish a descent path from the specified symmetric loading. The following direct calculation treats all cases uniformly.

Let \(z=l_s^2\), \(r_s=1-z\), \(U=b+h=1+Nz\), and \(V=b-h=Nr_s-2na\). Perturb
\[
\alpha=\sqrt z+\epsilon,\quad\beta=\sqrt z-\epsilon,
\qquad \psi_A=1-\alpha^2,\quad \psi_B=1-\beta^2.
\tag{7.1}
\]
For sufficiently small \(\epsilon\), both uniquenesses remain positive. On normalized common and difference block-mean vectors, the population covariance is \(\operatorname{diag}(U,V)\), and the model covariance is
\[
\begin{pmatrix}
U-\epsilon^2&2n\sqrt z\,\epsilon\\
2n\sqrt z\,\epsilon&r_s+N\epsilon^2
\end{pmatrix}.
\tag{7.2}
\]
The contrast covariances are \(\psi_A\) and \(\psi_B\), with population contrast variance \(a\). Expanding their trace and log-determinant contributions gives
\[
F(\epsilon)-F(0)=E\epsilon^2+O(\epsilon^4),
\]
\[
E=
\frac{2n\{a(1+3z)-r_s(1+z)\}}{r_s^3}
+\frac{N(r_s-V)}{r_s^2}
+\frac{4n^2zV}{Ur_s^2}
=\boxed{
\frac{4nm(1+z)\{(2-L)z-L\}}
{(1-z)^3(1+Nz)}}.
\tag{7.3}
\]
Here the remainder is even by block interchange. One can check the simplification by multiplying the first expression by \(r_s^3U/(4n)\): after using \(V=Nr_s-2na\), its numerator is
\[
a\{[m+(1-n)z]U-2n^2zr_s\}
+r_s\{-[m-nz]U+nNzr_s\}
=m(1+z)[a(1+z)-r_s].
\]
Finally,
\[
(2-L)z-L=\frac{mL(2-L)}N(\phi-\phi_*).
\tag{7.4}
\]
Thus \(E<0\) exactly below threshold, including \(\phi=0\).

There is also a strictly ascending direction. A common covariance scaling \(\Sigma\mapsto t\Sigma\) is feasible, with \(l\mapsto\sqrt t\,l\), \(\Psi\mapsto t\Psi\). Stationarity implies \(\operatorname{tr}(\Sigma^{-1}R)=p\), and
\[
F(t\Sigma)-F(\Sigma)=p(t^{-1}-1+\log t)
=\tfrac p2(t-1)^2+O((t-1)^3)>0
\]
for small \(t\ne1\). The symmetric point is therefore a saddle. \(\square\)

**Symbolic verification of (7.3).** The coefficient \(E\) was checked exactly in SymPy (script `check_saddle_7_3_sympy.py`), for symbolic \(m\), \(L\) and \(\phi\), independently of the expansion above. The check uses the exact block decomposition of \(F\) along the path (7.1): \(2n\) within-block contrast directions, on which \(R\) acts as \(a\) and the model as \(\psi_A\) or \(\psi_B\), plus the two-dimensional block-mean subspace, on which \(R=\begin{pmatrix}b&h\\h&b\end{pmatrix}\) and the model is \(\operatorname{diag}(\psi_A,\psi_B)+m(\alpha,\beta)^T(\alpha,\beta)\). It writes \(L=zN/(n+m\phi)\), so every quantity is rational in \(\sqrt z\). SymPy returns:
- first derivative at \(\epsilon=0\): \(0\) (stationarity);
- third derivative at \(\epsilon=0\): \(0\) (block-interchange symmetry, so the remainder is \(O(\epsilon^4)\));
- \(\tfrac12F''(0)\) minus the closed form in (7.3): \(0\); minus the three-term form in (7.3): \(0\);
- the identity (7.4): \(0\);
- at \(\phi=0\): \((2-L)z-L=-L\,(1+(m-1)L)/(2m-1)<0\).

**A shorter argument at \(\phi=0\), avoiding (7.3).** At \(\phi=0\) the block-mean part of \(R\) is \(b\,I_2\), so for \(\Psi=\psi_sI\) the leading eigenvalue of \(B=R/\psi_s\) is \(b/\psi_s\), with the two-dimensional eigenspace spanned by \(e_A=1_A/\sqrt m\) and \(e_B=1_B/\sqrt m\). With \(c^2=2mz\), consider the circle of loadings
\[
l(t)=c\,(\cos t\,e_A+\sin t\,e_B),\qquad 0\le t\le\pi/2,\qquad \Psi=\psi_sI,
\]
which passes through the symmetric point at \(t=\pi/4\). Because \(\Sigma=\psi_sI+ll^T\), we have \(\log\det\Sigma=p\log\psi_s+\log(1+c^2/\psi_s)\) and \(\operatorname{tr}(\Sigma^{-1}R)=\operatorname{tr}R/\psi_s-l^TRl/(\psi_s(\psi_s+c^2))\), with \(l^TRl=c^2(b+h\sin2t)=c^2b\) at \(\phi=0\). So \(F\) is constant on the circle. Every point of the circle has the optimal loading for \(\Psi=\psi_sI\): a leading eigendirection with the same norm as at the stationary symmetric point. Hence \(\partial F/\partial l=0\) there, and by the computation in Lemma 3.2 the uniqueness gradient is \(\partial F/\partial\psi_i=(\Sigma_{ii}-1)/\psi_s^2\). For \(t\ne\pi/4\) the diagonal entries are \(\psi_s+c^2\cos^2t/m\) on block A and \(\psi_s+c^2\sin^2t/m\) on block B; they cannot both equal one, so the uniqueness gradient is nonzero. A small move of \(\Psi\) then lowers \(F\) strictly below the symmetric value. Points of the circle come arbitrarily close to the symmetric point, so it is not a local minimum. Together with the ascending scaling direction, it is a saddle. (Checked numerically: \(F\) is constant on the circle to \(10^{-15}\), and the uniqueness gradient matches \((\Sigma_{ii}-1)/\psi_s^2\) to \(10^{-15}\).)


### Lemma 7.2 — Supercritical square-root scaling and regularity. **Proved.**

Equation (1.6) gives
\[
(\alpha-\beta)^2=L+(2-L)q^2-2q
=(2-L)(q_*-q)(1-q).
\]
Since \(q=(mL/b)\phi\), this proves (1.9) and (1.10). At the threshold,
\[
\alpha=\beta=\sqrt{q_*},\qquad
\psi_A=\psi_B=1-q_*=\frac{2a}{2-L},
\]
which agrees with the symmetric solution. There is no jump.

For a Landau expansion, hold \(m,L\) fixed and use the profiled coordinate \(d\). At \(q=q_*\),
\[
w(d)=q_*+\frac{1-q_*^2}{2q_*}d^2+O(d^4),
\]
and (6.8) gives
\[
\mathcal H_{ww}(q_*)
=\frac{nm(2-L)^4}{4a^2N}>0.
\tag{7.5}
\]
Consequently
\[
\mathcal H(d;\phi)-\mathcal H(0;\phi)
=K_2(\phi-\phi_*)d^2+K_4d^4
+O\!\left((\phi-\phi_*)^2d^2+|\phi-\phi_*|d^4+d^6\right),
\]
\[
\boxed{K_2=\frac{nm^2(2-L)^3}{2aNb}>0,\qquad
K_4=\frac{nm(2-L)^2}{2NL^2}>0.}
\tag{7.6}
\]
This is a supercritical pitchfork in the profiled coordinate. Its global character follows from (6.8), not just this expansion.

Let the asymmetric formula in (1.4) be analytically extended locally across the threshold. In (6.7), the part depending only on \(q\) is common to the symmetric and asymmetric formulas; their difference is the \(w\)-dependent part evaluated at \(w=q\) minus its value at \(q_*\). Thus both the values and first derivatives agree at the threshold. Their second derivatives with respect to \(\phi\) differ by
\[
\boxed{
F_s''(\phi_*)-F_a''(\phi_*)
=\frac{n m^3L^2(2-L)^4}{4a^2Nb^2}>0.}
\tag{7.7}
\]
This proves continuity, real-analyticity away from \(\phi_*\), and exactly \(C^1\), not \(C^2\), regularity there. The formula near zero has a real-analytic extension across zero; analyticity at this endpoint is understood in that sense. \(\square\)

## 8. The zero-correlation endpoint and the full flat family

### Lemma 8 — Endpoint global minimum and equality cases. **Proved.**

At \(\phi=0\), write \(R=\operatorname{diag}(R_A,R_B)\), where the identical block correlation matrices have determinant \(a^n b\). First assume all uniquenesses are positive. Then \(B=\Psi^{-1/2}R\Psi^{-1/2}\) is block diagonal.

If its largest eigenvalue exceeds one and is selected from block \(A\), the concentrated discrepancy splits as
\[
F_{\mathrm{conc}}
=F_{A,\text{one-factor}}(\Psi_A)
+F_{B,\text{diagonal}}(\Psi_B).
\tag{8.1}
\]
This formula remains valid if the leading eigenvalue is tied across blocks: removing one copy of that eigenvalue gives the same sum. The first term is nonnegative. The second is
\[
\sum_{i\in B}(1/\psi_i+\log\psi_i-1)-\log\det R_B
\ge-\log\det R_B=F_0,
\tag{8.2}
\]
with equality only if \(\Psi_B=I\). The same argument applies with \(A,B\) interchanged. If all eigenvalues are at most one, concentration gives a full diagonal model and a lower bound \(2F_0>F_0\). To check strictness, view \(F_0\) as a function of \(L\):
\[
F_0\big|_{L=0}=0,\qquad
\frac{dF_0}{dL}=n\left(\frac1{1-L}-\frac1{1+nL}\right)>0
\quad(L>0).
\]

The bound \(F_0\) is attained by fitting one block exactly and setting the other to the identity. For equality in (8.1), the selected block must be fitted exactly and \(\Psi_B=I\). In an exact active block, all loadings are nonzero, and with positive uniquenesses every \(0<\psi_i<1\). The entries of its whitened matrix are then strictly larger than those of \(R_A\). Evaluation on the positive uniform leading eigenvector of \(R_A\) shows its largest eigenvalue is strictly greater than \(b\), the largest eigenvalue in the unretained block. Thus at equality the leading eigenspace cannot mix blocks. The fitted covariance must be
\[
\operatorname{diag}(R_A,I_m)
\quad\text{or}\quad
\operatorname{diag}(I_m,R_B).
\tag{8.3}
\]

For \(m\ge3\), equations \(l_i l_j=L\) for every distinct retained pair imply all retained loadings have the same sign and magnitude \(\sqrt L=\lambda\). For \(m=2\), they require only \(l_1l_2=L\). Nonnegative uniquenesses and unit diagonal give \(|l_i|\le1\), hence exactly the family (1.8), with global sign included.

Finally, no missing boundary minimum invalidates the lower bound: Lemma 4 excludes all boundary local minima at this endpoint when \(m\ge3\). When \(m=2\), any boundary global minimum must be the unique star minimum on its stratum, which is already an endpoint of (1.8) and has value \(F_0\). This proves both the minimal value and all equality cases. \(\square\)

### Corollary 8.1 — Population CFI. **Proved.**

Here population CFI means the discrepancy ratio
\[
\mathrm{CFI}_{\mathrm{pop}}=1-\frac{F_{\min}}{F_{\mathrm{ind}}},
\]
where the independence-model optimum is \(I\), so \(F_{\mathrm{ind}}=-\log\det R\). By (2.1),
\[
F_{\mathrm{ind}}=2F_0-\log(1-q^2).
\]
In particular,
\[
\boxed{F_{\min}(0)=F_0,\quad F_{\mathrm{ind}}(0)=2F_0,
\quad\mathrm{CFI}_{\mathrm{pop}}(0)=\tfrac12.}
\tag{8.4}
\]
This is the population discrepancy-ratio statement, not a finite-sample claim about degrees-of-freedom corrections.

## 9. Equivalent threshold inequality

### Lemma 9 — Threshold forms and parameter range. **Proved.**

Let \(c=L/a\), \(\kappa_+=mc(1+\phi)\), and \(\kappa_-=mc(1-\phi)\). Since \(2mc>0\),
\[
N(\kappa_+-\kappa_-)\ge\kappa_-(2+\kappa_++\kappa_-)
\]
is equivalent, successively, to
\[
N\phi\ge(1-\phi)(1+mc),
\]
\[
\phi\ge\frac{1+mc}{N+1+mc}
=\frac{a+mL}{m(2a+L)}
=\frac{b}{m(2-L)}=\phi_*.
\tag{9.1}
\]
Also
\[
\phi_* =\frac{L+a/m}{2-L},\qquad
m(2-L)-b=Na>0.
\]
Thus \(0<\phi_*<1\), it decreases strictly with integer \(m\), and its limit is \(L/(2-L)\).

## 10. Reflection and equal cross-loadings

### Corollary C1 — Reflection. **Proved.**

For \(-1<\phi<0\), put \(D=\operatorname{diag}(I_m,-I_m)\). Then
\[
DR(\lambda,\phi,m)D=R(\lambda,-\phi,m),
\]
and
\[
D\Sigma D=(Dl)(Dl)^T+\Psi\in\mathcal M.
\]
The products \((D\Sigma D)^{-1}(DRD)\) and \(\Sigma^{-1}R\) are conjugate by \(D\), so they have the same trace and determinant. Hence the discrepancy is invariant and the theorem holds with \(|\phi|\). In the above-threshold regime the two blocks have equal-magnitude loadings of opposite signs. Below threshold, their magnitudes are unequal and their relative sign is the sign of \(\phi\).

### Corollary C2 — Equal cross-loadings, with explicit normalization (Proposition 3 of the article). **Proved with the normalization stated below.**

Suppose latent traits have variance one and covariance \(\rho\in[-1,1]\), so that \(\rho\) is the generating trait correlation. Before standardization, let an item in \(A\) have primary and secondary loadings \((\lambda_p,\lambda_s)\), and an item in \(B\) have loadings \((\lambda_s,\lambda_p)\). Assume residuals are mutually uncorrelated, uncorrelated with the traits, and have common variance \(\tau^2>0\). Define the signal variance \(\sigma_S^2\) and the between-block signal covariance \(\sigma_{AB}\) by
\[
\sigma_S^2=\lambda_p^2+\lambda_s^2+2\lambda_p\lambda_s\rho,\qquad
\sigma_{AB}=2\lambda_p\lambda_s+\rho(\lambda_p^2+\lambda_s^2).
\]
Every raw item has variance \(\tau^2+\sigma_S^2\); distinct within-block items have covariance \(\sigma_S^2\), and between-block items covariance \(\sigma_{AB}\). For \(\sigma_S^2>0\), standardization therefore produces exactly the matrix in this theorem with
\[
\boxed{\lambda_{\mathrm{eff}}^2=\frac{\sigma_S^2}{\tau^2+\sigma_S^2},
\qquad\phi_{\mathrm{eff}}=\frac{\sigma_{AB}}{\sigma_S^2}.}
\tag{10.1}
\]
The formula \(\sigma_S^2/(1+\sigma_S^2)\) is correct when **raw residual variance is one**, followed by standardization. If \(\lambda_p,\lambda_s\) are instead loadings of items already having variance one, then \(\tau^2=1-\sigma_S^2\) and the correct formula is **\(\lambda_{\mathrm{eff}}^2=\sigma_S^2\)**. Unit item variance alone does not imply \(\sigma_S^2/(1+\sigma_S^2)\).

The identities
\[
\sigma_S^2+\sigma_{AB}=(1+\rho)(\lambda_p+\lambda_s)^2\ge0,\qquad
\sigma_S^2-\sigma_{AB}=(1-\rho)(\lambda_p-\lambda_s)^2\ge0
\]
show \(|\phi_{\mathrm{eff}}|\le1\). With \(\sigma_S^2>0\), \(\tau^2>0\), and \(|\phi_{\mathrm{eff}}|<1\), apply Theorem 1 and C1 using the threshold computed at \(\lambda_{\mathrm{eff}}\):

| Effective correlation | Minimizer |
|---|---|
| \(\phi_{\mathrm{eff}}\ge\phi_*(\lambda_{\mathrm{eff}},m)\) | Equal-magnitude, same-sign block loadings: composite regime |
| \(\lvert\phi_{\mathrm{eff}}\rvert<\phi_*(\lambda_{\mathrm{eff}},m)\) | Unequal-magnitude block loadings: tilted regime, with the zero-correlation endpoint exception |
| \(\phi_{\mathrm{eff}}\le-\phi_*(\lambda_{\mathrm{eff}},m)\) | Equal-magnitude, opposite-sign block loadings: bipolar regime |

If \(\phi_{\mathrm{eff}}=\pm1\), the standardized population is itself an exact one-factor covariance and the minimum is zero. If \(\sigma_S^2=0\), the same inequalities force \(\sigma_{AB}=0\), giving the identity matrix. These degenerate cases lie outside the strict assumptions of Theorem 1 and are handled directly.

## 11. ULS: local threshold and a full global theorem

### Theorem U — Global unweighted least-squares fit (Theorem 2 of the article). **Proved.**

Let
\[
G(l)=\sum_{i<j}(R_{ij}-l_i l_j)^2,\qquad l\in[-1,1]^{2m},
\quad C_m=\frac{m(m-1)}2.
\]
Then
\[
\boxed{
G_{\min}=
\begin{cases}
C_mL^2-m^2(m-1)L^2\phi^2,&0\le\phi\le1/(2m),\\[3pt]
\dfrac{m^2(m-1)}{2m-1}L^2(1-\phi)^2,
&1/(2m)\le\phi<1.
\end{cases}}
\tag{11.1}
\]
For \(0<\phi<1/(2m)\), every global minimizer is block-uniform and, up to sign and block interchange,
\[
\boxed{
\alpha^2=\frac L2\left(1+\sqrt{1-4m^2\phi^2}\right),\qquad
\beta^2=\frac L2\left(1-\sqrt{1-4m^2\phi^2}\right),}
\tag{11.2}
\]
with \(\alpha,\beta>0\). Above threshold, including equality, it is symmetric with square \(L(n+m\phi)/N\). At \(\phi=0\), the equality cases are exactly one block fitted and the other zero, including (1.8) for \(m=2\). All positive-correlation global minimizers lie strictly inside the loading bounds.

Online Resource 2 extends the unrestricted global reduction to all positive
within-block and cross-block weights (Proposition 1 of the article, the
weighted least-squares result in Section 5.1). It includes threshold equality and all
zero-correlation equality cases, including the continuum at the critical
weight ratio.

The symmetric point has a negative Hessian eigenvalue exactly when \(\phi<1/(2m)\). The ULS asymmetric branch obeys
\[
\boxed{|\alpha-\beta|
=\lambda\sqrt{2m}\sqrt{\frac1{2m}-\phi}.}
\tag{11.3}
\]

### Lemma U1 — Existence, signs, and a relaxation. **Proved.**

The original domain is compact, so its minimum exists. Replacing loadings by their absolute values cannot increase \(G\), since \((r+z)^2-(r-z)^2=4rz\ge0\) for \(r,z\ge0\). For \(\phi>0\), all nonzero minimizing loadings have the same sign. The zero vector is not a minimum: assigning two coordinates a small common positive loading decreases their pair residual while leaving the other pair products unchanged. A nonnegative, nonzero minimizer cannot have a zero coordinate either, because there
\[
\partial_iG=-2\sum_{j\ne i}R_{ij}l_j<0.
\]
Thus it is strictly positive after choosing its global sign.

For nonnegative block loadings \(x_i,y_i\), write
\[
S_A=\sum_i x_i,\quad Q_A=\sum_i x_i^2,\quad
\xi=\frac{S_A^2-Q_A}{m(m-1)},
\]
and analogously \(S_B,Q_B,\eta\). The numbers \(\xi,\eta\) are mean within-block pair products. Set \(K_m=m(m-1)=2C_m\). Cauchy–Schwarz gives
\[
Q_A\ge m\xi,\quad Q_B\ge m\eta,\quad
S_A^2=Q_A+K_m\xi,\quad S_B^2=Q_B+K_m\eta.
\]
Convexity of the squared residual within each block, and direct expansion of the cross-block sum, yield
\[
G\ge f(\xi,\eta,Q_A,Q_B),
\]
\[
f=C_m[(L-\xi)^2+(L-\eta)^2]+Q_AQ_B
-2L\phi\sqrt{(Q_A+K_m\xi)(Q_B+K_m\eta)}+m^2L^2\phi^2
\tag{11.4}
\]
on the enlarged closed domain \(\xi,\eta\ge0\), \(Q_A\ge m\xi\), \(Q_B\ge m\eta\). Its upper bounds have deliberately been removed. If both inequalities are equalities, block-uniform loadings \(x_i=\sqrt\xi,y_i=\sqrt\eta\) attain this bound.

### Lemma U2 — The relaxed minimum exists for \(\phi>0\). **Proved.**

The bound \(Q_A+K_m\xi\le mQ_A\), and its analogue in \(B\), imply
\[
f\ge C_m[(L-\xi)^2+(L-\eta)^2]
+(\sqrt{Q_AQ_B}-mL\phi)^2.
\tag{11.5}
\]
On a bounded sublevel, \(\xi,\eta,Q_AQ_B\) are bounded. An unbounded sequence in that sublevel has a subsequence, after block interchange if necessary, with \(Q_A\to\infty\), \(Q_B\to0\), and \(\eta\to0\). Put \(T=Q_AQ_B\). Since \(\xi\) is bounded,
\[
\sqrt{(Q_A+K_m\xi)(Q_B+K_m\eta)}
\le\sqrt{mT(1+K_m\xi/Q_A)}.
\]
Taking a further subsequence for bounded quantities and minimizing the limiting quadratic in \(\sqrt T\) gives
\[
\liminf f\ge C_mL^2+
\inf_{T\ge0}\{T-2L\phi\sqrt{mT}+m^2L^2\phi^2\}
=E,
\]
\[
E=C_mL^2+m(m-1)L^2\phi^2.
\tag{11.6}
\]
There is a finite feasible point below \(E\). If \(0<\phi\le1/(2m)\), choose uniform-block loadings with \(\alpha^2+\beta^2=L\), \(\alpha\beta=mL\phi\); their value is the first expression in (11.1), and the gap below \(E\) is \(m(m-1)(m+1)L^2\phi^2>0\). If \(\phi\ge1/(2m)\), the symmetric candidate has the second value in (11.1), and
\[
E-G_s=\frac{m(m-1)L^2}{2m-1}
\left[-\tfrac12+2m\phi+(m-1)\phi^2\right]>0.
\]
These candidate values follow by direct substitution in \(G\); their global optimality is not being assumed. Thus no minimizing sequence can escape to infinity. Compactness of a bounded minimizing sequence proves that the relaxed minimum is attained.

### Lemma U3 — Both relaxed norm inequalities must be equalities. **Proved.**

On the relaxed domain, write \(S_A=\sqrt{Q_A+K_m\xi}\) and \(S_B=\sqrt{Q_B+K_m\eta}\).

At a relaxed minimizer for \(\phi>0\), \(Q_A,Q_B>0\). If, for example, \(Q_A=0\), then \(\xi=0\); introducing \((\xi,Q_A)=(t^2,mt^2)\) when \(Q_B>0\) changes \(f\) by
\[
-2mL\phi\,t\sqrt{Q_B+K_m\eta}+O(t^2)<0
\]
for small \(t>0\). If both norms are zero, small positive loadings in both blocks improve the value. Also \(\xi,\eta>0\): at \(\xi=0,Q_A>0\), increasing \(\xi\) slightly is feasible and
\[
\frac1{K_m}\frac{\partial f}{\partial\xi}
=\xi-L-L\phi\frac{S_B}{S_A}<0.
\]

Suppose \(Q_A>m\xi\). The variables \(Q_A,\xi\) are locally free, and their stationarity equations give
\[
0=Q_B-L\phi\frac{S_B}{S_A},\qquad
0=\xi-L-L\phi\frac{S_B}{S_A}.
\]
Hence \(\xi=L+Q_B\) and \(Q_A>mQ_B\). Both blocks cannot have strict norm inequalities, since the same argument would give \(Q_B>mQ_A\). Therefore \(Q_B=m\eta\).

Variation in \(\eta\) along that boundary is two-sided and yields
\[
0=K_m(\eta-L)+mQ_A-m^2L\phi\frac{S_A}{S_B}.
\]
Since \(S_B^2=m^2\eta\) and \(L\phi S_B/S_A=m\eta\), we have \(L\phi S_A/S_B=S_A^2/m\). The last equation reduces to \(\eta=L+m\xi\). Combined with \(\xi=L+m\eta\), it gives
\[
(1-m^2)\eta=(m+1)L,
\]
which is impossible. Thus \(Q_A=m\xi\), and symmetry gives \(Q_B=m\eta\).

The relaxed minimizer is consequently attained by uniform loadings in each block. Once the resulting loads are shown below to satisfy the original bounds, the relaxed minimum equals the original minimum. Any original minimizer must then also minimize the relaxation. Equality in \(Q_A\ge m\xi\), or \(S_A^2\le mQ_A\), forces all its block-A loadings equal; likewise for \(B\). This proves the unrestricted global reduction, including the loading-boundary cases.

### Lemma U4 — Exact cone minimization and endpoint equality. **Proved.**

For block-uniform loadings,
\[
G=C_m[(L-\alpha^2)^2+(L-\beta^2)^2]
+m^2(L\phi-\alpha\beta)^2.
\]
Set \(u=\alpha^2+\beta^2\), \(v=\alpha\beta\). Nonnegative loadings map onto the cone \(u\ge2v\ge0\), and
\[
\boxed{
G=C_m(u-L)^2+m(v-mL\phi)^2
+C_mL^2-m^2(m-1)L^2\phi^2.}
\tag{11.7}
\]
This is a strictly convex quadratic in \((u,v)\). Its unconstrained minimizer is \((L,mL\phi)\), feasible exactly for \(\phi\le1/(2m)\). For \(\phi>1/(2m)\) the unconstrained minimizer is infeasible, and the constrained minimizer is identified by the Karush-Kuhn-Tucker conditions. These are sufficient here, because the objective is convex and the constraints \(u-2v\ge0\), \(v\ge0\) are linear. Write \(Q(u,v)=C_m(u-L)^2+m(v-mL\phi)^2\) and \(v_s=L(n+m\phi)/N>0\), and take the candidate on the edge \(u=2v\) at \((u,v)=(2v_s,v_s)\). The constraint \(v\ge0\) is inactive there. Stationarity of \(Q+\varpi(2v-u)\) requires
\[
\partial_uQ-\varpi=0,\qquad \partial_vQ+2\varpi=0 .
\]
The first gives \(\varpi=2C_m(2v_s-L)\). Substituting it into the second gives \(2m(v-mL\phi)+4C_m(2v-L)=0\), which holds exactly at \(v=v_s\). Hence
\[
\varpi=2C_m(2v_s-L)=\frac{m(m-1)L\,(2m\phi-1)}{2m-1}\ \ge\ 0\quad\Longleftrightarrow\quad\phi\ge\frac1{2m}.
\]
So for \(\phi\ge1/(2m)\) the edge point is the global minimizer over the cone, and it is unique by strict convexity. At \(\phi=1/(2m)\) the multiplier vanishes and the edge point coincides with the unconstrained minimizer \((L,L/2)\). For \(\phi<1/(2m)\) the multiplier would be negative, consistent with the unconstrained minimizer being feasible. (The multiplier formula was checked symbolically in SymPy.) This proves (11.1), (11.2), and uniqueness up to block interchange. Every resulting squared loading is at most \(L<1\), validating the relaxation argument for the bounded original problem.

For \(\phi=0\), use (11.4) directly:
\[
G\ge C_m[(L-\xi)^2+(L-\eta)^2]+Q_AQ_B
\ge C_m[(L-\xi)^2+(L-\eta)^2]+m^2\xi\eta.
\]
Equation (11.7) at zero correlation shows this is at least \(C_mL^2\), with equality only when \((\xi,\eta)=(L,0)\) or its interchange. In the first case \(Q_A>0\), so equality requires \(Q_B=0\). Within-block equality requires every retained product \(x_ix_j=L\). For \(m\ge3\), all retained loads equal \(\lambda\). For \(m=2\), this gives exactly \((t,L/t)\), \(L\le t\le1\), up to simultaneous sign. Thus all endpoint equality cases have been accounted for.

### Lemma U5 — Complete symmetric Hessian and square-root law. **Proved.**

At \(l_i=\sqrt{q_s}\), where \(q_s=L(n+m\phi)/N\), the Hessian has entries
\[
H_{ii}=2Nq_s,\qquad H_{ij}=4q_s-2R_{ij}\quad(i\ne j).
\]
The block invariant subspaces give the following eigenvalues directly:

| Subspace | Multiplicity | Hessian eigenvalue |
|---|---:|---|
| Within-block sum-zero contrasts | \(2n\) | \(2[(2m-3)q_s+L]\) |
| Common direction \((1_A,1_B)\) | 1 | \(4L(n+m\phi)\) |
| Block difference \((1_A,-1_B)\) | 1 | \(4Ln(2m\phi-1)/N\) |

The first two are strictly positive. The last is negative exactly below \(1/(2m)\), zero at that threshold, and positive above it. The symmetric point is interior to the loading bounds, so this is the full local stability classification. On the asymmetric branch, \(\alpha^2+\beta^2=L\) and \(\alpha\beta=mL\phi\); hence
\[
(\alpha-\beta)^2=L-2mL\phi,
\]
which proves (11.3).

The global ULS transition is continuous and has no first-order change of global minimizer. The preceding proof does not classify every nonglobal ULS local minimum in the full constrained domain.

## 12. Independent numerical verification (not part of the proof)

The accompanying `verify_symmetry_breaking.py` evaluates the original covariance discrepancy, concentrates independently through numerical eigendecomposition, runs optimization over all loadings and positive uniquenesses from random starts, checks boundary derivatives, compares stationarity and spectrum identities, tests the endpoint family, and checks reflection, cross-loading normalization, and ULS formulas. Its floating-point results do not enter any proof above.

Selected ML values, evaluated from the exact solution and independently comparable against the original objective, are:

| \(m\) | \(\lambda\) | \(\phi\) | Larger block loading | Smaller block loading | \(F_{\min}\) | Symmetric \(F\) | Population CFI |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 0.7 | 0 | 0.700000 | 0 | 0.274568433 | 0.412924380 | 0.500000 |
| 2 | 0.7 | 0.2 | 0.692882 | 0.189850 | 0.257113249 | 0.310046931 | 0.546211 |
| 2 | 0.7 | 0.473, literally | 0.619918 | 0.501842 | 0.172775294 | 0.173070665 | 0.734572 |
| 2 | 0.7 | \(\phi_*-0.02=0.4733774834\ldots\) | 0.619521 | 0.502564 | 0.172604179 | 0.172888775 | 0.734904 |
| 2 | 0.7 | 0.513 | 0.575251 | 0.575251 | 0.153985914 | 0.153985914 | 0.770170 |
| 4 | 0.8 | 0 | 0.800000 | 0 | 1.993370126 | 2.814131657 | 0.500000 |
| 4 | 0.8 | 0.3 | 0.789336 | 0.333209 | 1.778314273 | 1.971951635 | 0.561822 |
| 4 | 0.8 | 0.6 | 0.702648 | 0.702648 | 1.040350814 | 1.040350814 | 0.758657 |

The thresholds are \(0.4933774834437086\ldots\) for \((m,\lambda)=(2,0.7)\), and \(0.536764705882353\ldots\) for \((4,0.8)\).

**Rounding note.** The reference values \(0.172604\) and \(0.172889\) used by the script are reproduced at \(\phi=\phi_*-0.02\), not at literal \(\phi=0.473\). The script includes both inputs explicitly. The other reference values agree with the formulas to the stated rounding.

An earlier 144-population, 25-start search was reported without its exact eight correlation values. The script's optional `--full-grid` therefore uses a documented independent 144-population grid; it does not claim to reproduce that earlier grid exactly. Finite optimizer uniqueness bounds and floating-point tolerances are recorded by the script. Boundary points are checked separately through the exact formulas and direct covariance evaluation.

### Executed verification

The following command was run successfully, using Python 3.12 with NumPy and SciPy; exact versions and every recorded numerical result are in `verification_results.json`:

```bash
python verify_symmetry_breaking.py \
  --full-grid --starts 25 --jobs 4 --seed 20261008 \
  --output verification_results.json
```

**Result: 4,224 checks passed, zero failures.** This run included 128 random identity-check populations, the full independent grid of 144 ML populations with 25 random starts each (3,600 unconstrained-loading, untied-uniqueness optimizations), and 12 ULS populations with eight random starts each. No analytic starting solution was included in the ML optimization pool. The grid was
\[
m\in\{2,3,4\},\quad
\lambda\in\{.4,.5,.6,.7,.8,.9\},\quad
\phi\in\{0,.05,.15,.30,.45,.60,.80,.95\}.
\]

For the best numerical ML fit in each population, the largest absolute discrepancy error relative to the exact formula was \(5.3291\times10^{-15}\); the largest fitted-covariance entry error, after allowing block interchange, was \(6.7805\times10^{-8}\). Numerical parameter errors are larger than objective errors, as expected near flat directions. The largest ULS objective error in this run was \(5.8565\times10^{-15}\).

The ML optimizer parameterizes all \(2p\) variables as \((l,\log\Psi)\), with numerical bounds \(\log\psi_i\in[-20,5]\). These bounds keep its search inside the positive-uniqueness region; they do not constitute a search over the exact Heywood boundary. The boundary derivatives and exceptional endpoint covariances are evaluated separately. Away from the known two-item endpoint exception, the reported best fits did not approach the numerical uniqueness bounds. These are implementation checks, not mathematical optimality certificates.

## 13. Gaps, assumptions, and exceptional cases

1. **The main ML theorem has no remaining proof gap.** Existence, unrestricted reduction, global optimality, uniqueness, all positive-correlation local minima, the endpoint, and the supercritical transition are proved for the full stated parameter range.
2. **The saddle algebra (7.3) is verified symbolically** (SymPy, symbolic \(m,L,\phi\)); a shorter argument at \(\phi=0\) that avoids it is given after Lemma 7.1.
3. **No factor-analysis result is used only by citation.** Fixed-uniqueness concentration, score equations, diagonal reproduction, and Heywood handling are derived in Sections 3–4. The Heywood inward-derivative formula (4.3) is derived in full in Lemma 4, Steps 2 to 4, and holds for any correlation matrix. No interval arithmetic or certified root-isolation fallback is needed.
4. **The two-item endpoint exception is essential.** At \(m=2,\phi=0\), there are two minimizing covariance matrices, with a continuous family of parameter representations for each. Each family's endpoints have a zero uniqueness. For \(\phi>0\), this nonidentification disappears. The positive-correlation branch approaches the uniform representative \(t=\lambda\) as \(\phi\downarrow0\).
5. **Local minima and parameter topology.** The no-coexistence proof applies to local minima in the original \((l,\Psi)\) parameterization, including its feasible boundary, and hence in particular excludes any competing minimizing covariance. No classification of all nonminimum stationary points is needed.
6. **Corollary C2 needs a stated measurement normalization.** The formula \(\lambda_{\mathrm{eff}}^2=\sigma_S^2/(1+\sigma_S^2)\) is valid for raw residual variance one. For already standardized loadings, \(\lambda_{\mathrm{eff}}^2=\sigma_S^2\). Residual independence or zero residual covariances, common residual variance, and no residual–trait covariances must be assumed. The degenerate cases \(\sigma_S^2=0\) and \(|\phi_{\mathrm{eff}}|=1\) are handled separately.
7. **CFI means the population discrepancy ratio.** No finite-sample correction or assertion about sample-estimated CFI is included.
8. **The ULS global theorem and symmetric local threshold are proved.** A complete classification of other, nonglobal ULS local minima is **Not proved** here. It is not needed for the full global value or the continuous global transition.
9. **Numerical verification is not proof.** Random starts, finite grids, finite differences, and floating-point residual checks support implementation correctness only. They are not certificates of global optimality; the analytic arguments establish that.

## Appendix to Online Resource 1. An independent algebraic classification of positive ML stationary points

**Status: Proved.** This appendix provides a second route to the stationary branch equations, independent of the one-dimensional profile.

By Lemma 5.2, a positive interior stationary point has block loadings \(\alpha,\beta\). Set \(u=\alpha^2\), \(v=\beta^2\), \(z=\alpha\beta\), and \(t=z-L\phi\). The row scores (3.6) are
\[
n\frac{\alpha(L-u)}{1-u}
=m\frac{\beta t}{1-v},\qquad
n\frac{\beta(L-v)}{1-v}
=m\frac{\alpha t}{1-u}.
\tag{A.1}
\]
If \(t\le0\), these equations force \(u,v\ge L\), which implies \(z\ge L>L\phi\), a contradiction. Thus \(t>0\), and \(u,v<L\).

Eliminating \(t\) from (A.1) gives
\[
u(L-u)(1-v)^2-v(L-v)(1-u)^2=0,
\]
whose factorization is
\[
\boxed{(u-v)\{L-u-v+(2-L)uv\}=0.}
\tag{A.2}
\]
If \(u=v\), (A.1) yields \(u=L(n+m\phi)/N\), the symmetric solution. Otherwise
\[
u+v=L+(2-L)z^2,
\quad (L-u)(L-v)=a^2z^2.
\]
Multiplying (A.1) gives \(n^2(L-u)(L-v)=m^2t^2\), so \(t=(na/m)z\). Therefore
\[
z-L\phi=\frac{na}{m}z
\quad\Longrightarrow\quad
z=\frac{mL\phi}{b}=q.
\]
The discriminant is nonnegative exactly when
\[
u+v-2z=L+(2-L)z^2-2z=(1-z)\{L-(2-L)z\}\ge0.
\]
Because \(0<z<1\), an asymmetric solution is possible exactly for \(z<q_*\), with equality giving the symmetric threshold. This reproduces (1.6)–(1.7) and independently rules out additional positive stationary branches. It is consistent with, and not needed to replace, the global monotonicity proof in Section 6.

## References

Horn, R. A., & Johnson, C. R. (2013). *Matrix analysis* (2nd ed.). Cambridge University Press. Chapter 8, positive and nonnegative matrices (Perron's theorem).
