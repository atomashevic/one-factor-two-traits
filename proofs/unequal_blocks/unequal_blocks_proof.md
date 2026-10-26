---
title: 'Online Resource 3 for “Fitting One Factor to Two Traits: Population Solutions and Fit Indices”'
subtitle: 'One factor fitted to two traits beyond exchangeability'
---

## Scope

This Online Resource treats two-trait populations whose blocks differ in size or in population loading. Section 1 proves Lemma 1 of the article (constant fitted loadings within blocks), including the exclusion of zero-uniqueness minima and the exact scope of the interchangeable-items extension. Section 2 treats the zero-correlation endpoint described in Section 5.2 of the article, reduces the positive stationary problem to an explicit quartic with an exact selection rule for arbitrary block sizes (Appendix D of the article), proves Proposition 2 for equal block sizes with unequal loadings, and analyses the unequal-size example of Section 5.2. Section 3 proves Theorem 3, the cusp unfolding near the exchangeable threshold, with explicit constants. Section 4 characterizes the local stability threshold for mirrored heterogeneous loadings, Section 5 reports the numerical verification, Section 6 lists the remaining gaps, and Section 7 collects summary statements of the main results. Each claim is labelled as proved, proved with computer algebra, or numerical only.

**Notation.** Symbols follow the article where they denote the same objects: \(m_A,m_B\) are the block sizes and \(\nu_{AB}=m_A+m_B-1\); \(\psi_A,\psi_B\) are block uniquenesses; \(\mathcal G\) denotes a profiled discrepancy; \(k_A,k_B,c_B,c_A,e_0\) are the coefficients of the stationary quartic; \(\mathrm{TC}(R_A)=-\tfrac12\log\det R_A\) is the total correlation of the within-block correlation matrix \(R_A\); and \(C_{\mathrm{tilt}}\), \(E_{\mathrm{fold}}\) are the constants of Theorem 3. For traceability with `verify_unequal_blocks.py`, the script writes \(m_A,m_B\) as `n`, `k` and \(\nu_{AB}\) as `q`; the coefficients \(k_A,k_B,c_B,c_A,e_0\) correspond to `A`, `B`, `C`, `D`, `E`; the block uniquenesses \(\psi_A,\psi_B\) are `x`, `y`; and the profile \(\mathcal G\) is `G` (for example, `profile_g`). In Section 3, which has equal block sizes, \(n=m-1\) and \(N=2m-1\) (the article's \(\nu\); the article uses \(N\) for sample size); in Section 2.3, \(N\) instead denotes an auxiliary function of \(v\), written \(H\) in Appendix D of the article. Symbols are otherwise local to their displayed definitions. The scalar coordinate \(v\) in Section 2 is a ratio of precision-weighted fitted loadings, not a trait weight. The normalized coordinate used only in the equal-size uniqueness proof differs from the original cusp coordinate by a fixed shift; that change is displayed explicitly.

## Summary of results

**Proved.** Lemma 1 of the article holds under the stated two-trait assumptions, for all block sizes at least two and all strictly interior population loadings. Every local minimum is interior, positive up to a common sign, and constant within each block. The proof handles both choices of the Heywood block, a two-item block, and every finite boundary representation. At zero trait correlation, the factor fits exactly the block with the larger Gaussian multi-information (total correlation). An equality of the two information values gives two minimizing covariances, with a flat parameter family only when a retained block has two items.

**Proved.** The positive stationary problem reduces to an explicit quartic, with at most three positive roots. Every positive root is admissible, and the within-block contrast Hessians are positive definite. Thus simple roots are full minima or one-index saddles according to the scalar profile. With equal block sizes and unequal loadings, the global minimum is unique and favors the stronger block for every positive trait correlation (Proposition 2 of the article). The weaker-block minimum exists below one fold precisely when the strict criterion (S11) holds. Both cusp exponents in Theorem 3 are proved, with closed-form coefficients and controlled remainders.

**Proved with computer algebra.** Exact SymPy identities verify the quartic expansion, its profile derivative, Hessian determinants, and cusp coefficients. An exact discriminant and Sturm certificate establish that the unequal-size example has a unique global branch throughout \([1/2,1]\), with a transverse equal-share crossing at \(0.757539241683999\ldots\). A computer-assisted interval certificate (Section 2.6) extends the global comparison to \(0<\phi\le.51\), so the global minimum in the example is unique for every \(0<\phi<1\).

**Counterexample and qualification.** Interchangeability alone does not imply the lemma: a positive definite, interchangeable four-item correlation matrix has a strict Heywood local minimum. A valid common-parts extension assumes the sufficient inequality \(c^2\le r_A r_B\). Also, the unsigned tilt is not smooth at symmetry; it is the absolute value of a smooth signed imbalance.

**Numerical only.** Arbitrary precision computations refine the .68/.72 fold to \(.241633835500955\ldots\) and the .6/.8 fold to \(.016546837716879\ldots\). All reference values in `reference_section31.json` agree within \(10^{-3}\). Two earlier rounded estimates, .755 for the unequal-size crossing and .02 for the strong-asymmetry fold, each differ from the refined values by more than that threshold; the article reports the refined values. Numerical verification is separate from the proofs.

**Scope limits.** This resource does not enumerate every zero-loading or sign-changing saddle, does not replace the exact unequal-size root-comparison rule with a closed partition of the whole parameter space, and does not prove a universal bound for the heterogeneous mirrored-design approximation. These limits are listed again in Section 6.

## Sources and verification conventions

The relevant arguments of Online Resource 1 (its Lemmas 2.2, 3.1 to 3.3, 4, 5.1 to 5.2, 7.2, and 8) are rederived here for unequal blocks, so the present resource is self-contained. The file `reference_section31.json` contains earlier numerical reference values and serves only as a comparison target. Its companion reference program was neither imported nor executed and its implementation was not reused.

Apart from the interval certificate of Section 2.6, which is the separate script `checks/unequal_example_certificate.py`, all exact and numerical verification is consolidated in `verify_unequal_blocks.py`. In particular, `scalar_symbolic_checks`, `cusp_symbolic_checks`, `extension_symbolic_checks`, and `unequal_example_certificate` identify the exact algebra being checked. `verify_results.json` records the numerical inputs, tolerances, roots, Hessians, reference differences, and software versions.

The article's followed-trait share is \(\max(S_A,1-S_A)\), where \(S_A=w_A/(w_A+w_B)\). Statements about a smooth crossing of one half refer to \(S_A\), or to its signed imbalance, rather than to the maximum.

## 1. The unrestricted block reduction (Lemma 1)

### 1.0. Scope and notation

**Proved.** Under the hypotheses of Lemma 1 of the article (block A with \(m_A\ge2\) items loading \(\lambda_A\), block B with \(m_B\ge2\) items loading \(\lambda_B\), \(0<\lambda_A,\lambda_B<1\), and \(0<\phi<1\)), the lemma holds for every finite local minimum, including local minima that are not global. Its proof uses neither equal block sizes nor equal within-block loadings. More strongly, every positive interior stationary point is block-uniform. The extension to arbitrary interchangeable correlation matrices is false. A precise common-parts extension is proved in Section 1.7, together with an exact counterexample to the unrestricted extension.

Put \(n_A=m_A-1\), \(n_B=m_B-1\), and write
\[
R=\operatorname{diag}(a_A I_{m_A},a_B I_{m_B})+
\begin{pmatrix}
L_A1_A1_A^T&c1_A1_B^T\\
c1_B1_A^T&L_B1_B1_B^T
\end{pmatrix}.
\]
Here \(a_A=1-L_A>0\), \(a_B=1-L_B>0\), and \(c=\sqrt{L_AL_B}\phi>0\). The contrast eigenvalues of \(R\) are \(a_A\) and \(a_B\), with multiplicities \(n_A\) and \(n_B\). Its restriction to the orthonormal block-mean basis is
\[
\begin{pmatrix}b_A&h\\h&b_B\end{pmatrix},
\qquad h=c\sqrt{m_Am_B}.
\]
In particular,
\[
\begin{split}
b_Ab_B-h^2
&=a_Aa_B+a_Am_BL_B+a_Bm_AL_A
  +m_Am_BL_AL_B(1-\phi^2)>0,\\
\det R&=a_A^{n_A}a_B^{n_B}(b_Ab_B-h^2)>0.
\end{split}
\]
Thus \(R\succ0\). The discrepancy is considered only at
\(\Sigma=ll^T+\Psi\succ0\), where \(\Psi\) is diagonal and nonnegative. A finite local minimum means a local minimum of this feasible parameter problem. The conclusions also hold for local minima in the feasible covariance model, because any sufficiently small parameter perturbation gives a sufficiently small covariance perturbation.

### 1.1. Existence, with boundary representations included

**Proved.** For every positive-definite population covariance matrix \(R\), the one-factor discrepancy attains its infimum over this model.

Let \(\mu_1,\ldots,\mu_p\) be the eigenvalues of
\(R^{1/2}\Sigma^{-1}R^{1/2}\). Then
\[
F(\Sigma)=\sum_{i=1}^p g(\mu_i),
\qquad g(t)=t-\log t-1.
\]
The function \(g\) is nonnegative and tends to infinity at both ends of \((0,\infty)\). On a finite sublevel of \(F\), all \(\mu_i\) therefore lie in a fixed compact subinterval of \((0,\infty)\). Consequently \(\Sigma\) is uniformly bounded above and uniformly positive definite on that sublevel.

The model is closed relative to the positive-definite cone. Indeed, if
\(l^{(k)}l^{(k)T}+\Psi^{(k)}\to\Sigma\succ0\), then
\[
(l_i^{(k)})^2\le\Sigma_{ii}^{(k)},\qquad
0\le\psi_i^{(k)}\le\Sigma_{ii}^{(k)}.
\]
The parameters have a convergent subsequence whose limit is a feasible representation of \(\Sigma\). A minimizing sequence in the finite sublevel containing \(I\) therefore has a feasible positive-definite limit. Continuity of \(F\) proves attainment. No exchangeability or block-size condition enters this argument.

### 1.2. Concentration in the loading vector

**Proved.** Fix a positive diagonal \(\Psi\), and put
\[
S=\Psi^{-1/2}R\Psi^{-1/2},\quad
v=\Psi^{-1/2}l,\quad r=v^Tv.
\]
The rank-one inverse and determinant identities give
\[
F=\operatorname{tr}S-\log\det S-p
 -\frac{v^TSv}{1+r}+\log(1+r).
\tag{1.1}
\]
For fixed \(r\), a minimum over the direction of \(v\) selects a leading eigenvector of \(S\). If \(\zeta_1\ge\cdots\ge\zeta_p\) are its eigenvalues, the derivative in \(r\) after this selection is
\[
\frac{1+r-\zeta_1}{(1+r)^2}.
\]
It follows that
\[
\min_lF=
\begin{cases}
\displaystyle\sum_{j=2}^p g(\zeta_j),&\zeta_1>1,\\
\displaystyle\sum_{j=1}^p g(\zeta_j),&\zeta_1\le1,
\end{cases}
\tag{1.2}
\]
and the optimum has \(r=\max(\zeta_1-1,0)\).

This also constrains every nonzero local minimum over \(l\), not only a global one. If its direction were associated with a smaller eigenvalue, an arbitrarily small rotation toward a larger-eigenvalue direction, keeping \(r\) fixed, would strictly decrease (1.1). Radial stationarity then gives \(r=\zeta_1-1>0\). Neither exchangeability nor any block-size assumption is used.

### 1.3. Interior score equations and positivity

**Proved.** At any interior stationary point, let
\[
K=\Sigma^{-1},\quad H=K-KRK,
\quad u=\Psi^{-1}l,\quad \gamma=l^T\Psi^{-1}l.
\]
Differentiation in \(l\) and the uniquenesses gives
\[
Hl=0,\qquad H_{ii}=0.
\tag{1.3}
\]
Since \(Kl=u/(1+\gamma)\), multiplication of the first equation by \(\Sigma\) yields
\[
Ru=(1+\gamma)l.
\tag{1.4}
\]
Also \(\Sigma u=(1+\gamma)l\). Consequently
\((\Sigma-R)u=0\). Substituting
\(K=\Psi^{-1}-uu^T/(1+\gamma)\) into
\(H=K(\Sigma-R)K\) gives
\[
H=\Psi^{-1}(\Sigma-R)\Psi^{-1}.
\]
The diagonal equations in (1.3) therefore imply
\[
\Sigma_{ii}=1,\qquad \psi_i=1-l_i^2.
\tag{1.5}
\]
Conversely, (1.4) and (1.5) imply both equations (1.3). Subtracting the diagonal term in (1.4) gives the equivalent row form
\[
\sum_{j\ne i}u_j(R_{ij}-l_i l_j)=0.
\tag{1.6}
\]
These identities hold for any correlation matrix \(R\), without exchangeability or equal block sizes.

**Proved.** Every interior local minimum for the present \(0<\phi<1\) design has strictly positive loadings after a simultaneous sign change.

If \(l=0\), uniqueness stationarity gives \(\Psi=I\). Since every off-diagonal entry of \(R\) is positive,
\(1^TR1/p>1\), so \(\lambda_{\max}(R)>1\). Equation (1.1) then decreases for a sufficiently small nonzero loading in its leading direction. Thus \(l=0\) cannot be a local minimum.

If \(l\ne0\), concentration forces its whitened direction to be a leading eigenvector of \(S\). Every entry of \(S\) is strictly positive. The Rayleigh quotient of a mixed-sign vector is strictly increased by replacing it with its coordinatewise absolute values. A nonnegative leading eigenvector has strictly positive coordinates, as follows from its eigenvalue equation. Finally, a leading eigenspace of dimension at least two would contain a nonzero vector orthogonal to a positive leading eigenvector, hence a mixed-sign leading eigenvector. This is impossible. Thus the leading eigenspace is simple and its eigenvectors have one common strict sign. Equation (1.5) now gives
\[
0<l_i<1,\qquad u_i>0.
\tag{1.7}
\]
This argument requires strict positivity of the entries of \(R\), but not equal sizes or loadings.

### 1.4. Complete exclusion of Heywood local minima

**Proved.** A positive-definite one-factor covariance has at most one zero uniqueness. If two coordinates had zero uniqueness, their two-by-two principal submatrix would be rank at most one. Equivalently, a nonzero vector supported on those coordinates and orthogonal to the corresponding loading subvector would have zero quadratic form. Positive definiteness excludes this possibility. If \(\psi_j=0\), then \(l_j\ne0\) for the same reason.

Suppose a finite local minimum has \(\psi_j=0\). All other uniquenesses are positive. On this boundary stratum use coordinates
\[
t=l_j^2>0,\qquad \delta_i=l_i/l_j\quad(i\ne j).
\]
Up to an additive constant independent of these coordinates, block inversion gives
\[
F=\log t+\frac1t+
\sum_{i\ne j}\left[
\log\psi_i+\frac{1-2\delta_iR_{ij}+\delta_i^2}{\psi_i}
\right].
\tag{1.8}
\]
All these coordinates vary freely within the stratum. Its only stationary point, and unique strict minimum in these coordinates, is
\[
t=1,\qquad \delta_i=R_{ij},\qquad
\psi_i=1-R_{ij}^2.
\tag{1.9}
\]
The positive-definiteness of \(R\) gives \(|R_{ij}|<1\), so the displayed uniquenesses are positive. Choose the global sign so that \(l_j=1\), hence \(l_i=R_{ij}\).

At (1.9), \(\Delta=R-\Sigma\) has zero diagonal and zero row and column \(j\), while
\[
\Delta_{ik}=R_{ik}-R_{ij}R_{kj}
\quad(i\ne k,\ i,k\ne j).
\]
Also
\[
K_{ji}=-\frac{R_{ij}}{1-R_{ij}^2}\quad(i\ne j).
\]
The feasible inward perturbation that increases only \(\psi_j\) therefore has derivative
\[
\left.\frac{d}{d\eta}F(\Sigma+\eta e_je_j^T)\right|_{0+}
=-\sum_{\substack{i\ne k\\i,k\ne j}}
w_iw_k\bigl(R_{ik}-R_{ij}R_{kj}\bigr),
\qquad
w_i=\frac{R_{ij}}{1-R_{ij}^2}.
\tag{1.10}
\]
The sum is over ordered pairs. This rederives the boundary criterion of Proposition 4 of the article (Appendix B) and verifies all of its hypotheses here.

If \(j\in A\), the three possible off-diagonal conditional covariances are
\[
\begin{array}{c|c}
\text{locations of }i,k& R_{ik}-R_{ij}R_{kj}\\\hline
A\mathbin{\backslash}\{j\},\ A\mathbin{\backslash}\{j\}&
L_A-L_A^2=L_Aa_A\\
B,\ B&L_B-c^2=L_B(1-L_A\phi^2)\\
A\mathbin{\backslash}\{j\},\ B&c-L_Ac=ca_A.
\end{array}
\tag{1.11}
\]
Every available entry is strictly positive, and all \(w_i\) are strictly positive because \(\phi>0\). In particular, the two distinct items in \(B\) give a strictly positive summand. Thus (1.10) is strictly negative. This remains true when \(m_A=2\), even though there is then no pair in \(A\mathbin{\backslash}\{j\}\). If \(j\in B\), interchange \(A\) and \(B\): the corresponding three entries are
\(L_Ba_B\), \(L_A(1-L_B\phi^2)\), and \(ca_B\), and a pair in \(A\) proves strictness. No boundary point can be a finite local minimum.

### 1.5. At most two loading values and the exact product identity

**Proved.** At a positive interior stationary point, each block has at most two distinct loadings. If two distinct loadings \(x<y\) occur in block \(A\), then \(xy>L_A\).

Let \(U_A=\sum_{i\in A}u_i\), \(U_B=\sum_{i\in B}u_i\), and \(H_0=1+\gamma\). Equation (1.4) in block \(A\) is
\[
H_0l_i-a_A\frac{l_i}{1-l_i^2}=C_A,
\qquad C_A=L_AU_A+cU_B>0.
\tag{1.12}
\]
The function
\[
f_A(t)=H_0t-a_A\frac{t}{1-t^2},\qquad 0<t<1,
\]
has
\[
f_A''(t)=-\frac{2a_At(3+t^2)}{(1-t^2)^3}<0.
\]
Since its derivative is strictly decreasing, a horizontal line meets its graph in at most two points. The same argument applies to block \(B\), with \(a_B\) in place of \(a_A\).

Choose two actual item indices with distinct loadings \(x<y\) in \(A\), and define
\[
P=(1-x^2)(1-y^2),\qquad
\gamma_{\rm rest}=\sum_{k\text{ outside these two indices}}
\frac{l_k^2}{1-l_k^2}.
\]
The elementary divided-difference identity
\[
\frac{y/(1-y^2)-x/(1-x^2)}{y-x}
=\frac{1+xy}{(1-x^2)(1-y^2)}
\]
and \(f_A(x)=f_A(y)\) yield
\(H_0P=a_A(1+xy)\). Direct expansion also gives
\[
H_0P=1-x^2y^2+\gamma_{\rm rest}P.
\]
Subtracting proves exactly
\[
\boxed{
\gamma_{\rm rest}(1-x^2)(1-y^2)
=(1+xy)(xy-L_A).}
\tag{1.13}
\]
Every summand in \(\gamma_{\rm rest}\) is positive. There are at least the \(m_B\ge2\) items of the other block outside the selected pair, even if the selected pair is the whole of a size-two block \(A\). Consequently \(\gamma_{\rm rest}>0\), and (1.13) gives \(xy>L_A\). The analogous result in block \(B\) is \(zw>L_B\).

### 1.6. Strict residual-sign contradiction and conclusion

**Proved.** Every positive interior stationary point is constant within each block.

Suppose \(A\) is mixed, with its two loading values \(x<y\). By (1.13), \(xy>L_A\), hence \(y>\lambda_A\). In the score row (1.6) of an item with loading \(y\), every within-block residual is strictly negative: its fitted product with any other item in \(A\) is at least \(xy>L_A\). At least one such item exists. Since every weight \(u_j\) is positive and the weighted row sum is zero, at least one cross-block residual must be strictly positive. With
\(z=\min_{j\in B}l_j\), this implies
\[
yz<c.
\tag{1.14}
\]
In particular,
\[
z<\frac c y<\frac c{\lambda_A}
=\lambda_B\phi<\lambda_B.
\tag{1.15}
\]
In the row of an item with loading \(z\), all cross-block residuals are strictly positive, because every loading in \(A\) is at most \(y\) and \(zy<c\). If \(B\) were uniform, all of its within-block residuals would equal \(L_B-z^2>0\). The entire row sum would then be strictly positive, contradicting (1.6). Thus \(B\) is also mixed, with its two values \(z<w\). Equation (1.13) in block \(B\) gives
\[
zw>L_B.
\tag{1.16}
\]
Now consider a row with loading \(w\). All its within-block residuals are strictly negative. A cross-block residual must again be strictly positive. Since \(x\) is the smallest loading in \(A\), this implies
\[
wx<c.
\tag{1.17}
\]
But
\[
L_AL_B<(xy)(zw)=(yz)(wx)<c^2
=L_AL_B\phi^2<L_AL_B,
\]
a contradiction. If \(B\) is initially the mixed block, the identical argument starts with \(B\). Thus neither block is mixed.

Combining Sections 1.3 through 1.6 proves Lemma 1: every finite local minimum is interior and, after one global sign change, has loadings \(\alpha\in(0,1)\) on \(A\), \(\beta\in(0,1)\) on \(B\), and uniquenesses \(1-\alpha^2\), \(1-\beta^2\).

**Proved, dependency audit.** Existence, concentration, and the interior scores need only \(R\succ0\) with unit diagonal for the score normalization. Positivity additionally uses strictly positive entries and \(\lambda_{\max}(R)>1\). None uses equality of sizes or loadings. The condition that both block sizes are at least two supplies an opposite-block pair in the strict Heywood derivative and supplies positive terms outside any chosen within-block pair in (1.13). A mixed block itself necessarily contains at least two items, which makes its high-loading row contain a strictly negative within-block term. The proof remains valid for arbitrarily small positive \(\phi\), for \(\phi\) arbitrarily close to one, and for any fixed \(0<\lambda_A,\lambda_B<1\). It does not assert uniform quantitative bounds as a loading tends to zero or one. At a zero loading strict positivity fails; at a loading equal to one the population matrix can be singular. The endpoint \(\phi=0\) needs the separate treatment below.

### 1.7. Scope of the interchangeable-items extension and a counterexample

**Proved.** The block-uniformity lemma extends to every correlation matrix with unit diagonal, constant within-block off-diagonal correlations
\(r_A,r_B\in(0,1)\), and constant nonzero cross-block correlation \(c\), provided
\[
c^2\le r_Ar_B.
\tag{1.18}
\]
For \(c>0\), replace \(L_A,L_B\) by \(r_A,r_B\) throughout the proof. The strict inequality \(c^2<r_Ar_B\) is not necessary. At equality, (1.15) still gives \(z<\sqrt{r_B}\), since \(y>\sqrt{r_A}\), and the final contradiction still has
\((xy)(zw)>r_Ar_B\) and \((yz)(wx)<c^2=r_Ar_B\). The conditional covariances used in the Heywood argument remain strictly positive because \(r_A,r_B<1\).

The common-parts construction described in the article implies (1.18). Suppose standardized items have the form
\[
X_i=U_A+e_i\ (i\in A),\qquad
X_j=U_B+e_j\ (j\in B),
\]
where all variables are centered, residuals are mutually uncorrelated and uncorrelated with both common parts, and the residual variances are strictly positive. If
\(\operatorname{Var}(U_A)=r_A>0\),
\(\operatorname{Var}(U_B)=r_B>0\), and
\(\operatorname{Cov}(U_A,U_B)=c\), then Cauchy-Schwarz gives
\(c^2\le r_Ar_B\). Equality holds exactly when the common parts are perfectly correlated or perfectly negatively correlated, equivalently when
\(U_B=(c/r_A)U_A\) almost surely. The inequality is strict when their correlation has absolute value less than one.

This includes equal cross-loading constructions after proper item standardization. If the common parts are \(v_A^T\theta\) and \(v_B^T\theta\), with trait covariance \(\Omega\), the same inequality is Cauchy-Schwarz for \(\Omega^{1/2}v_A\) and \(\Omega^{1/2}v_B\). Equality occurs exactly when these nonzero vectors are linearly dependent.

If \(c<0\), put \(D=\operatorname{diag}(I_{m_A},-I_{m_B})\). Then \(DRD\) has cross correlation \(|c|\),
\[
D\Sigma D=(Dl)(Dl)^T+\Psi,
\qquad
F_{DRD}(D\Sigma D)=F_R(\Sigma).
\]
This is an invertible, continuous transformation of the whole feasible parameter problem, including its boundary, so it preserves every local and global minimum. Reflection therefore fully handles the negative cross-correlation case in the stated extension: the two fitted blocks have opposite loading signs before reflection. At \(c=0\), strict positivity is unavailable and the endpoint theorem below applies instead.

**Counterexample.** Interchangeability alone does not suffice, even when all entries of \(R\) are positive. Consider
\[
R=\begin{pmatrix}
1&1/10&1/2&1/2\\
1/10&1&1/2&1/2\\
1/2&1/2&1&1/10\\
1/2&1/2&1/10&1
\end{pmatrix}.
\tag{1.19}
\]
This matrix has interchangeable size-two blocks and eigenvalues
\(9/10,9/10,21/10,1/10\), so it is positive definite. It violates (1.18), because \(c^2=1/4>1/100=r_Ar_B\).

The star representation
\[
l=(1,1/10,1/2,1/2)^T,\qquad
\Psi=\operatorname{diag}(0,99/100,3/4,3/4)
\tag{1.20}
\]
is a strict local minimum of the full feasible problem. To prove this exactly, first note that it is the unique strict minimum on the face \(\psi_1=0\), by (1.8). Its two cross-pair conditional residuals are \(9/20\), while its conditional residual between items 3 and 4 is \(-3/20\). The weights are \(w_2=10/99\), \(w_3=w_4=2/3\). Equation (1.10) gives
\[
\left.\frac{dF}{d\psi_1}\right|_{0+}
=-\left[
4\frac{10}{99}\frac23\frac9{20}
 +2\left(\frac23\right)^2\left(-\frac3{20}\right)
\right]
=-\left(\frac4{33}-\frac2{15}\right)
=\frac2{165}>0.
\tag{1.21}
\]
For completeness, let \(z\) denote the free face coordinates in (1.8), and \(\eta=\psi_1\ge0\). The objective is smooth near (1.20), because its covariance is positive definite. By continuity, \(\partial_\eta F(z,\eta)>1/165\) throughout a sufficiently small neighborhood. Therefore
\[
F(z,\eta)\ge F(z,0)+\eta/165\ge F(z_0,0),
\]
with strict inequality unless \((z,\eta)=(z_0,0)\). This proves a strict local minimum, including directions combining an inward move with changes along the face.

It is also a strict local minimum in covariance space. All off-diagonal entries of its fitted covariance are positive, and a nearby feasible covariance has, after its global sign is fixed, the continuous recovery formulas
\(l_1^2=\Sigma_{12}\Sigma_{13}/\Sigma_{23}\),
\(l_i=\Sigma_{1i}/l_1\) for \(i\ne1\), and
\(\psi_i=\Sigma_{ii}-l_i^2\). Thus nearby covariances cannot escape the parameter neighborhood used above.

It is even a nonglobal boundary minimum. Its discrepancy is
\(\log(275/84)\), while the feasible uniform covariance with loading square \(11/30\) and uniqueness \(19/30\) has the smaller discrepancy
\(\log(6859/2187)\). The strict comparison follows from
\(6859\cdot84=576156<601425=275\cdot2187\).

## 2. The exact two-parameter problem

### 2.1. The zero-correlation endpoint

#### 2a.1. Global value and selected block

**Proved.** At \(\phi=0\), write the total correlations (Gaussian multi-information) of the within-block correlation matrices \(R_A\) and \(R_B\) as
\[
\begin{split}
\mathrm{TC}(R_A)&=-\tfrac12\log\det R_A=-\tfrac12\left(n_A\log a_A+\log b_A\right),\\
\mathrm{TC}(R_B)&=-\tfrac12\log\det R_B=-\tfrac12\left(n_B\log a_B+\log b_B\right).
\end{split}
\]
Then
\[
\boxed{F_{\min}(0)=2\min\{\mathrm{TC}(R_A),\mathrm{TC}(R_B)\}.}
\tag{2a.1}
\]
If \(\mathrm{TC}(R_A)>\mathrm{TC}(R_B)\), the unique minimizing covariance is
\(\operatorname{diag}(R_A,I_{m_B})\). If \(\mathrm{TC}(R_B)>\mathrm{TC}(R_A)\), it is
\(\operatorname{diag}(I_{m_A},R_B)\). Thus the fitted factor retains the block with larger Gaussian multi-information, and the minimum discrepancy is the cost of leaving the other block diagonal.

First assume all uniquenesses are positive. At \(\phi=0\), the whitened population \(S\) is block diagonal. If its largest eigenvalue exceeds one, concentration can select a leading eigenvector supported entirely in one block, say \(A\). This remains possible when the largest eigenvalue is tied between blocks. Its discrepancy splits as
\[
F_{\rm conc}
=F_A(l_A,\Psi_A)+F_{B,\rm diag}(\Psi_B),
\]
where
\[
F_A(l_A,\Psi_A)\ge0,
\qquad
F_{B,\rm diag}(\Psi_B)
=\sum_{i\in B}g(1/\psi_i)+2\,\mathrm{TC}(R_B)\ge2\,\mathrm{TC}(R_B).
\tag{2a.2}
\]
If the leading eigenvector is selected from \(B\), the same argument gives \(F_{\rm conc}\ge2\,\mathrm{TC}(R_A)\). Hence every positive-uniqueness fit has discrepancy at least \(2\min\{\mathrm{TC}(R_A),\mathrm{TC}(R_B)\}\). If all whitened eigenvalues are at most one, concentration gives a fully diagonal model, whose discrepancy is at least \(2\,\mathrm{TC}(R_A)+2\,\mathrm{TC}(R_B)\), strictly larger than the claimed minimum. To check this strictness, note that \(\mathrm{TC}(R_A),\mathrm{TC}(R_B)>0\): for a block with \(m\ge2\) items and within-block correlation \(L\), the quantity \(C(L)=-\log\det R=2\,\mathrm{TC}(R)\) of its correlation matrix \(R\) satisfies
\[
C(0)=0,\qquad
\frac{dC(L)}{dL}=(m-1)
\left(\frac1{1-L}-\frac1{1+(m-1)L}\right)>0
\quad(0<L<1).
\]

The bound also holds at any positive-definite boundary representation, by replacing \(\Psi\) with \(\Psi+\delta I\), applying the proved bound for \(\delta>0\), and taking \(\delta\downarrow0\). Finally, fitting either block exactly and making the other block diagonal attains \(2\,\mathrm{TC}(R_B)\) or \(2\,\mathrm{TC}(R_A)\), respectively. This proves (2a.1).

#### 2a.2. All equality cases, including unequal sizes and two-item families

**Proved.** If \(\mathrm{TC}(R_A)=\mathrm{TC}(R_B)\), there are exactly two minimizing covariance matrices,
\[
\operatorname{diag}(R_A,I_{m_B})
\quad\text{and}\quad
\operatorname{diag}(I_{m_A},R_B).
\tag{2a.3}
\]
There is no flat family of different minimizing covariance matrices connecting them. The only additional flat families are the familiar parameter nonidentifiabilities of an exactly fitted two-item retained block.

To classify interior equality, every global minimum must already be concentrated in \(l\) for its fixed \(\Psi\). If the leading eigenspace belongs to only one block, its loading vector has support only in that block. Equality in (2a.2) then forces that block's covariance to equal its population covariance and the other uniqueness matrix to equal the identity.

A tie of leading eigenvalues between blocks is impossible at a global minimum. Under a tie, concentrated pure-\(A\) and pure-\(B\) loading vectors have the same discrepancy. Equations (2a.2) then require this common discrepancy to be at least both \(2\,\mathrm{TC}(R_A)\) and \(2\,\mathrm{TC}(R_B)\). If the two costs differ, that contradicts (2a.1). If they are equal, both lower bounds must be equalities. Equality for the pure-\(A\) vector forces \(\Psi_B=I\), while equality for the pure-\(B\) vector forces \(\Psi_A=I\). But equality for the pure-\(A\) vector would then require
\(R_A=I+l_Al_A^T\). Its unit diagonal forces \(l_A=0\), which contradicts its strictly positive within-block correlation. Thus no minimizing leading eigenspace mixes the blocks.

For a boundary global minimum, Section 1.4 still gives at most one zero uniqueness and requires the star form (1.9). If its index lies in \(A\) and \(m_A\ge3\), a pair in \(A\mathbin{\backslash}\{j\}\) has positive weights and conditional residual \(L_Aa_A\), so (1.10) is strictly negative even at \(\phi=0\). Such a boundary point cannot minimize. If \(m_A=2\), the star covariance fits \(R_A\) exactly and gives the identity on \(B\), hence is already among the covariances above. The same argument applies with the blocks interchanged. This excludes every other boundary equality case.

For a retained block with \(m\ge3\), its exact off-diagonal equations
\(l_i l_j=L>0\) imply that all retained loadings have one common sign. Taking three distinct indices gives
\[
l_i^2=\frac{(l_il_j)(l_il_k)}{l_jl_k}=L.
\]
Consequently the representation is the uniform loading \(\lambda=\sqrt L\) with uniqueness \(1-L\), up to the global sign. For a retained block of size two, the complete family is instead
\[
\boxed{
(l_1,l_2)=(t,L/t),\qquad L\le t\le1,
\qquad
(\psi_1,\psi_2)=(1-t^2,1-L^2/t^2).}
\tag{2a.4}
\]
Include simultaneous negatives. The other block has zero loadings and unit uniquenesses. The endpoints \(t=L\) and \(t=1\) each have one zero uniqueness. All members of (2a.4) represent exactly the same fitted covariance. If \(\mathrm{TC}(R_A)=\mathrm{TC}(R_B)\) with unequal sizes, the same classification applies independently to either retained block: one of the two covariance minima may have a parameter family while the other has an identifiable representation.

#### 2a.3. Limit from positive trait correlations

**Proved.** If \(\mathrm{TC}(R_A)\ne\mathrm{TC}(R_B)\), every global minimizing covariance for \(\phi>0\) converges as \(\phi\downarrow0\) to the unique covariance specified in Section 2a.1. With the positive block-uniform loading convention, if \(\mathrm{TC}(R_A)>\mathrm{TC}(R_B)\),
\[
\alpha(\phi)\longrightarrow\lambda_A,
\qquad\beta(\phi)\longrightarrow0,
\]
and the analogous conclusion holds if \(\mathrm{TC}(R_B)>\mathrm{TC}(R_A)\). This is true even if the retained block has size two: the positive-\(\phi\) block-uniform fits select the uniform member of its endpoint parameter family.

For proof, \(R(\phi)\) stays uniformly positive definite in a neighborhood of zero, and a fixed feasible exact-one-block fit provides a uniform finite upper bound on the optimum. The coercivity argument of Section 1.1 therefore places all global minimizing covariances in a common compact subset of the positive-definite cone. Every convergent subsequence has an endpoint limit. Continuity of the objective and comparison with each fixed endpoint minimizer show that every such limit is an endpoint global minimizer. When \(\mathrm{TC}(R_A)\ne\mathrm{TC}(R_B)\), its covariance is unique, so the whole family converges. Block uniformity, the unit fitted diagonal, and the positive sign convention determine the displayed loading limits. If \(\mathrm{TC}(R_A)=\mathrm{TC}(R_B)\), every cluster covariance is one of the two covariances (2a.3); this endpoint argument alone does not determine which continuation is globally selected for small positive \(\phi\).

#### 2a.4. Useful endpoint local-stability criterion

**Proved.** At \(\phi=0\), consider the uniform exact-\(A\), zero-\(B\) stationary fit
\[
l_A=\lambda_A1_A,\quad\Psi_A=a_AI,
\qquad l_B=0,\quad\Psi_B=I.
\]
Put \(\kappa_A=m_AL_A/a_A\). The Hessian block for perturbations in the inactive block's loadings is
\[
\boxed{H_{l_Bl_B}=2\left(I-\frac{R_B}{1+\kappa_A}\right).}
\tag{2a.5}
\]
Its smallest eigenvalue is
\(2[1-b_B/(1+\kappa_A)]\). These loading perturbations have zero Hessian cross terms with all other parameter directions. If \(m_A\ge3\), the Hessian on the remaining directions is positive definite. Thus the fit is a strict local minimum if
\[
1+\kappa_A>b_B,
\tag{2a.6}
\]
and is a saddle if the inequality is reversed. At equality, this second-order calculation is degenerate and does not alone classify the point. If the retained block has size two, its exact-fit parameter family gives a zero Hessian direction, so the positive-definite conclusion for the full Hessian needs the stated \(m_A\ge3\) condition.

To derive (2a.5), hold the uniquenesses fixed and append an inactive-block loading perturbation \(v\). In (1.1), the initial whitened norm is \(\kappa_A\), its Rayleigh numerator is \(\kappa_A(1+\kappa_A)\), and the perturbed quantities add \(\|v\|^2\) and \(v^TR_Bv\), respectively. The quadratic change in the objective is exactly
\[
v^T\left(I-\frac{R_B}{1+\kappa_A}\right)v.
\]
At \(\phi=0\), the gradient in \(l_B\) vanishes whenever \(l_B=0\), for every nearby value of the other parameters. Its Hessian cross terms therefore vanish. The active-block exact-fit discrepancy has positive quadratic form in every nonzero covariance perturbation. For \(m_A\ge3\), the differential of the active one-factor parametrization is injective: its off-diagonal equations imply \(\delta l_i+\delta l_j=0\) for every distinct pair after division by \(\lambda_A\), which forces all \(\delta l_i=0\), and the diagonal equations then force all \(\delta\psi_i=0\). Finally, the inactive diagonal-model uniqueness Hessian is positive definite. This proves the remaining statements.

For the particular endpoint \(m_A=m_B=4\), \(\lambda_A=.6\), \(\lambda_B=.8\), (2a.6) reads \(3.25>2.92\). The weaker-block stationary fit is therefore a strict local minimum at \(\phi=0\), and the implicit function theorem continues it as a strict interior local minimum for all sufficiently small positive \(\phi\). In particular, the possibility that this minimum exists only at \(\phi=0\) is ruled out analytically. This does not by itself locate or classify its eventual fold.

### 2.2. Reduced discrepancy and global scale profiling

**Proved.** Let \(\Psi=\operatorname{diag}(\psi_AI_{m_A},\psi_BI_{m_B})\) be block-uniform, and let \(\zeta_+\ge\zeta_-\) be the eigenvalues of the matrix \(M(\psi_A,\psi_B)\) displayed below. Whenever \(\zeta_+>1\), concentration over the loadings gives the reduced discrepancy
\[
F(\psi_A,\psi_B)=(m_A-1)\,g(a_A/\psi_A)+(m_B-1)\,g(a_B/\psi_B)+g(\zeta_-),
\]
and for this two-trait population no additional eigenvalue-ordering restriction is needed. For arbitrary \(\psi_A,\psi_B>0\), the whitened contrast eigenvalues are \(a_A/\psi_A\), repeated \(m_A-1\) times, and \(a_B/\psi_B\), repeated \(m_B-1\) times. The other eigenvalues are those of
\[
M(\psi_A,\psi_B)=\begin{pmatrix}b_A/\psi_A&h/\sqrt{\psi_A\psi_B}\\h/\sqrt{\psi_A\psi_B}&b_B/\psi_B\end{pmatrix},
\qquad h=c\sqrt{m_A m_B}.
\]
Because \(h>0\), the larger eigenvalue of \(M\) is simple and exceeds both diagonal entries, which in turn exceed their respective contrast eigenvalues. Thus it is always the leading eigenvalue of the whole whitened matrix. Concentration (Section 1.2) gives the reduced discrepancy whenever this eigenvalue exceeds one.

For \(\psi_A=e^{s+d}\), \(\psi_B=e^{s-d}\), put \(\nu_{AB}=m_A+m_B-1\), and let
\[
r_\pm(d)=\frac{b_Ae^{-d}+b_Be^d
\mathbin{\pm}\sqrt{(b_Ae^{-d}-b_Be^d)^2+4h^2}}2,
\]
\[
T(d)=(m_A-1)a_Ae^{-d}+(m_B-1)a_Be^d+r_-(d).
\]
The exact active concentrated discrepancy is
\[
F(s,d)=e^{-s}T(d)+\nu_{AB}s+(m_A-m_B)d
-(m_A-1)\log a_A-(m_B-1)\log a_B-\log r_-(d)-\nu_{AB}.
\tag{R.1}
\]
Every nonleading whitened eigenvalue before scaling is strictly below \(r_+\). Hence \(T<\nu_{AB} r_+\). The unique active scale optimum is
\[
e^{\widehat s(d)}=T(d)/\nu_{AB},
\qquad F_{ss}(\widehat s,d)=\nu_{AB}>0,
\]
and it satisfies \(e^{-\widehat s}r_+>1\). On the inactive side \(s\ge\log r_+\), the diagonal objective has derivative
\[
p-e^{-s}\operatorname{tr}S_0>0,
\]
because \(\operatorname{tr}S_0<p r_+\). Reducing \(s\) improves every such point. Thus the active optimum is the global scale optimum, including both regimes. The resulting exact scalar profile is
\[
\boxed{\mathcal G(d)=\nu_{AB}\log\frac{T(d)}{\nu_{AB}}+(m_A-m_B)d-\log r_-(d)
-(m_A-1)\log a_A-(m_B-1)\log a_B.}
\tag{R.2}
\]
This also gives the direct transcendental scalar equation
\[
\nu_{AB}\frac{T'(d)}{T(d)}+(m_A-m_B)-\frac{r_-'(d)}{r_-(d)}=0.
\]
The quartic below is a rational parametrization of this same profile, with explicit reconstruction and no extraneous positive roots.

### 2.3. One scalar equation, a quartic, and an exact selection rule

**Proved.** With \(\nu_{AB}=m_A+m_B-1\) as above, introduce the coefficients

\[
 k_A=(m_A-1)a_A=m_A-b_A,\qquad k_B=(m_B-1)a_B=m_B-b_B,
 \quad c_B=m_Bc,\quad c_A=m_Ac,
\]
\[
 e_0=m_Ab_B+m_Bb_A-b_Ab_B-c_Bc_A=m_Am_B(1-c^2)-k_Ak_B.
\]
For \(v>0\) define

\[
 U=b_A+c_Bv,\quad V=b_B+c_A/v,\quad
 N=m_AV+m_BU-UV=e_0+c_Bk_Bv+c_Ak_A/v.
\]
Here \(v\) is the ratio of precision-weighted block loadings,

\[
 v=\frac{\beta/\psi_B}{\alpha/\psi_A}.
\]
The scale-profiled uniquenesses and largest whitened eigenvalue are

\[
 \psi_A=\frac{N}{\nu_{AB}V},\qquad \psi_B=\frac{N}{\nu_{AB}U},\qquad
 \tau=\frac{\nu_{AB}UV}{N}.
\]
The exact profiled discrepancy is

\[
\begin{split}
 \mathcal G(v)={}&\nu_{AB}\log N-(m_B-1)\log U-(m_A-1)\log V-\nu_{AB}\log\nu_{AB}\\
 &-\log(b_Ab_B-m_Am_Bc^2)
 -(m_A-1)\log a_A-(m_B-1)\log a_B.
\end{split}                                                    \tag{S1}
\]
The scalar coordinates in (R.2) and (S1) are related by
\[
d(v)=\frac12\log\frac{\psi_A}{\psi_B}=\frac12\log\frac{U}{V},
\qquad d'(v)=\frac12\left(\frac{c_B}U+\frac{c_A}{v^2V}\right)>0.
\]
For positive trait correlation, \(c_A,c_B>0\), and \(d(v)\) tends to
\(-\infty\) as \(v\downarrow0\) and to \(+\infty\) as
\(v\to\infty\). Thus this is a smooth change of coordinate onto
\(\mathbb R\). At stationary points it preserves root multiplicity and
the sign of the scalar second derivative.

Its stationary equation is the quartic

\[
 P(v)=v\{v^2V(\nu_{AB}V-N)-U(\nu_{AB}U-N)\}=0.              \tag{S2}
\]
The expression on the right is a polynomial although its displayed factors contain \(1/v\). To give its coefficients explicitly, put

\[
 J=\nu_{AB}b_A-e_0,\quad K=\nu_{AB}b_B-e_0,\quad f=\nu_{AB}-k_B,\quad g=\nu_{AB}-k_A.
\]
Then, in descending powers of \(v\),

\[
\begin{split}
 P_4&=-b_Bc_Bk_B,\\
 P_3&=b_BK-c_Ac_Bk_B-c_B^2f,\\
 P_2&=b_Bc_Ag+c_AK-b_Ac_Bf-c_BJ,\\
 P_1&=c_A^2g-b_AJ+c_Bc_Ak_A,\\
 P_0&=b_Ac_Ak_A.
\end{split}                                                    \tag{S3}
\]
In particular, \(P_4<0\) and \(P_0>0\).

*Derivation.* The loading score gives

\[
 \tau\psi_A=b_A+m_Bc v=U,\qquad \tau\psi_B=b_B+m_Ac/v=V.
\]
The diagonal scores give

\[
 \tau=1+m_A\alpha^2/\psi_A+m_B\beta^2/\psi_B=m_A/\psi_A+m_B/\psi_B-\nu_{AB}.
\]
These identities give the displayed \(\psi_A\), \(\psi_B\) and \(\tau\). Conversely, the definition of \(v\) requires

\[
 v^2=\frac{1-\psi_B}{1-\psi_A}\left(\frac{\psi_A}{\psi_B}\right)^2
     =\frac{1-\psi_B}{1-\psi_A}\left(\frac UV\right)^2.
\]
Substituting \(\psi_A\) and \(\psi_B\) gives (S2).

For a derivation independent of the diagonal score, use uniquenesses \(tU\) and \(tV\). Their whitened block-mean matrix has largest eigenvalue \(1/t\), with a positive eigenvector, and smaller eigenvalue \((b_Ab_B-c_Bc_A)/(tUV)\). Concentrating the factor and then \(t\) gives \(t=N/(\nu_{AB}UV)\), and the logarithmic terms give (S1).

**Proved. Positive-root admissibility.** Every positive root of (S2) has \(0<\psi_A,\psi_B<1\) and \(\tau>1\). First, \(N>0\): \(k_A,k_B>0\) and

\[
 e_0=m_Am_B(1-c^2)-k_Ak_B>0,
\]
for \(c^2<L_AL_B\) and \(k_Ak_B<m_Am_Ba_Aa_B\), while \(1-L_AL_B>a_Aa_B\). Since \(U>b_A>1\) and \(V>b_B>1\),

\[
 0<\frac N{UV}=\frac{m_A}U+\frac{m_B}V-1<\nu_{AB},
\]
so \(\tau>1\). Equation (S2) implies that \(1-\psi_A\) and \(1-\psi_B\) have the same sign. If both were nonpositive, the trace identity \(\tau=m_A/\psi_A+m_B/\psi_B-\nu_{AB}\) would give \(\tau\le1\), a contradiction. A zero of either forces both to vanish and gives the same contradiction. Hence \(\psi_A,\psi_B<1\). Set \(\alpha=\sqrt{1-\psi_A}\), \(\beta=\sqrt{1-\psi_B}\). The quartic identity recovers the positive ratio \(v\), and the two loading equations follow. The diagonal scores follow from the displayed trace identity. Thus every positive root is a genuine positive block-uniform stationary point, and every such point is represented.

**Proved with computer algebra.** Let \(\bar V=vV\) and \(\bar N=vN\). Direct differentiation gives

\[
 \mathcal G'(v)=-\frac{cm_Am_B\,P(v)}{vU\bar V\bar N}.       \tag{S4}
\]
`verify_unequal_blocks.py::scalar_symbolic_checks` checks this rational identity, as well as the coefficient identity (S3), exactly with SymPy. The denominator is positive for \(v>0\).

**Proved. Candidate count and exact global selection.** There are at most three positive roots counted with multiplicity. There is always at least one positive root and at least one negative root, since \(P_4<0<P_0\). Thus, away from degeneracies, there are either one or three positive roots. The profile is coercive at both ends of \(v>0\). Its derivative starts negative and ends positive. Consequently the simple positive roots, in increasing order, have types

\[
 \text{minimum},\quad \text{saddle},\quad \text{minimum},
\]
with only the first type when there is one root. Section 2.4 proves that these are their types in the full parameter space. The global solution is obtained by evaluating (S1) on the one or two minimum roots. This is an exact finite algebraic and logarithmic selection rule for arbitrary unequal block sizes.

An exact parameter-dependent count can be stated without a large case partition: form the Sturm sequence of the explicit quartic (S3), and subtract its sign variation at positive infinity from its variation at 0. This gives the number of distinct positive roots, with the standard gcd step identifying multiplicities. The discriminant by itself does not determine the number of positive roots without a sign/root-location condition.

**Proved. Degenerate positive roots.** A double positive root has a cubic leading term in the local profile, so it is not a local minimum at the exact collision. A triple positive root is the only positive root: writing \(P(v)=-a(v-v_0)^3(v-r_-)\), where \(a>0\) and \(r_-<0\), the sign identity (S4) gives a positive quartic leading term in \(\mathcal G-\mathcal G(v_0)\). It is therefore a strict global minimum with one zero full-Hessian eigenvalue. There is no positive root of multiplicity four, because the quartic has at least one negative root. This completes the classification at degeneracies as well as at simple roots.

A fold satisfies

\[
 P(v;\phi)=P_v(v;\phi)=0,\qquad v>0,\quad0<\phi<1. \tag{S5}
\]
Equivalently, the resultant \(\operatorname{Res}_v(P,P_v)\) vanishes, together with the positive double-root and parameter-domain conditions. A zero of the resultant alone can concern another real root or a boundary parameter and must not automatically be called a fold. For a nondegenerate fold as \(\phi\) varies, also require \(P_{vv}\ne0\) and \(P_\phi\ne0\) at the double root. These transversality conditions are proved below for the equal-size weaker-branch fold.

**Scope limitation.** This is a complete classification of positive block-uniform stationary points and, given Lemma 1 (Section 1), of all local minima. It does not enumerate every saddle of the unrestricted \(2p\)-parameter likelihood. In particular, zero-loading and sign-changing stationary points can lie outside this scalar family.

### 2.4. The full Hessian, including within-block contrasts

**Proved.** All within-block contrast sectors are positive definite at every positive root of (S2).

First, the row scores give

\[
 (m_A-1)\frac{\alpha}{\psi_A}(L_A-\alpha^2)
   +m_B\frac{\beta}{\psi_B}(c-\alpha\beta)=0,
\]
\[
 (m_B-1)\frac{\beta}{\psi_B}(L_B-\beta^2)
   +m_A\frac{\alpha}{\psi_A}(c-\alpha\beta)=0.
\]
If \(\alpha\beta\le c\), these imply \(\alpha^2\ge L_A\) and \(\beta^2\ge L_B\), hence \(\alpha\beta\ge\sqrt{L_AL_B}>c\), a contradiction. Therefore

\[
 \alpha\beta>c,\qquad 0<\alpha^2<L_A,\qquad0<\beta^2<L_B. \tag{S6}
\]

Take a vector \(e\) supported on A with \(\sum_i e_i=0\) and \(\sum_i e_i^2=1\). Put \(\delta l=te\) and \(\delta\psi=we\). In these coordinates the exact Hessian quadratic form is

\[
 H_A(t,w)=\frac{2t^2}{\psi_A}+\frac{w^2}{\psi_A^2}
 -\frac{2a_A}{\tau}\left(\frac t{\psi_A}-\frac{\alpha w}{\psi_A^2}\right)^2.
                                                               \tag{S7}
\]
Its determinant is

\[
 \det H_A=\frac2{\psi_A^3}
 \left[1-\frac{a_A(1+\alpha^2)}{\tau \psi_A^2}\right]. \tag{S8}
\]
The needed strict positivity follows from the identity

\[
 \tau \psi_A^2-a_A(1+\alpha^2)
 =(L_A-\alpha^2)(1+\alpha^2)
 +(m_A-2)\alpha^2\psi_A+m_B\beta^2\psi_A^2/\psi_B>0.              \tag{S9}
\]
This includes \(m_A=2\). The first principal minor is also positive. The same calculation applies to the \(m_B-1\) contrast copies in B.

*Derivation of (S7).* Use \(u=\Psi^{-1}l\) and \(\gamma=l^T\Psi^{-1}l\) of Section 1.3, and put \(\tau=1+\gamma\) and \(\mathcal Q=u^TRu\). In the notation of (1.1), \(r=\gamma\) and \(\mathcal Q=v^TSv\). At a stationary point \(\tau\) is the largest whitened eigenvalue of Section 2.3. Since \(\operatorname{tr}S=\sum_i1/\psi_i\) and \(-\log\det S=\sum_i\log\psi_i-\log\det R\), (1.1) reads
\[
F=\sum_i\Bigl(\frac1{\psi_i}+\log\psi_i\Bigr)+\log\tau-\frac{\mathcal Q}\tau-\log\det R-p.
\]
Move along \(l+\eta te\), \(\psi+\eta we\), and write primes for derivatives in \(\eta\) at \(\eta=0\). The loadings and uniquenesses are constant on the support of \(e\), so each derivative below reduces to a multiple of \(\sum_ie_i=0\) or \(\sum_ie_i^2=1\). With \(\kappa=t/\psi_A-\alpha w/\psi_A^2\),
\[
\gamma'=0,\qquad \gamma''=2\psi_A\kappa^2,\qquad u'=\kappa e,\qquad u_i''=-\frac{2w\kappa}{\psi_A}e_i^2 .
\]
Because \(\tau'=\gamma'=0\), the second derivative of \(\log\tau-\mathcal Q/\tau\) is \(\gamma''/\tau-\mathcal Q''/\tau+\mathcal Q\gamma''/\tau^2\). At a stationary point the loading equation (1.4) is \(Ru=\tau l\), so \(\mathcal Q=\gamma\tau\), and this second derivative equals \(\gamma''-\mathcal Q''/\tau\). A within-block contrast satisfies \(Re=a_Ae\), hence
\[
\mathcal Q''=2u'^TRu'+2(Ru)^Tu''=2a_A\kappa^2-\frac{4\tau\alpha w\kappa}{\psi_A}.
\]
The diagonal terms contribute \(2w^2/\psi_A^3-w^2/\psi_A^2\). Using \(2\psi_A\kappa^2+4\alpha w\kappa/\psi_A=2t^2/\psi_A-2\alpha^2w^2/\psi_A^3\), the terms add to
\[
H_A(t,w)=\frac{2t^2}{\psi_A}+\frac{2(1-\alpha^2)w^2}{\psi_A^3}-\frac{w^2}{\psi_A^2}-\frac{2a_A}{\tau}\kappa^2 ,
\]
and the diagonal relation (1.5), \(1-\alpha^2=\psi_A\), gives (S7). The determinant of (S7) is \((2/\psi_A^3)[1-a_A(\psi_A+2\alpha^2)/(\tau\psi_A^2)]\). It equals (S8) because \(\psi_A+2\alpha^2=1+\alpha^2\); `verify_unequal_blocks.py::scalar_symbolic_checks` checks this identity exactly. Section 5 compares (S7), transformed to log-uniqueness coordinates, with numerically differentiated second derivatives of the full-covariance discrepancy at 65 digits, for one contrast in each block.

*Sector decomposition.* Permuting the items within A and within B, with loadings and uniquenesses permuted together, leaves \(R\), \(F\) and every block-uniform point unchanged. So at such a point the Hessian \(\mathcal H\) of \(F\) in \((l,\psi)\) commutes with these permutations. Write \(\mathbb R^p=W\oplus C_A\oplus C_B\), where \(W\) is spanned by \(1_A\) and \(1_B\) and \(C_A,C_B\) are the within-block contrast spaces, of dimensions \(m_A-1\) and \(m_B-1\). Loading and uniqueness perturbations split in the same way. Invariance under all permutations of A forces each \(m_A\times m_A\) block of \(\mathcal H\) on the items of A (loading–loading, loading–uniqueness and uniqueness–uniqueness) to have the form \(xI+yJ\), which acts on \(C_A\) as the scalar \(x\). Invariance under the permutations of both blocks forces each block coupling the items of A with those of B to be a multiple of \(1_A1_B^T\), which annihilates \(C_A\) and \(C_B\). Hence \(\mathcal H\) has no cross terms between the four-dimensional block-mean sector \(W\times W\) and the contrast sectors, nor between the A and B contrast sectors. In an orthonormal contrast basis it equals \(I_{m_A-1}\otimes\mathcal H_A\) on \(C_A\times C_A\), with a \(2\times2\) matrix \(\mathcal H_A\), and likewise for B. Equivalently, \(C_A\) and \(C_B\) are absolutely irreducible representations of these permutations (the sign representation when a block has two items). They are inequivalent to each other and to the trivial representation, which makes up \(W\) twice, and Schur's lemma applies to the isotypic components. At a positive root of (S2), \(\mathcal H_A\) is the matrix of (S7), so by (S8) and (S9) the contrast sectors are positive definite.

*Inertia of the block-mean sector.* Use coordinates \((\alpha,\beta,\psi_A,\psi_B)\) on \(W\times W\). This linear change of coordinates preserves inertia. At fixed uniquenesses write \(x=\Psi^{-1/2}l\), the vector \(v\) of Section 1.2, and \(\zeta_1>\zeta_2\ge\dots\) for the eigenvalues of \(S\), with eigenvectors \(q_1,q_2,\dots\). At a positive root the fitted \(x\) is the positive Perron vector of \(S\) (Section 1.3). So \(\zeta_1=\tau>1\) (Section 2.3), \(\zeta_1\) is simple, and the other eigenvalues are strictly smaller. The loading-dependent part of (1.1), \(\log(1+|x|^2)-x^TSx/(1+|x|^2)\), has Hessian \(2I-2S/\tau+4xx^T/\tau^2\) at \(x=\sqrt{\tau-1}\,q_1\). Its eigenvalues are
\[
\frac{4(\tau-1)}{\tau^2}\ \text{along}\ q_1,\qquad \frac{2(\tau-\zeta_k)}{\tau}\ \text{along}\ q_k,\ k\ge2,
\]
all positive. The loading block of the block-mean Hessian is a restriction of this form, so it is positive definite. Two standard facts finish the argument.

- *Haynsworth's inertia additivity:* a symmetric matrix with a nonsingular leading block has the inertia of that block plus that of its Schur complement.
- *Schur complement and concentration:* if \(f_x(\hat x(y),y)\equiv0\) with \(f_{xx}\) nonsingular, the Hessian of \(y\mapsto f(\hat x(y),y)\) is \(f_{yy}-f_{yx}f_{xx}^{-1}f_{xy}\).

The block-mean sector therefore has the inertia of the loading block, \((2,0,0)\), plus the inertia, at the stationary point, of the Hessian of the loading-concentrated discrepancy \(F(\psi_A,\psi_B)\) of Section 2.2. The concentrated loading is the simple leading eigenvector, which the within-block permutations fix, so it is block-uniform: concentrating over \((\alpha,\beta)\) and over all of \(l\) gives the same function. Pass to \(\psi_A=e^{s+d}\), \(\psi_B=e^{s-d}\); at a stationary point an invertible change of coordinates changes the Hessian by congruence. There \(F_{ss}=\nu_{AB}>0\) at the scale optimum, and by the second fact the Schur complement in \(d\) is \(\mathcal G''(d)\) of (R.2). The block-mean sector therefore has inertia \((4,0,0)\), \((3,1,0)\) or \((3,0,1)\) as \(\mathcal G''(d)\) is positive, negative or zero. Because \(d'(v)>0\), the sign of \(\mathcal G''(d)\) at the root is that of the second derivative of (S1).

*Conclusion.* At a simple scalar minimum the full \(2p\times2p\) Hessian is positive definite. At a simple scalar maximum, the saddle between two minima, it has exactly one negative eigenvalue. At a nondegenerate fold it has one zero eigenvalue and all others positive. Replacing \(\psi\) by \(\log\psi\) is again a congruence at stationary points and preserves these conclusions.

`checks/or3_full_hessian.py`, which `build.py --checks` runs, builds the analytic Hessian at 30 digits at all 18 positive roots of ten designs, four of them saddles. In an orthonormal basis of block means and random contrasts it confirms three things:

- the cross terms vanish;
- each contrast sector equals the matrix of (S7), to within \(4\times10^{-26}\), the accuracy of the roots;
- the stated inertia.

The script does not test the fold case.

### 2.5. Equal sizes with unequal loadings: uniqueness, preference, and one fold

Assume \(m_A=m_B=m\) and \(\lambda_A<\lambda_B\). This is the setting of Proposition 2 of the article.

**Proved. Unique global minimum.** The global minimum is unique up to the common loading sign for every \(0<\phi<1\).

For this proof use a different scale coordinate,

\[
 \psi_A=b_Ae^{s+d},\qquad \psi_B=b_Be^{s-d},\qquad
 \rho=\frac{mc}{\sqrt{b_Ab_B}}.
\]
Put

\[
 u=\frac{(m-1)a_A}{b_A},\quad
 z=\frac{(m-1)a_B}{b_B},\quad
 K=u+z,\quad\delta=u-z>0,\quad \nu_{AB}=2m-1.
\]
The smaller eigenvalue before multiplication by \(e^{-s}\) is

\[
 r(d)=\cosh d-\sqrt{\sinh^2d+\rho^2}.
\]
After profiling \(s\), the \(d\)-dependent part is

\[
 \mathcal G(d)=\nu_{AB}\log\{[u e^{-d}+z e^d+r(d)]/\nu_{AB}\}-\log r(d).
\]
Since \(r\) is even and \(u>z\), \(\mathcal G(d)<\mathcal G(-d)\) for every \(d>0\). Also \(\mathcal G'(0)<0\). Hence every global minimizer has \(d>0\).

Set \(t=\tanh d\), \(R=\sqrt{\rho^2+(1-\rho^2)t^2}\), and

\[
 a=\nu_{AB}(K+1)-1,\qquad b=\nu_{AB}-1-K.
\]
Both \(a\) and \(b\) are positive: \(0<K<2(m-1)=\nu_{AB}-1\). Differentiation reduces stationarity to

\[
 \delta=F(t;\rho)
 :=\frac{t(aR-b)}{\nu_{AB}R+t^2},\qquad -1<t<1.       \tag{S10}
\]
The derivative has the explicit positive prefactor

\[
 \mathcal G_t=\frac{\nu_{AB}R+t^2}
 {(K+1-\delta t-R)(1-t^2)R}\,[F(t;\rho)-\delta].
\]
Indeed \(K+1-\delta t-R\) is \(T/\cosh d\) and is positive. Thus the sign of \(F_t\) at a stationary root is the sign of the scalar Hessian, and the zero multiplicities agree.

For \(t>0\) wherever \(F\ge0\), \(aR-b\ge0\) and

\[
 \operatorname{sign} F_t
 =\operatorname{sign}\left\{(aR-b)(\nu_{AB}R-t^2)
   +tR'(at^2+\nu_{AB}b)\right\}>0.
\]
Here \(\nu_{AB}R>t^2\), \(R'=(1-\rho^2)t/R>0\), and all remaining factors in the second term are positive. Moreover \(F(1;\rho)=K>\delta\). Therefore \(F(t;\rho)=\delta\) has exactly one positive solution, and its scalar Hessian is strictly positive. This proves the unique global minimum and its smooth, indeed analytic, dependence on all interior parameters.

**Proved. The global minimum favors B in both senses.** It has \(\beta>\alpha\) and \(\beta/\lambda_B>\alpha/\lambda_A\), equivalently B's composite share exceeds 1/2.

As \(\phi\) approaches 0, the global minimum fits B, since for equal block sizes its larger loading gives strictly larger multi-information. Along the unique continuous branch, \(\alpha=\beta\) cannot occur: the two loading scores would force \(b_A=b_B\). Likewise, if \(\alpha/\lambda_A=\beta/\lambda_B=t>0\), the two row scores imply

\[
 \frac{L_A}{1-t^2L_A}=\frac{L_B}{1-t^2L_B},
\]
which is impossible when \(L_A\ne L_B\). (The difference of the two row scores factors as \((L_A-L_B)[(m-1)(1-t^2)-m(\phi-t^2)]\) times a positive factor; the second factor together with the sum of the row scores would force \(\phi=1\) or \(t^2>1\), both excluded.) Neither inequality can change sign. This also shows that a larger raw loading and a larger normalized composite weight, although logically different measures, agree on this global branch.

**Proved. Exact existence criterion for a weaker-block local minimum.** There is a weaker-block local minimum for some \(\phi>0\) if and only if

\[
 L_A>\frac{(m-1)L_B}{m+(m-1)L_B}.              \tag{S11}
\]
When (S11) holds, there is a unique fold \(\phi_f\in(0,1)\). The weaker local minimum and the intervening saddle exist exactly for \(0<\phi<\phi_f\), merge at \(\phi_f\), and are absent for \(\phi>\phi_f\). If (S11) fails, including equality, the global solution is the only positive stationary point for every \(\phi>0\).

*Proof.* Negative-\(d\) roots correspond by oddness to \(F(t;\rho)=-\delta\) for \(t>0\). At fixed \(t>0\),

\[
 \frac{\partial F}{\partial R}
 =\frac{t(at^2+\nu_{AB}b)}{(\nu_{AB}R+t^2)^2}>0,
\]
so \(F\) increases strictly with \(\rho\). At \(\rho=0\) and \(t>0\),

\[
 F(t;0)=\frac{at-b}{\nu_{AB}+t},
 \qquad\inf_{t>0}F(t;0)=-b/\nu_{AB}.
\]
For every \(\rho>0\) this lower bound is strict. Thus extra roots are possible if and only if \(\delta<b/\nu_{AB}\). Elementary algebra converts this into (S11). If the inequality holds, continuity gives a negative minimum below \(-\delta\) for small positive \(\rho\). At \(\phi=1\) the only positive block-uniform stationary point is the exact fit: the two row-score signs show that \(c-\alpha\beta\) cannot be either positive or negative, hence \(\alpha=\lambda_A\) and \(\beta=\lambda_B\). That point has \(d>0\). Consequently the negative-\(d\) pair is absent by \(\phi=1\). The minimum of \(F\), whenever negative, rises strictly with \(\rho\). It crosses \(-\delta\) exactly once.

There are at least two roots below this crossing because the endpoints satisfy \(F(0;\rho)=0\) and \(F(1;\rho)=K\). There cannot be more: (S2) has at most three positive roots and the global \(d>0\) root already accounts for one. At the crossing the root has even multiplicity. Its multiplicity is exactly two because the quartic also has the separate global positive root and at least one negative root. Thus the stationary root has multiplicity exactly two. Also \(\partial F/\partial\rho>0\) at this interior point, so the parameter direction is transverse. This is a nondegenerate fold. Section 2.4 makes the weaker scalar minimum a full local minimum.

The same criterion can be written

\[
 \frac{mL_A}{1-L_A}>(m-1)L_B,
\]
which is the condition stated in Proposition 2 of the article and the strict endpoint stability condition when A alone is fitted at \(\phi=0\). It is a stability criterion, not the criterion deciding which endpoint fit is global.

**Proved. The weaker-block minimum emphasizes A in both senses.** It has \(\alpha>\beta\) and \(\alpha/\lambda_A>\beta/\lambda_B\), equivalently A's composite share exceeds 1/2. Every stationary point other than the global minimum has \(d<0\) in the coordinate of this section, so \(\psi_A/\psi_B=(b_A/b_B)e^{2d}<b_A/b_B<1\), because \(L_A<L_B\). Diagonal reproduction, \(\psi_A=1-\alpha^2\) and \(\psi_B=1-\beta^2\), then gives \(\alpha>\beta\). Since \(\lambda_A<\lambda_B\), also \(\alpha/\lambda_A>\beta/\lambda_A>\beta/\lambda_B\).

**Proved. Exchangeable special case.** When \(\lambda_A=\lambda_B\), \(\delta=0\). Equation (S10) gives the symmetric root \(t=0\) and, below the threshold, two roots satisfying \(aR=b\). The threshold \(\rho_*=b/a\) simplifies to \(L/(2-L)\), hence

\[
 \phi_* = \frac{1+(m-1)L}{m(2-L)}.
\]
This recovers the exchangeable pitchfork from the same scalar geometry.

### 2.6. Unequal sizes: share, crossing, and the high-correlation limit

**Proved. Positive and well-defined composite weights.** Write \(f=\alpha/\lambda_A\) and \(g=\beta/\lambda_B\). By (S6), \(f,g<1\) and \(fg>\phi\). Hence \(f>\phi g\) and \(g>\phi f\). The weights

\[
 w_A=\frac{f-\phi g}{1-\phi^2},\qquad
 w_B=\frac{g-\phi f}{1-\phi^2}
\]
are both positive. Thus absolute values are unnecessary, and the A share is

\[
 S_A=\frac{f-\phi g}{(1-\phi)(f+g)}.           \tag{S12}
\]
It exceeds 1/2 exactly when \(f>g\). The article's followed-trait share is \(\max(S_A,1-S_A)\).

**Proved. Exact candidate for a share crossing.** Suppose \(m_A\ne m_B\) and \(S_A=1/2\) at a stationary point. Write \(\alpha/\lambda_A=\beta/\lambda_B=\sqrt t\), and set

\[
 Q=\sqrt{\frac{m_B(m_B-1)}{m_A(m_A-1)}},\qquad
 r=\sqrt{\frac{(m_A-1)(m_B-1)}{m_Am_B}}.
\]
The row equations imply

\[
 \frac{L_A/(1-tL_A)}{L_B/(1-tL_B)}=Q,
\]
so necessarily

\[
 t_0=\frac{QL_B-L_A}{(Q-1)L_AL_B},\qquad
 \phi_0=(1+r)t_0-r.                            \tag{S13}
\]
If \(0<t_0<1\) and \(0<\phi_0<1\), the displayed \(\alpha\) and \(\beta\) do satisfy the row scores, so this is also sufficient for a stationary share crossing. Whether that stationary point is the global minimum must be determined by (S1) and (S2). In particular, an algebraic share crossing should not be confused with equality of two local-minimum likelihoods.

**Proved. High-correlation limit.** At \(\phi=1\) the population covariance belongs to the fitted one-factor model, with loadings \(\lambda_A\) and \(\lambda_B\). Its parametrization is locally identifiable for \(p\ge3\) and nonzero loadings: if all off-diagonal first variations vanish, \(\delta l_i/l_i+\delta l_j/l_j=0\) for every distinct pair, forcing all loading variations to vanish; then the diagonal variations force \(\delta\psi=0\). Therefore the exact-fit Hessian is positive definite. The implicit function theorem gives an analytic stationary branch near \(\phi=1\). Global minimizing covariances tend to the exact-fit covariance, since their discrepancies tend to zero and the population remains positive definite. For the positive block-uniform minimizers, their within-block covariances are \(\alpha^2\) and \(\beta^2\). Each block has at least two items, so covariance convergence gives \(\alpha\to\lambda_A\) and \(\beta\to\lambda_B\); the uniquenesses then converge as well. Thus the parameters eventually enter this locally unique branch. Thus the global minimum is unique for all \(\phi\) sufficiently close to 1.

Put

\[
 r_A=\frac{L_A}{1-L_A},\quad r_B=\frac{L_B}{1-L_B},\quad
 j_A=2(m_A-1)r_A,\quad j_B=2(m_B-1)r_B,
\]
\[
 \kappa_A=m_Ar_A,\qquad\kappa_B=m_Br_B.
\]
There are expansions

\[
 \frac\alpha{\lambda_A}=1-u(1-\phi)+O((1-\phi)^2),\quad
 \frac\beta{\lambda_B}=1-z(1-\phi)+O((1-\phi)^2),
\]
where

\[
 u=\frac{\kappa_B/j_A}{1+\kappa_B/j_A+\kappa_A/j_B},\qquad
 z=\frac{\kappa_A/j_B}{1+\kappa_B/j_A+\kappa_A/j_B}. \tag{S14}
\]
To derive (S14), insert the expansions into the row scores and compare the first-order terms:

\[
 j_Au=\kappa_B(1-u-z),\qquad
 j_Bz=\kappa_A(1-u-z).
\]
Consequently

\[
 \lim_{\phi\uparrow1}S_A=\frac{1+z-u}{2}.       \tag{S15}
\]
The limiting preference is for A if and only if

\[
 m_A(m_A-1)r_A^2>m_B(m_B-1)r_B^2,                     \tag{S16}
\]
or equivalently

\[
 \kappa_A\sqrt{\frac{m_A-1}{m_A}}
 >\kappa_B\sqrt{\frac{m_B-1}{m_B}}.
\]
Thus \(\kappa\) alone is not the exact general criterion. The weights are unidentified at \(\phi=1\) itself, because the two traits then coincide. Their limit along the unique fitted branch exists and is given by (S15). Numerical evaluation from raw normalized loadings becomes ill-conditioned as \(1/(1-\phi)\); the analytic expansion avoids that cancellation.

**Proved with computer algebra. A smooth share crossing in the unequal-size example.** For \(m_A=2\), \(\lambda_A=17/20\), \(m_B=3\) and \(\lambda_B=3/4\), the positive stationary point is unique throughout \(1/2\le\phi\le1\). Therefore the global minimum is unique there. Its share crosses 1/2 smoothly at

\[
 t_0=\frac{3088-512\sqrt3}{2601},\qquad
 \phi_0=\frac{2576}{2601}-\frac{1049\sqrt3}{7803}
 =0.75753924168399945137635637418433165\ldots.       \tag{S17}
\]
The derivative of the A share along this branch is

\[
 S_A'(\phi_0)=\frac{315514305-180683667\sqrt3}{5739392}
 =0.44621689164007189828219258607613467\ldots>0.    \tag{S18}
\]
Hence the crossing is transverse and has a positive definite full Hessian. It is not a jump between competing minima. Rounded to three decimals, the crossing is at .758, the value reported in the article; an earlier rounded estimate of .755 was low by about .00254.

The exact certificate is in `verify_unequal_blocks.py::unequal_example_certificate`. In this rational design the quartic discriminant is

\[
 \operatorname{Disc}_v(P)=
 \frac{4913(459\phi^2-689)^2}{70368744177664000000000000000000}
 H(\phi),
\]
where

\[
\begin{split}
H(\phi)={}&80519601508683987527503641\phi^{10}\\
 &+2682578090235904365337528572\phi^8\\
 &-7230857150333348166462593850\phi^6\\
 &+5608822138312392897886943748\phi^4\\
 &-1711849644321614301603840471\phi^2\\
 &+177529911547540532703125000.
\end{split}
\]
SymPy verifies this polynomial identity exactly. Its exact Sturm count gives no roots of either discriminant factor in \([1/2,1]\), and exactly one positive quartic root at \(\phi=1/2\). Since the leading and constant coefficients stay nonzero on this interval, no root enters through zero or infinity. The positive-root count is therefore one on the whole interval.

For (S18), let \(f=\alpha/\lambda_A\) and \(g=\beta/\lambda_B\). Differentiate the two normalized row scores with respect to \(\log f\), \(\log g\) and \(\phi\), evaluate at \(f=g=\sqrt{t_0}\), and solve the resulting two-by-two linear system. The script checks that the result simplifies exactly to (S18). Positivity follows from the integer comparison \(315514305^2>3(180683667)^2\). The low-\(\phi\) multi-information criterion establishes the preference for B near zero; (S17) and (S18) rigorously establish the continuous crossing on the unique branch near .758. The interval certificate below extends the comparison to \(0<\phi\le.51\).

**Proved with a computer-assisted interval certificate. Global selection in the example for \(0<\phi\le.51\).** Take \(m_A=2\), \(\lambda_A=17/20\), \(m_B=3\), \(\lambda_B=3/4\) and any \(0<\phi\le51/100\). Let \(v_3(\phi)\) be the largest positive root of (S2). Then

\[
\mathcal G(v)-\mathcal G\bigl(v_3(\phi)\bigr)\ \ge\ .0538
\qquad\text{for every other positive root } v \text{ of (S2)},
\tag{S19}
\]

and the stationary point at \(v_3(\phi)\) satisfies \(\alpha/\lambda_A<\beta/\lambda_B\), so \(S_A<1/2\) by (S12). By Section 1, a global minimum exists and every local minimum is a positive block-uniform stationary point. By Section 2.3, these points correspond one to one to the positive roots of (S2), and the discrepancy at each is (S1). Hence, for every such \(\phi\), the global minimum of the unrestricted discrepancy is unique up to a common sign and favors the longer block B. With the paragraphs above, the global minimizer of the example is unique for every \(0<\phi<1\). It favors block B for \(\phi<\phi_0\) and block A for \(\phi>\phi_0\). At \(\phi_0\), given in (S17), the change is a smooth crossing, not a jump.

The root structure is exact, although the certificate does not use it. The discriminant factor \(H\) has exactly one zero in \((0,1)\), \(\phi_c=.480237334866\ldots\), and the factor \(459\phi^2-689\) has none. Since \(P_4<0<P_0\), the number of distinct positive roots is constant on \((0,\phi_c)\) and on \((\phi_c,1)\). Exact Sturm counts give three positive roots at \(\phi=1/4\) and one at \(\phi=1/2\). Below \(\phi_c\), therefore, a nonglobal minimum favoring block A and a saddle coexist with the global minimum. The A-favoring minimum ends in a fold at \(\phi_c\). In floating point, its profile excess over the global minimum falls from \(\log(\det R_A/\det R_B)=.1614\) as \(\phi\downarrow0\) to .0674 at \(\phi_c\).

*Coordinates regular at \(\phi=0\).* As \(\phi\downarrow0\), the two minimum roots of (S2) tend to 0 and to infinity, so \(v\) cannot be used near the endpoint. Put \(\xi=cv\) and \(\eta=c/v\), so that \(\xi\eta=c^2\). (These are local symbols; \(\eta\) in Section 1 is unrelated.) Then \(U=b_A+m_B\xi\), \(V=b_B+m_A\eta\) and \(N=e_0+m_Bk_B\xi+m_Ak_A\eta\) are affine in \((\xi,\eta)\), and \(e_0\) depends on \(c\) only through \(c^2\). Define

\[
\mathcal W=\frac{N^{\nu_{AB}}}{U^{m_B-1}V^{m_A-1}},\qquad
\Pi=\xi\,V\,T_B-\eta\,U\,T_A,
\]
\[
T_B=\nu_{AB}m_Bk_BU-(m_B-1)m_BN,\qquad
T_A=\nu_{AB}m_Ak_AV-(m_A-1)m_AN.
\]

Then \(\mathcal G=\log\mathcal W\) plus terms that depend on \(\phi\) only, and

\[
\Pi(cv,c/v;c)=-\frac{c\,m_Am_B\,P(v)}{v^2}=v\,\mathcal G'(v)\,NUV.
\tag{S20}
\]

At fixed \(\phi\), the positive roots of (S2) are the zeros of \(\Pi\) on the curve \(\xi\eta=c^2\). Profile values at fixed \(\phi\) have the same order as the corresponding values of \(\mathcal W\). Two charts cover the curve: \(\xi\ge c\), where \(v\ge1\), and \(\eta\ge c\), where \(v\le1\). Since \(P_4=-c\,b_Bm_Bk_B\) and \(P_0=c\,b_Am_Ak_A\), Cauchy's root bound gives, at every positive root and uniformly down to \(c=0\),

\[
\xi\le c+\frac{\max_{j\le3}|P_j|}{b_Bm_Bk_B},\qquad
\eta\le c+\frac{\max_{j\ge1}|P_j|}{b_Am_Ak_A}.
\tag{S21}
\]

*Certificate.* The interval \([0,.51]\) is tiled by 204 subintervals \(\Phi\) of width \(1/400\) with rational endpoints. Put \(C=\lambda_A\lambda_B\Phi\). Each \(\Phi\) must pass three tests.

(i) An interval \(X\) in the \(\xi\)-chart satisfies \(h_c(X_{\rm lo})<0<h_c(X_{\rm hi})\) and \(h_c'>0\) on \(X\) for all \(c\in C\), where \(h_c(\xi)=\Pi(\xi,c^2/\xi;c)\). Thus \(X\) contains exactly one root, \(cv_3\). Along the curve, \(\mathcal G\) decreases on \(X\) before this root and increases after it. Hence \(\mathcal W\) at the root is at most \(\mathcal W\) at the midpoint of \(X\).

(ii) On \(X\times C\), \((1-\psi_B)/L_B>(1-\psi_A)/L_A\), where \(\psi_A=N/(\nu_{AB}V)\) and \(\psi_B=N/(\nu_{AB}U)\). This inequality is \(\beta/\lambda_B>\alpha/\lambda_A\).

(iii) The rest of the \(\xi\)-chart and the whole \(\eta\)-chart, up to the bounds (S21), are covered by subintervals. On each of them, either \(0\notin\Pi\), so there is no root, or the infimum of \(\mathcal W\) exceeds the supremum over \(C\) of \(\mathcal W\) at the midpoint of \(X\), so any root there has a larger profile value. Within a chart, \(\xi\), \(\eta\) and \(c\) are enclosed as independent intervals, and the relation \(\eta=c^2/\xi\le c\) (or its mirror image) is used only to bound the other coordinate. This can only widen the enclosures.

All 204 subintervals pass the three tests, using 31,474 sub-boxes in all. No root lies above \(X\), so \(v_3\) is the largest positive root. On \([.4875,.51]\) it is the only positive root, consistent with the value of \(\phi_c\). The smallest certified excess, .0538, occurs just above \(\phi_c\). The decisions use only \(+,-,\times,\div\) and integer powers. Every decision is made twice and independently: once in exact rational interval arithmetic (Python `fractions`, no rounding) and once in outward-rounded `mpmath.iv` arithmetic at 113 bits. The two runs make identical decisions. Floating-point roots serve only to place \(X\).

*The endpoint.* The first subinterval, \(\Phi=[0,1/400]\), contains \(c=0\). There \(U\), \(V\), \(N\), \(\Pi\) and \(\mathcal W\) remain polynomial or rational with positive denominators. Tests (i) to (iii) on this subinterval therefore hold simultaneously for every \(\phi\in(0,1/400]\), and no limiting argument is needed for positive \(\phi\). At \(\phi=0\) itself, Section 2.1 (2a.1) applies. Since \(\det R_B=833/2048<\det R_A=76479/160000\), we have \(\mathrm{TC}(R_B)>\mathrm{TC}(R_A)\), and the fit retains block B. Section 2a.3 gives convergence as \(\phi\downarrow0\). As an exact check, consider the roots at \(c=0\): \(\xi=8777/8400\) (block B retained) and \(\eta=3625/1776\) (block A retained). Their \(\mathcal W\) ratio is exactly \(\det R_A/\det R_B\), the limit of the excess in (S19).

*Checks and dependencies.* SymPy verifies (S20), the coefficients (S3), the derivative (S4) and the partial derivatives used in (i). It also confirms agreement with `verify_unequal_blocks.py`. A 50-digit evaluation at four values of \(\phi\) confirms that every positive root is a stationary point of the raw discrepancy, with value (S1). As a negative control, the same procedure was applied to two items at .90 versus four at .70, where two minima exchange global status at \(\phi\approx.3257\). It certifies \([.20,.21]\) and fails on every subinterval of \([.325,.35]\), as a sound procedure must. The certificate depends on four things: Section 1 (existence, interiority and block uniformity); Section 2.3 (root correspondence, admissibility and profile value); (S12); and the correctness of exact rational arithmetic and of mpmath's directed rounding. No floating-point value, multistart result or sampled grid enters a decision. The script is `checks/unequal_example_certificate.py`. Its run as `python3 -I checks/unequal_example_certificate.py --verbose` takes about 50 s with Python 3.13.16, mpmath 1.3.0 and SymPy 1.14.0, and `python3 build.py --checks` writes its output to `build/checks/unequal-example-certificate.log`.

The high-correlation limit obtained from (S15) is

\[
 S_A(1^-)=0.54935648928344790117901033438\ldots.
\]
Thus the asymptotic preference for A is modest and is governed by (S16).

**Numerical, from the exact quartic. A discontinuous switch: two items at .90 versus four at .70.** This is the second example of Section 5.2 of the article, with block A the shorter block. For every \(0<\phi<.536151\) the quartic (S2) has three positive roots: a minimum favouring A (the smallest root), a saddle, and a minimum favouring B (the largest root). At \(\phi=.536151\) the minimum favouring B meets the saddle, and above it only the minimum favouring A remains. The profile values (S1) of the two minima are equal at

\[
 \phi_M=0.325654740638\ldots .
\]
Below \(\phi_M\) the minimum favouring B is global; above it, the minimum favouring A. At \(\phi_M\) the three stationary points are

| Stationary point | \(v\) | \(\alpha\) | \(\beta\) | \(S_A\) |
|---|---|---|---|---|
| Minimum favouring A | 0.06738 | 0.8937 | 0.2763 | 0.9237 |
| Saddle | 0.4572 | 0.7197 | 0.5070 | 0.5486 |
| Minimum favouring B | 3.011 | 0.3771 | 0.6912 | 0.1027 |

As \(\phi\) increases through \(\phi_M\), the global solution therefore jumps between the two minima. Its coefficient share for A's trait jumps from .103 to .924, and its loadings jump as well. Evaluated directly at the three points, the full six-item discrepancy takes the same value, 0.9640812038, at both minima. Its Hessian in the loadings and uniquenesses is positive definite at both minima and has exactly one negative eigenvalue at the saddle. `checks/unequal_jump_example.py` computes these values at 50 digits, and `build.py --checks` runs it. They are numerical values, not a certificate. The same design is the negative control of the interval certificate above.

### 2.7. Boundaries and remaining limitations

**Proved.** All formulas above hold for arbitrary \(m_A,m_B\ge2\), including \(m_A=2\) or \(m_B=2\). The strict contrast inequality (S9) explicitly covers those cases. They remain valid for loadings arbitrarily close to 0 or 1 while the loadings remain in the open interval. They do not assert uniform conditioning in these limits. At a zero loading the quartic degree and positivity arguments degenerate; at a loading of 1 a block with at least two items makes the population correlation matrix singular, so the discrepancy defined with \(\log\det R\) is not the same finite problem.

The \(\phi=0\) endpoint requires separate treatment because \(v\) can tend to 0 or infinity and zero fitted block loadings become admissible. The equal-multi-information endpoint and its possible flat parameter representations are treated in Section 2.1. For unequal sizes, (S1) to (S4) are a complete exact rule to find and compare all local minima, but this resource does not supply a closed case partition of every possible coexistence region in the full five-parameter design space. It also does not enumerate all unrestricted saddle points.

## 3. The cusp unfolding and exact fold computation (Theorem 3)

### 3.0. Status

**Proved.** The symmetry-breaking transition has a nondegenerate cusp unfolding for every fixed integer \(m\ge2\) and every fixed \(0<\lambda<1\). The constants in the coordinate \(d=\tfrac12\log(\psi_A/\psi_B)\) are given below. The leading orders can be sharpened to
\[
d(\phi_*,\varepsilon)=D\varepsilon^{1/3}+O(\varepsilon),\qquad
\phi_*-\phi_f(\varepsilon)=E_{\mathrm{fold}}\varepsilon^{2/3}+O(\varepsilon^{4/3}).
\]
The calculation establishes the local cusp. The separate global equal-size classification establishes that this fold is the unique disappearance point of the weaker-block minimum, so together they prove both parts of Theorem 3 of the article.

**Proved with computer algebra.** The independent script `verify_unequal_blocks.py` checks eight stated rational or polynomial identities, including all three coefficients, the coordinate conversion, the scalar quartic, its triple root at the symmetric cusp, and an exact discriminant factorization for the .68/.72 design. Computer algebra is used for verification of identities that are also derived below; no numerical differentiation is used to assert the coefficient formulas.

**Numerical only.** The fold locations and tilt values below were calculated with mpmath at 80 decimal digits. The defining quartic and the independent exponential-coordinate profile were both checked. The .6/.8 design has a positive fold at \(0.0165468377168788\ldots\). Thus the weaker minimum is present on a positive interval; it does not exist only at \(\phi=0\).

**Proved (smoothness qualification).** The unsigned tilt used in the article (the followed-trait share minus one half) is not a two-sided smooth function at the symmetric point because it contains a maximum. The signed tilt is smooth. The unsigned tilt has the stated scaling on either oriented branch, with the conversion constant given below.

### 3.1. Exact scale profiling

**Proved.** Put
\[
n=m-1,\quad N=2m-1,\quad L=\lambda^2,\quad a=1-L,
\quad b=1+nL,
\]
and let the two loadings be \(\lambda_A=\lambda-\varepsilon\) and \(\lambda_B=\lambda+\varepsilon\). Parameters are restricted by
\[
|\varepsilon|<\min(\lambda,1-\lambda),\qquad 0<\phi<1.
\]
Write \(a_j=1-\lambda_j^2\), \(b_j=1+n\lambda_j^2\), and \(h=m\lambda_A\lambda_B\phi\). Define
\[
T(d)=b_Ae^{-d}+b_Be^d,
\]
\[
r(d)=\frac{T(d)-\sqrt{(b_Ae^{-d}-b_Be^d)^2+4h^2}}2,
\qquad
Q(d)=na_Ae^{-d}+na_Be^d+r(d).
\]
For \(\psi_A=e^{s+d},\psi_B=e^{s-d}\), the concentrated discrepancy in the active region is exactly
\[
F(s,d)=e^{-s}Q(d)+Ns-n\log a_A-n\log a_B-\log r(d)-N.
\]
It has the unique scale optimum
\[
\widehat s(d)=\log\frac{Q(d)}N,
\qquad F_{ss}(\widehat s(d),d)=N>0,
\]
and the exact scalar profile is
\[
\boxed{\mathcal G(d;\phi,\varepsilon)=N\log\frac{Q(d)}N
-n\log a_A-n\log a_B-\log r(d).}
\tag{C.1}
\]
The leading eigenvalue is strictly greater than every contrast eigenvalue, and \(Q/N\) is strictly smaller than the leading eigenvalue. Hence the scale optimum is in the active region. To exclude a better inactive scale, let \(S_0\) denote the whitened population matrix before multiplication by \(e^{-s}\), and let \(T_+\) be its leading eigenvalue. In the inactive region \(s\ge\log T_+\), the diagonal-model discrepancy has derivative \(p-e^{-s}\operatorname{tr}S_0>0\), because \(\operatorname{tr}S_0<pT_+\). Its minimum on that region is therefore at the active-region boundary, which is worse than the interior active optimum. Thus (C.1) is the actual scale profile, rather than a formal continuation of an inapplicable formula.

**Proved.** Since \(h>0\) in a neighborhood of the critical point, \(r\), \(Q\), and \(\mathcal G\) are real analytic there. Exchanging the two blocks gives the exact symmetry
\[
\mathcal G(d;\phi,\varepsilon)=\mathcal G(-d;\phi,-\varepsilon).
\tag{C.2}
\]
At \(\varepsilon=0\), put \(C=2na+b\). Then
\[
r=b\cosh d-\sqrt{b^2\sinh^2d+h^2},\qquad
Q=C\cosh d-\sqrt{b^2\sinh^2d+h^2},
\quad h=mL\phi.
\tag{C.3}
\]
The critical values are
\[
\phi_* = \frac{b}{m(2-L)},\quad
h_* = \frac{Lb}{2-L},\quad
Q_* = \frac{2aN}{2-L},\quad
\psi_{A*}=\psi_{B*} = \frac{2a}{2-L},\quad
\alpha_* = \beta_* = \frac{\lambda}{\sqrt{2-L}}.
\tag{C.4}
\]

### 3.2. The three exact coefficients

**Proved.** With \(\delta=\phi-\phi_*\), the expansion is
\[
\mathcal G(d;\phi_*+\delta,\varepsilon)
=\mathcal G(0;\phi_*+\delta,\varepsilon)-k_1\varepsilon d
+\frac{k_2}2\delta d^2+\frac{k_4}4d^4+R(d,\delta,\varepsilon),
\tag{C.5}
\]
where
\[
\boxed{k_1=\frac{2mn\lambda(2-L)}{ab},\qquad
k_2=\frac{nm^2(2-L)^3}{aNb},\qquad
k_4=\frac{2nm(2-L)^2}{NL^2}.}
\tag{C.6}
\]
All three constants are strictly positive for \(m\ge2\) and \(0<L<1\).

#### Mixed derivative

**Proved.** At \(d=\varepsilon=0\), with \(h=mL\phi\) temporarily arbitrary, the first derivatives of both \(r\) and \(Q\) in either \(d\) or \(\varepsilon\) vanish. Moreover,
\[
r_{d\varepsilon}=2n\lambda\left(1-\frac bh\right)
=-\frac{2n\lambda(b-h)}h,
\]
\[
Q_{d\varepsilon}=-\frac{2n\lambda(b+h)}h.
\]
These follow by differentiating \(a_{A,\varepsilon}=2\lambda\), \(a_{B,\varepsilon}=-2\lambda\), \(b_{A,\varepsilon}=-2n\lambda\), \(b_{B,\varepsilon}=2n\lambda\), and \(h_\varepsilon=0\). Therefore
\[
\mathcal G_{d\varepsilon}
=N\frac{Q_{d\varepsilon}}{C-h}-\frac{r_{d\varepsilon}}{b-h}
=-\frac{4mn\lambda(nL+h)}{h(C-h)}.
\]
At \(h=h_*\), the identities \(nL+h_*=LN/(2-L)\) and \(C-h_*=2aN/(2-L)\) give \(\mathcal G_{d\varepsilon}=-k_1\). In particular, the perturbation enters at first order and has the sign favoring block B when \(\varepsilon>0\).

#### Quadratic and quartic coefficients

**Proved.** At \(d=0\), let \(r_j=\partial_d^j r\) and \(Q_j=\partial_d^j Q\). Equation (C.3) gives
\[
r_0=b-h,\quad Q_0=C-h,
\quad r_2=b-\frac{b^2}h,\quad Q_2=C-\frac{b^2}h,
\]
\[
r_4=b-\frac{4b^2}h+\frac{3b^4}{h^3},\qquad
Q_4=C-\frac{4b^2}h+\frac{3b^4}{h^3}.
\]
Thus
\[
\mathcal G_{dd}=N\frac{Q_2}{Q_0}-\frac{r_2}{r_0},
\]
\[
\mathcal G_{dddd}=N\left(\frac{Q_4}{Q_0}-3\frac{Q_2^2}{Q_0^2}\right)
-\left(\frac{r_4}{r_0}-3\frac{r_2^2}{r_0^2}\right).
\]
Substituting \(h_*\) gives \(\mathcal G_{dd}=0\), \(\mathcal G_{dd\phi}=mL\partial_h \mathcal G_{dd}=k_2\), and \(\mathcal G_{dddd}=6k_4\). The last normalization is essential: the coefficient of \(d^4\) in (C.5) is \(k_4/4\), not \(k_4\). The result agrees with Lemma 7.2 of Online Resource 1: \(k_2=2K_2\) and \(k_4=4K_4\).

**Proved with computer algebra.** In `verify_unequal_blocks.py`, exact identity checks `C1` through `C4` respectively verify \(\mathcal G_{d\varepsilon}=-k_1\), \(\mathcal G_{dd}=0\), \(\mathcal G_{dd\phi}=k_2\), and \(\mathcal G_{dddd}=6k_4\), using symbolic \(L,m,h,\lambda\). SymPy reduces each difference to the zero rational function, without numerical substitution.

### 3.3. Remainder control and the local cusp theorem

**Proved.** On a sufficiently small fixed neighborhood of the critical point, the remainder in (C.5) satisfies
\[
R=O\left(
d^6+|\delta|d^4+\delta^2d^2
+|\varepsilon||d|^3+|\varepsilon\delta d|
+\varepsilon^2d^2+|\varepsilon|^3|d|
\right).
\tag{C.7}
\]
The corresponding bounds hold after one or two derivatives with respect to \(d\), by the convergent analytic expansion. In particular,
\[
\mathcal G_d=-k_1\varepsilon+k_2\delta d+k_4d^3
+O\left(d^5+|\delta||d|^3+\delta^2|d|
+|\varepsilon|d^2+|\varepsilon\delta|
+\varepsilon^2|d|+|\varepsilon|^3\right).
\tag{C.8}
\]
To see why this list suffices, (C.2) permits a Taylor monomial \(d^i\delta^j\varepsilon^k\) only when \(i+k\) is even. Consequently the terms \(\varepsilon d^2\) and \(\varepsilon^2d\) in \(\mathcal G\) vanish identically. The term \(\varepsilon d^3\) is permitted, but it is of order \(\varepsilon^2\) under \(d=O(\varepsilon^{1/3})\), whereas the retained terms are of order \(\varepsilon^{4/3}\). Under the simultaneous cusp scaling \(\delta=O(\varepsilon^{2/3})\), all other terms displayed in (C.7) are also of higher order. This checks the omitted terms in the potential and in the stationary equation.

#### Response at the symmetric threshold

**Proved.** For \(\varepsilon>0\), define \(\rho=\varepsilon^{1/3}\) and put \(d=\rho z\), \(\delta=0\). Dividing (C.8) by \(\rho^3\) gives
\[
0=-k_1+k_4z^3+O(\rho^2).
\]
The limiting root \(z=D=(k_1/k_4)^{1/3}>0\) is simple because its derivative is \(3k_4D^2>0\). The implicit function theorem gives
\[
\boxed{d(\phi_*,\varepsilon)=D\varepsilon^{1/3}+O(\varepsilon),
\qquad
D=\left(\frac{N\lambda^5}{ab(2-L)}\right)^{1/3}.}
\tag{C.9}
\]
The zero is a local minimum since \(\mathcal G_{dd}=3k_4D^2\varepsilon^{2/3}+O(\varepsilon^{4/3})>0\). For sufficiently small positive \(\varepsilon\), it is the unique stationary point in a sufficiently small neighborhood at \(\phi=\phi_*\). This follows either from the local cusp form or by combining the dominant balance in (C.8) with uniqueness of its scaled root. The global equal-size classification identifies it as the unique global minimum.

#### Fold and its remainder

**Proved.** To locate a fold, solve \(\mathcal G_d=\mathcal G_{dd}=0\). Put \(d=\rho z\), \(\delta=\rho^2 t\), \(\varepsilon=\rho^3\). After division by \(\rho^3\) and \(\rho^2\), respectively, the limiting system is
\[
-k_1+k_2tz+k_4z^3=0,\qquad k_2t+3k_4z^2=0.
\]
Its weaker-block solution is
\[
z_f=-\left(\frac{k_1}{2k_4}\right)^{1/3},\qquad
t_f=-\frac{3k_4}{k_2}\left(\frac{k_1}{2k_4}\right)^{2/3}.
\]
The determinant of the Jacobian with respect to \((z,t)\) is
\[
-6k_2k_4z_f^2\ne0.
\]
By (C.7), both scaled equations have corrections of order \(\rho^2\). The implicit function theorem therefore gives
\[
d_f(\varepsilon)=-\left(\frac{k_1}{2k_4}\right)^{1/3}\varepsilon^{1/3}
+O(\varepsilon),
\]
\[
\boxed{\phi_*-\phi_f(\varepsilon)=E_{\mathrm{fold}}\varepsilon^{2/3}
+O(\varepsilon^{4/3}),\qquad
E_{\mathrm{fold}}=\frac{3k_4}{k_2}\left(\frac{k_1}{2k_4}\right)^{2/3}.}
\tag{C.10}
\]
Equivalently,
\[
E_{\mathrm{fold}}=\frac{3\,2^{1/3}N^{2/3}(ab)^{1/3}}
{m\lambda^{2/3}(2-L)^{5/3}}.
\]
The fold is nondegenerate because \(\mathcal G_{ddd}=6k_4d_f+O(\varepsilon)\ne0\) and \(\mathcal G_{d\phi}=k_2d_f+O(\varepsilon)\ne0\). Locally, the weaker minimum and the intervening saddle exist below this fold and annihilate at it; the positive-\(d\) minimum continues through it.

**Proved, with scope stated explicitly.** Equations (C.9) and (C.10) alone are local results. The assertion that the weaker minimum exists for every \(0<\phi<\phi_f\), with no additional folds elsewhere, uses the separate global scalar classification for equal block sizes. It must not be inferred from a Taylor expansion alone. Together the local and global results prove Theorem 3 of the article for every fixed \(m,\lambda\) and sufficiently small \(\varepsilon>0\).

**Proved.** The ratio of the fold shift to the squared response at \(\phi_*\) has the limit
\[
\boxed{\frac{\phi_*-\phi_f(\varepsilon)}{d(\phi_*,\varepsilon)^2}
\longrightarrow
\frac{3k_4}{2^{2/3}k_2}
=\frac{3\,2^{1/3}ab}{mL^2(2-L)}.}
\tag{C.11}
\]
This ratio is independent of \(\varepsilon\). It still depends on \(m\) and \(\lambda\), so calling it parameter-free without qualification would be inaccurate.

### 3.4. Conversion to the tilt

**Proved.** For positive fitted block loadings, solve the definition of \(w\) explicitly. The signed B tilt is
\[
\tau_B:=\frac{w_B-w_A}{2(w_A+w_B)}
=\frac{1+\phi}{2(1-\phi)}
\frac{\beta/\lambda_B-\alpha/\lambda_A}
{\beta/\lambda_B+\alpha/\lambda_A}.
\tag{C.12}
\]
It is analytic near the symmetric critical point, where \(w_A=w_B>0\), and the unsigned tilt is \(|\tau_B|\). At a symmetric solution \(\alpha=\beta=\ell_s(\phi)>0\), with equal population loadings, the coefficient converting a fitted loading difference into tilt is
\[
\boxed{\mathrm{tilt}
=\frac{1+\phi}{4\ell_s(\phi)(1-\phi)}|\beta-\alpha|
+o(|\beta-\alpha|).}
\tag{C.13}
\]
The signed version has the corresponding signed linear term. At the critical point, \(\ell_s(\phi_*)=\lambda/\sqrt{2-L}\), so the coefficient is
\[
C_{\alpha\beta}=\frac{(2m+1-L)\sqrt{2-L}}{4Na\lambda}.
\]

**Proved.** At the critical point the diagonal relations \(\psi_A=1-\alpha^2\) and \(\psi_B=1-\beta^2\) hold along the scale-profiled curve up to \(O(d^3)\); see the order bound below. Differentiating \(d=\tfrac12\log(\psi_A/\psi_B)\) there gives
\[
\beta-\alpha=\frac{2a}{\lambda\sqrt{2-L}}d+o(d).
\]
Therefore, on the cusp stationary branches,
\[
\tau_B=C_d d+O(d^3+|\delta d|+|\varepsilon|),\qquad
\boxed{C_d=\frac{a(1+\phi_*)}{2L(1-\phi_*)}
=\frac{2m+1-L}{2NL}>0.}
\tag{C.14}
\]
*Order bound.* Work in a neighbourhood of the critical point on which the expansions converge, with \(|d|\le d_0\). There \(Q>0\) and \(h>0\), so the scale optimum \(\widehat s(d;\delta,\varepsilon)=\log(Q(d)/N)\) of Section 3.1 is real analytic. Let \(\alpha(d;\delta,\varepsilon)\) and \(\beta(d;\delta,\varepsilon)\) be the concentrated block loadings at \(\psi_A=e^{\widehat s+d}\), \(\psi_B=e^{\widehat s-d}\). By Section 1.2 the concentrated whitened loading is \(\sqrt{\zeta_1-1}\) times the leading eigenvector of \(S\), taken with positive entries. Since \(\zeta_1>1\) is simple there, both loadings are real analytic in \((d,\delta,\varepsilon)\), and on every stationary branch they are the fitted loadings. Substituting them, with \(\lambda_A=\lambda-\varepsilon\) and \(\lambda_B=\lambda+\varepsilon\), into (C.12) defines a real-analytic function \(\vartheta(d,\delta,\varepsilon)\) that equals \(\tau_B\) on stationary branches. Both blocks have \(m\) items, so exchanging them maps \((d,\delta,\varepsilon)\) to \((-d,\delta,-\varepsilon)\), as in (C.2). The exchange also swaps \(\alpha\) with \(\beta\) and \(\lambda_A\) with \(\lambda_B\), so it reverses the sign of (C.12):
\[
\vartheta(-d,\delta,-\varepsilon)=-\vartheta(d,\delta,\varepsilon).
\]
A Taylor monomial \(d^i\delta^j\varepsilon^k\) of \(\vartheta\) can therefore occur only when \(i+k\) is odd. The terms with \(i+k=1\) are \((C_d+O(\delta))d\) and \(O(\varepsilon)\). Every term with \(i+k\ge3\) is \(O(|d|^3+d^2|\varepsilon|+|d|\varepsilon^2+|\varepsilon|^3)\), and for \(|d|\le d_0\) the last three are \(O(|\varepsilon|)\). This gives the remainder in (C.14).

The slope \(C_d=\vartheta_d(0,0,0)\) is evaluated at \(\varepsilon=\delta=0\). Along the scale-profiled curve the loading equation holds, and the derivation of \(H=\Psi^{-1}(\Sigma-R)\Psi^{-1}\) in Section 1.3 uses only that equation. So \(\partial F/\partial\psi_i=(\Sigma_{ii}-1)/\psi_i^2\) there. By the envelope theorem, the derivatives \(F_s\) and \(F_d\) of the loading-concentrated discrepancy of Section 3.1 equal the partial derivatives at the concentrated loading. Moreover \(F_s=0\) at the scale optimum, and \(F_d=\mathcal G_d\) because \(F_s(\widehat s(d),d)\equiv0\). The block residuals \(\Delta_A=\alpha^2+\psi_A-1\) and \(\Delta_B=\beta^2+\psi_B-1\) therefore satisfy
\[
\frac{m\Delta_A}{\psi_A}+\frac{m\Delta_B}{\psi_B}=F_s=0,\qquad
\frac{m\Delta_A}{\psi_A}-\frac{m\Delta_B}{\psi_B}=F_d=\mathcal G_d=O(d^3),
\]
the last by (C.8) at \(\varepsilon=\delta=0\). The diagonal relations \(\psi_A=1-\alpha^2\) and \(\psi_B=1-\beta^2\) thus hold up to \(O(d^3)\). Differentiating them gives the slope of \(\beta-\alpha\) displayed above, and (C.13) converts that slope into \(C_d\).

On the cusp branches \(d=O(\varepsilon^{1/3})\) and \(\delta=O(\varepsilon^{2/3})\), so the remainder is \(O(\varepsilon)\). It includes the change of the normalization loadings \(\lambda_A,\lambda_B\), and it is smaller than the cusp term, which is of order \(\varepsilon^{1/3}\).

**Proved.** The observable form of (C.9) is consequently
\[
\boxed{\mathrm{tilt}(\phi_*,\varepsilon)
=C_dD\varepsilon^{1/3}+O(\varepsilon).}
\tag{C.15}
\]
The coefficient \(C_dD\) is the constant \(C_{\mathrm{tilt}}\) of Theorem 3 of the article, and \(E_{\mathrm{fold}}\) in (C.10) is its fold constant.
The ratio using this observable is
\[
\frac{\phi_*-\phi_f(\varepsilon)}{\mathrm{tilt}(\phi_*,\varepsilon)^2}
\longrightarrow
\frac{12\,2^{1/3}abN^2}{m(2-L)(2m+1-L)^2}.
\tag{C.16}
\]

**Proved with computer algebra.** Check `C5` verifies the equality of the two expressions for \(C_d\).

### 3.5. Exact fold equation and the extinction of metastability

**Proved.** An explicit polynomial condition is available, without continuation in \(\phi\). For arbitrary equal-size loadings \(0<\lambda_A<\lambda_B<1\), put
\[
k_A=m-b_A=(m-1)a_A,\quad k_B=m-b_B=(m-1)a_B,
\quad e_0=m(b_A+b_B)-b_Ab_B-h^2,
\]
\[
U(v)=b_A+hv,\quad V(v)=b_B+h/v,\quad
W(v)=e_0+hk_Bv+hk_A/v.
\]
At a positive stationary point,
\[
v=\frac{\beta/\psi_B}{\alpha/\psi_A}>0,\qquad
\psi_A=\frac{W(v)}{NV(v)},\quad \psi_B=\frac{W(v)}{NU(v)}.
\]
Eliminating the remaining loading relation gives the quartic
\[
P(v)=v\{v^2V(v)[NV(v)-W(v)]-U(v)[NU(v)-W(v)]\}.
\]
Its explicit form is
\[
\boxed{\begin{aligned}
P(v)={}&-b_Bhk_Bv^4
+[b_B(Nb_B-e_0)-Nh^2]v^3\\
&+h(3m-2)(b_B-b_A)v^2
+[Nh^2-b_A(Nb_A-e_0)]v+b_Ahk_A.
\end{aligned}}
\tag{C.17}
\]
The reconstruction and the equivalence of its positive roots to stationary points are established in the scalar-classification section. In particular, this polynomial is not merely a necessary resultant with unfiltered spurious roots.

**Proved.** The fold is characterized by
\[
\boxed{P(v_f;\phi_f)=0,\qquad P_v(v_f;\phi_f)=0,
\qquad v_f>0,\quad 0<\phi_f<1,}
\tag{C.18}
\]
on the weaker-block branch. Equivalently,
\[
\operatorname{Res}_v(P,P_v)=0
\]
with the same physical-root conditions. Since all coefficients of (C.17) are explicit polynomials in \(\phi\), this is an explicit discriminant condition. A discriminant zero unaccompanied by the physical-root conditions would also count collisions of negative or complex roots and is insufficient by itself.

**Proved with computer algebra.** Check `C6` verifies the quartic expansion. Check `C7` verifies the particularly informative identity
\[
P(v;\phi_*,\lambda,\lambda)
=-\frac{naLb^2}{2-L}(v-1)^3(v+1).
\tag{C.19}
\]
Thus the three positive stationary points coalesce at the symmetric cusp. Check `C8` verifies the exact rational discriminant factorization for \(m=4,\lambda_A=17/25,\lambda_B=18/25\); the independent numerical routine solves (C.18), rather than importing roots from the reference code.

**Proved, using the global equal-size classification.** The weaker minimum exists for some \(\phi>0\) if and only if
\[
\boxed{\frac{m\lambda_A^2}{1-\lambda_A^2}>(m-1)\lambda_B^2,}
\tag{C.20}
\]
equivalently,
\[
\lambda_A^2>\frac{(m-1)\lambda_B^2}{m+(m-1)\lambda_B^2}.
\]
If this inequality holds, it exists exactly for \(0<\phi<\phi_f\), has a nondegenerate fold at \(\phi_f\), and is absent above it. At equality and when the inequality is reversed, there is no weaker-block local minimum for any \(\phi>0\). The strict inequality is the endpoint stability condition for the fit to block A: at \(\phi=0\), the fitted whitened leading eigenvalue is \(1+m\lambda_A^2/(1-\lambda_A^2)\), while the omitted block has leading eigenvalue \(1+(m-1)\lambda_B^2\).

**Proved.** Along \(\lambda_A=\lambda-\varepsilon,\lambda_B=\lambda+\varepsilon\), the extinction boundary solves
\[
m(\lambda-\varepsilon_c)^2
=(m-1)(\lambda+\varepsilon_c)^2[1-(\lambda-\varepsilon_c)^2].
\tag{C.21}
\]
The difference between the two sides of (C.20) is strictly decreasing in \(\varepsilon\), so there is at most one such boundary. If
\[
0<\lambda<\frac{1+\sqrt{(m-1)/(2m-1)}}2,
\]
the boundary lies inside the admissible interval \(0<\varepsilon<\min(\lambda,1-\lambda)\). At equality it occurs at the excluded endpoint \(\lambda_B=1\). For larger \(\lambda\), every admissible asymmetry satisfies (C.20), so there is no extinction boundary inside the parameter domain. This qualification matters when discussing loading values close to one.

### 3.6. Independent numerical values

**Numerical only.** For \(\lambda=.7,m=4\), the formulas give
\[
\begin{aligned}
\phi_*&=0.40894039735099337748\ldots,\\
k_1&=20.13812812574422481543\ldots,\\
k_2&=18.74161058755486000068\ldots,\\
k_4&=32.55929077170226691260\ldots,\\
D&=0.85201659644491631662\ldots,\\
E_{\mathrm{fold}}&=2.38340983028281959501\ldots,\\
C_d&=1.24052478134110787172\ldots,\\
C_{\mathrm{tilt}}=C_dD&=1.05694770200382476012\ldots.
\end{aligned}
\]
The limits in (C.11) and (C.16) are \(3.28323990256385665532\ldots\) and \(2.13349548700835912484\ldots\), respectively.

**Numerical only.** Exact double-root computations refine the earlier continuation values as follows. Values shown have twelve digits after the decimal point; the JSON retains 65 significant digits from an 80-digit calculation.

| \(\varepsilon\) | \(\lambda_A\) | \(\lambda_B\) | \(\phi_f\) | Reference continuation value |
|---:|---:|---:|---:|---:|
| .001 | .699 | .701 | .385205714901 | .3850 |
| .0025 | .6975 | .7025 | .365399984950 | .3650 |
| .005 | .695 | .705 | .340242180799 | .3400 |
| .01 | .69 | .71 | .301127990088 | .3010 |
| .02 | .68 | .72 | .241633835501 | .2415 |
| .04 | .66 | .74 | .155686124698 | .1555 |
| .10 | .60 | .80 | .016546837717 | approximately .02 (earlier rounded estimate) |

**Numerical only.** The extinction value for the same center and block size is
\[
\varepsilon_c=0.12068183423973835387551379885533577\ldots,
\]
with loadings \(\lambda_A=0.57931816576026164612\ldots\) and \(\lambda_B=0.82068183423973835388\ldots\). The .6/.8 design satisfies (C.20), since \(4(.6)^2/[1-(.6)^2]=2.25>3(.8)^2=1.92\), so its strictly positive fold is also consistent with the analytic existence criterion.

**Numerical only.** The tilt computations use the exact symmetric threshold, rather than a threshold rounded to five decimal places.

| \(\varepsilon\) | Independent exact-profile tilt | Leading approximation \(1.0569477020\varepsilon^{1/3}\) |
|---:|---:|---:|
| .0005 | .083600041519 | .083889994702 |
| .001 | .105113024826 | .105694770200 |
| .002 | .131997779905 | .133167065839 |
| .004 | .165423240498 | .167779989404 |
| .008 | .206619482281 | .211389540401 |
| .016 | .256620542921 | .266334131679 |

**Numerical only.** The reference tilt values agree within \(2\times10^{-8}\). Each fold-table value agrees within the stated continuation resolution of .0005. Among the independently reproduced values in `reference_section31.json`, none of these comparisons differs by more than \(10^{-3}\). The earlier rounded estimate .02 for the .6/.8 fold differs from the precise value .016546837717 by more than \(10^{-3}\); the article reports .017.

**Numerical only.** Direct checks at smaller asymmetries show convergence to the proved constants, rather than merely approximately matching exponents:

| \(\varepsilon\) | \(d/(D\varepsilon^{1/3})\) | \(\mathrm{tilt}/(C_dD\varepsilon^{1/3})\) | \((\phi_*-\phi_f)/(E_{\mathrm{fold}}\varepsilon^{2/3})\) |
|---:|---:|---:|---:|
| \(10^{-3}\) | .999268231505 | .994495987142 | .995828839354 |
| \(10^{-6}\) | .999989504445 | .999945427581 | .999961977235 |
| \(10^{-9}\) | .999999894729 | .999999454323 | .999999620144 |
| \(10^{-12}\) | .999999998947 | .999999994543 | .999999996201 |

**Numerical only.** The routine reports residuals for \(P\), \(P_v\), \(\mathcal G_d\), \(\mathcal G_{dd}\), and the original two loading score equations. The root tolerance is \(10^{-65}\), and the independently checked residual tolerance is \(10^{-55}\). Nonzero values of \(P_{vv}\), \(P_\phi\), \(\mathcal G_{ddd}\), and \(\mathcal G_{d\phi}\) are retained for every fold. These residual checks verify the numerical computation; the existence and uniqueness statements rely on the analytic classification, not on the residuals.

### 3.7. Edge cases and scope

**Proved.** The coefficient formulas include \(m=2\), since \(n=1\) and every denominator remains nonzero. No exceptional two-item flat family occurs at the cusp, because \(\phi_*>0\). Endpoint degeneracies at \(\phi=0\) are a separate matter.

**Proved.** The limits \(\lambda\to0\) and \(\lambda\to1\) are not uniform in this theorem. For each fixed interior loading the critical point and all coefficients are well defined. As \(\lambda\to0\), \(k_4\) and the coordinate-to-tilt conversion diverge. As \(\lambda\to1\), the critical point approaches \(\phi=1\), while \(k_1\) and \(k_2\) diverge and the admissible range \(\varepsilon<1-\lambda\) shrinks. These observations preclude substituting either boundary loading into the fixed-\(\lambda\) expansion.

**Proved.** The cusp theorem is an \(\varepsilon\to0\) statement. For an asymmetry of order one, the exact quartic, the physical-root reconstruction, and (C.20) remain applicable, but the truncated cusp expressions are not asserted to be accurate. For example, the leading fold approximation can become negative while the exact .6/.8 fold remains small and positive.

## 4. Mirrored heterogeneous designs

**Proved. Local stability condition.** Suppose both blocks have the same loading vector \(\lambda_1,\ldots,\lambda_m\), with corresponding items exchanged by block swap. At a positive swap-symmetric stationary fit write its fitted loadings as \((t,t)\), where \(t=(t_1,\ldots,t_m)\), and its uniquenesses as \((\psi,\psi)\), where \(\psi_i=1-t_i^2\). Define
\[
\gamma=2\sum_i\frac{t_i^2}{1-t_i^2},\qquad
W=\sum_i\frac{\lambda_i t_i}{1-t_i^2}.
\]
The symmetric stationary equations are
\[
(1+\gamma)t_i
-(1-\lambda_i^2)\frac{t_i}{1-t_i^2}
-(1+\phi)\lambda_i W=0,\qquad i=1,\ldots,m.
\tag{M.1}
\]
They follow from the general loading score before any within-block uniformity assumption. Let \(H\) be the full Hessian in \((l,\log\psi)\) coordinates and let \(B_-\) have orthonormal columns spanning antisymmetric loading and log-uniqueness perturbations. Swap symmetry makes the symmetric and antisymmetric Hessian sectors invariant. Therefore the exact antisymmetric local-stability boundary satisfies
\[
\lambda_{\min}(B_-^T H B_-)=0,
\tag{M.2}
\]
together with (M.1). Equivalently, one can use a zero determinant with the condition that the remaining antisymmetric eigenvalues are positive. A determinant zero without this spectral condition is insufficient.

For clarity, the Hessian used here is computed analytically. With \(K=\Sigma^{-1}\), \(E=KRK\), \(D=K-E\), and a covariance perturbation \(\dot\Sigma\),
\[
\dot D=-K\dot\Sigma K+K\dot\Sigma E+E\dot\Sigma K.
\]
For a parameter perturbation \((\dot l,\dot\eta)\), where \(\eta=\log\psi\),
\[
\dot\Sigma=l\dot l^T+\dot l l^T+\operatorname{diag}(\psi\dot\eta),
\]
and the Hessian action is
\[
H\binom{\dot l}{\dot\eta}
=\binom{2\dot D l+2D\dot l}
{\psi\,\operatorname{diag}(\dot D)+\psi\dot\eta\,\operatorname{diag}(D)}.
\tag{M.3}
\]
Thus (M.1) to (M.3) give an explicit finite system for the local stability threshold. This characterization does not assert that a local stability loss is necessarily the first global transition for every heterogeneous mirrored design.

**Definition.** The mean-c approximation, which substitutes an average signal-to-noise ratio into the exchangeable threshold (Appendix D of the article), uses
\[
\bar c=\frac1m\sum_{i=1}^m\frac{\lambda_i^2}{1-\lambda_i^2},
\qquad L_{\rm eff}=\frac{\bar c}{1+\bar c},\qquad
\phi_{\rm pred}=\frac{1+(m-1)L_{\rm eff}}{m(2-L_{\rm eff})}.
\]
Here \(\bar c\) is an average within-item signal-to-noise ratio, distinct from the cross-block correlation \(c\) used elsewhere in this resource.

**Numerical only.** The independent 60-digit computations give:

| Mirrored population loadings | Antisymmetric stability threshold | Mean-c prediction | Difference |
|---|---:|---:|---:|
| (.7, .7, .7, .7) | .408940397350993 | .408940397350993 | 0 |
| (.6, .65, .75, .8) | .441074899453591 | .433546682625291 | .007528216828300 |
| (.5, .6, .7, .8) | .413141495072835 | .398335732593659 | .014805762479176 |

All these values agree with the reference values in `reference_section31.json` at much better than \(10^{-3}\). The unequal-loading values were solved with mpmath, with joint stationary-equation and determinant residual below \(10^{-48}\), and with a single near-zero antisymmetric eigenvalue and all remaining antisymmetric eigenvalues positive.

**Conjecture, not established here.** The mean-c formula is a lower bound for every heterogeneous mirrored design. Two supporting designs do not prove that assertion. No universal improved closed-form approximation is established here. These questions remain open.

## 5. Independent verification and comparison with reference values

**Numerical only.** The verification command is:

```bash
python3 verify_unequal_blocks.py --output verification_results_unequal.json
```

The script is standalone apart from NumPy, SciPy, SymPy, and mpmath. It contains its reference-number baseline, so the reference program is not needed. A different reference JSON file can be compared with `--reference-json`. The reference program is not imported or executed. The randomized full-parameter optimizer uses an independently derived rank-one objective and analytic score. Dense covariance evaluation and direct derivatives are checked separately.

**Executed result: 133 checks passed, 0 failed, in 21.16 seconds.** This includes 520 unconstrained-loading, untied-uniqueness optimizations across 52 designs (10 random full-parameter starts each, default seed 20261008), 19 exact algebraic or Sturm checks, 80-digit fold and scalar calculations, 60-digit mirrored thresholds, direct 65-digit full-covariance contrast Hessian checks, and all 29 reference comparisons. The largest random-fit objective error was 1.4211e-14; the largest within-block loading or uniqueness range was 2.7012e-08. The complete recorded results, rather than rounded tables alone, are in `verification_results_unequal.json`.

The randomized experiment includes fixed cases with block sizes two, trait correlations .001 and .999, loadings .08, .12, .96, and .98, and additional independently seeded unequal designs. Optimization uses unrestricted fitted loadings and log uniquenesses in [-24,4]. It does not search the exact zero-uniqueness boundary: that boundary is handled analytically in Section 1 and through separate inward-derivative checks. The bounded numerical search is not an optimality proof.

The numerical within-block-spread check applies only after convergence to a stationary point with a positive Hessian whose smallest eigenvalue exceeds \(10^{-6}\), using score tolerance \(5\times10^{-8}\) and spread tolerance \(2\times10^{-5}\). All starts in the executed experiment satisfied these resolution conditions. This avoids interpreting unresolved nearly flat directions as mathematical counterexamples. The JSON retains every run, including optimizer status, score norm, curvature, and spread.

**Numerical only. Reference audit.** All 29 entries of `reference_section31.json` were independently compared. None differs by more than \(10^{-3}\). The largest substantive differences are the fold refinements, all below .0005, consistent with the reference continuation step. The reference maximum spread comes from a different 60-design random sample; the independently seeded sample of this resource is reported rather than presented as an exact replay of that sample.

Two earlier rounded approximations, which are not exact targets in the JSON file, differ from the refined values by more than \(10^{-3}\). The article reports the refined values, .758 and .017.

| Quantity | Earlier rounded value | Refined value | Absolute difference |
|---|---:|---:|---:|
| Equal-share crossing, 2 at .85 versus 3 at .75 | .755 | .757539241683999 | .002539241683999 |
| Weaker-minimum fold, .6 versus .8, four per block | about .02 | .016546837716879 | .003453162283121 |

**Numerical only. Other reported values.** At \(\phi=.05\) for .68/.72 with four items per block, the weak minimum has \(F=1.248023902062169\ldots\), the global minimum has \(F=.987835089106310\ldots\), and the weak minimum's smallest full Hessian eigenvalue in \((l,\log\psi)\) coordinates is \(.481811077269085\ldots\). This verifies that the weak minimum is a local minimum in the full problem. The intervening stationary point has exactly one negative full Hessian eigenvalue.

| Design, four items per block | Followed-trait share at literal \(\phi=.40894\) | Followed-trait share at \(\phi=.9\) |
|---|---:|---:|
| .68 versus .72 | .774667148809 | .572593850907 |
| .6 versus .8 | .924273880565 | .807692930642 |

The cusp tilt table uses the exact threshold, not literal .40894. This distinction is retained in the script. At \(\phi=.9\), the unequal-size example has \(S_A=.536666256981\ldots\), which rounds to .54. At \(\phi=0\), its block information values are \(\mathrm{TC}(R_A)=.369078810945\ldots\) and \(\mathrm{TC}(R_B)=.449792671996\ldots\), while \(\kappa_A=5.207207207207\ldots\) exceeds \(\kappa_B=3.857142857143\ldots\).

## 6. Remaining gaps and exact scope

1. **The main lemma has no remaining gap under the stated two-trait assumptions.** The extension to interchangeable items requires the additional condition \(c^2\le r_Ar_B\) proved sufficient in Section 1.7; interchangeability alone does not suffice.
2. **The complete positive stationary family and all local minima are classified.** The scalar polynomial does not enumerate every stationary saddle with zero or sign-changing loadings. Such a stationary point cannot be a local minimum, by the concentration and positivity arguments, but a count of every unrestricted stationary saddle is not given here.
3. **Arbitrary unequal designs have an exact selection procedure.** It consists of quartic root isolation, multiplicity detection, and explicit likelihood comparison. A closed inequality partition of every coexistence and equal-likelihood surface in the full design space is not given here.
4. **The unequal-size example is settled for \(0<\phi<1\).** Section 2.6 proves a unique positive stationary point and a transverse equal-share crossing on \([1/2,1]\). Its interval certificate shows that, for \(0<\phi\le.51\), the minimum favoring the longer block is the unique global minimum. Every other stationary value of the profile exceeds it by at least .0538. The global minimizer is therefore unique for every \(0<\phi<1\). It favors block B below \(\phi_0\) and block A above it, and it does not jump. The comparison below \(\phi=1/2\) is computer-assisted: it relies on `checks/unequal_example_certificate.py` and inherits the dependencies stated in Section 2.6. A nonglobal minimum favoring block A exists for \(\phi<\phi_c=.480237\ldots\).
5. **The cusp unfolding (Theorem 3 of the article) is proved for each fixed interior loading and sufficiently small asymmetry.** Exact quartic and extinction results also cover larger admissible asymmetries. The cusp approximation is not asserted to be uniform as \(\lambda\to0\), \(\lambda\to1\), or for an asymmetry of order one.
6. **The mirrored-design threshold computation is a local stability calculation.** A universal mean-c lower bound, an improved universal approximation, and identification of every such local threshold with the global transition remain unproved.

## 7. Summary statements

The statements below restate the main results of this resource in self-contained form, each with its status.

### A. Lemma 1: interior, block-uniform local minima

**Proved.** Let \(m_A,m_B\ge2\), \(0<\lambda_A,\lambda_B<1\), and \(0<\phi<1\). Every local minimum of the one-factor Gaussian ML discrepancy is interior and, after a common sign change, has strictly positive fitted loadings \(\alpha\) and \(\beta\) constant within the respective blocks, with uniquenesses \(1-\alpha^2\) and \(1-\beta^2\). More strongly, every strictly positive interior stationary point is block-uniform. This follows from concentration, strict positivity of the leading eigenvector, the strict inward derivative at every feasible Heywood candidate, and the within-block product and residual-sign contradiction.

### B. Scope of the interchangeable-items extension

**Proved.** The same block-uniformity result holds for constant within-block correlations \(r_A,r_B\in(0,1)\) and a constant nonzero cross-block correlation \(c\) satisfying \(c^2\le r_A r_B\). This inequality follows when items in each block share a common component and residuals have positive variances, are mutually uncorrelated, and are uncorrelated with both common components. Properly standardized equal cross-loading designs are included. Negative cross-block correlation is handled by reflecting one block. Interchangeability by itself is insufficient.

The counterexample (1.19), with its strict inward derivative \(2/165\), shows why this qualification is necessary.

### C. Equal block sizes with unequal loadings (Proposition 2 of the article)

**Proved.** With equal block sizes and \(\lambda_A<\lambda_B\), the global one-factor ML minimizer is unique, up to a common loading sign, for every \(0<\phi<1\). It satisfies \(\beta>\alpha\) and \(\beta/\lambda_B>\alpha/\lambda_A\), so the fitted factor gives the stronger block the larger composite weight. A second local minimum favoring the weaker block exists at positive correlations if and only if
\[
\frac{m\lambda_A^2}{1-\lambda_A^2}>(m-1)\lambda_B^2.
\]
When this condition holds, the weaker minimum and an intervening saddle exist exactly for \(0<\phi<\phi_f\), merge at a unique nondegenerate fold, and are absent for \(\phi>\phi_f\). At equality or when the inequality reverses, no weaker-block minimum exists for positive correlation.

### D. Exact quartic condition for positive stationary points and folds

**Proved.** Positive stationary solutions are in one-to-one correspondence with positive roots of the quartic (S2), with the reconstruction given before (S1). Every positive root is admissible. All within-block contrast directions have positive curvature, so the scalar classification gives the full local classification. The fold solves \(P=P_v=0\) at a positive double root, with the parameter-domain restrictions; a resultant zero without these restrictions is insufficient. For general unequal block sizes, global minimization reduces to comparing (S1) at the one or two minimum roots.

**Numerical only.** For four items per block with loadings .68 and .72, \(\phi_f=.241633835501\ldots\). For .6 and .8, \(\phi_f=.016546837717\ldots\). The fold table in Section 3.6 refines the earlier continuation grid values. These decimal evaluations use 80-digit arithmetic; their defining equations and existence classification are exact.

### E. The zero-correlation information criterion

**Proved.** At \(\phi=0\), let \(\mathrm{TC}(R_j)=-\tfrac12\log\det R_j\) for \(j=A,B\). Then \(F_{\min}=2\min\{\mathrm{TC}(R_A),\mathrm{TC}(R_B)\}\), and the factor fits the block with larger \(\mathrm{TC}(R_j)\) exactly while leaving the other block diagonal. Equal information gives exactly two minimizing covariances. A retained two-item block has the nonidentifiable representations \((l_1,l_2)=(t,L/t)\), \(L\le t\le1\), with uniquenesses \(1-t^2\) and \(1-L^2/t^2\); this is a parameter family of one covariance, not a flat family interpolating between the two retained blocks.

### F. Theorem 3: cusp unfolding with explicit amplitudes

**Proved.** Fix \(m\ge2\) and \(0<\lambda<1\). Let \(L=\lambda^2\), \(a=1-L\), \(b=1+(m-1)L\), \(\nu=2m-1\) (written \(N\) in Section 3), and \(\phi_*=b/[m(2-L)]\). Under loadings \(\lambda\mp\varepsilon\), the unique global minimum at \(\phi_*\) has
\[
\mathrm{tilt}
=C_{\mathrm{tilt}}\varepsilon^{1/3}+O(\varepsilon),
\qquad
C_{\mathrm{tilt}}=\frac{2m+1-L}{2\nu L}
\left[\frac{\nu\lambda^5}{ab(2-L)}\right]^{1/3}.
\]
For sufficiently small \(\varepsilon>0\), the weaker-block minimum exists exactly for \(0<\phi<\phi_f(\varepsilon)\), and
\[
\phi_*-\phi_f(\varepsilon)
=E_{\mathrm{fold}}\varepsilon^{2/3}+O(\varepsilon^{4/3}),
\qquad
E_{\mathrm{fold}}=\frac{3\,2^{1/3}\nu^{2/3}(ab)^{1/3}}
{m\lambda^{2/3}(2-L)^{5/3}}.
\]
The fold is nondegenerate. These expansions hold for each fixed interior \(\lambda\) and do not assert uniformity at the loading boundaries.

The proof uses both the local analytic expansion and the global equal-size classification of Section 2.5. The statement that the weak minimum exists throughout the interval below the fold does not follow from the Taylor expansion alone.

### G. Smoothness of the tilt

**Proved.** Define the signed B imbalance by
\[
\tau_B=\frac{1+\phi}{2(1-\phi)}
\frac{\beta/\lambda_B-\alpha/\lambda_A}
{\beta/\lambda_B+\alpha/\lambda_A}.
\]
It is smooth near symmetry, while the reported tilt is \(|\tau_B|\). At an equal-loading symmetric fit with fitted loading \(\ell_s\), the conversion from \(|\beta-\alpha|\) is \((1+\phi)/[4\ell_s(1-\phi)]\). At the critical point, the conversion from \(|d|\), where \(d=\tfrac12\log(\psi_A/\psi_B)\), is \((2m+1-L)/[2(2m-1)L]\). An amplitude ratio independent of \(\varepsilon\) may still depend on \(m\) and \(\lambda\).

### H. The unequal-size example

**Proved with computer algebra and a computer-assisted interval certificate.** For two items with loading .85 and three with loading .75, the global minimizer is unique for every \(0<\phi<1\): on \([.5,1]\) by the exact Sturm certificate and on \((0,.51]\) by the interval certificate of Section 2.6. It favors the longer block below \(\phi_0\) and the shorter block above it. The share of block A crosses one half at
\[
\phi_0=\frac{2576}{2601}-\frac{1049\sqrt3}{7803}
=.757539241683999\ldots.
\]
Its derivative there is positive, so the crossing is smooth and is not a jump between two minima. The unsigned followed-trait share and tilt have the absolute-value corner implied by their definitions. The limiting A share as \(\phi\uparrow1\) is \(.549356489283448\ldots\), although the separate trait weights are unidentified exactly at \(\phi=1\).

**Proved.** More generally, the limiting preference near \(\phi=1\) is governed by
\[
m_A(m_A-1)\left(\frac{L_A}{1-L_A}\right)^2
\quad\hbox{versus}\quad
m_B(m_B-1)\left(\frac{L_B}{1-L_B}\right)^2,
\]
with the larger quantity receiving the larger limiting share. Thus the scalar signal-to-noise measure \(\kappa=mL/(1-L)\) alone is not the exact unequal-size criterion. The complete limiting share is given in (S14) and (S15).

### I. Mirrored heterogeneous designs

**Numerical only.** In the two heterogeneous mirrored designs examined, the mean-c formula underestimates the antisymmetric local-stability threshold by .0075282168 and .0148057625. These examples do not establish a universal lower bound. The reported threshold is the solution of the symmetric stationary equations and the antisymmetric Hessian singularity condition.
