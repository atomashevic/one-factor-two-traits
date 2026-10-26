## 4. Heywood boundary points

### Lemma 4: Boundary local minima are excluded except at \(m=2,\phi=0\). **Proved.**

**Step 0: at most one zero uniqueness.** A positive-definite one-factor covariance can have at most one zero uniqueness. If \(\psi_j=\psi_k=0\) with \(j\ne k\), the principal submatrix of \(\Sigma\) on \(\{j,k\}\) equals that of the rank-one matrix \(ll^T\), which is singular.

**Step 1: coordinates on the stratum.** Fix \(j\) and consider the stratum \(\psi_j=0\). Positive definiteness gives \(\Sigma_{jj}=l_j^2>0\). Fixing the sign of \(l_j\), the map
\[
(l,\ \psi_{-j})\ \longmapsto\ (t,\beta,\psi_{-j}),\qquad
t=l_j^2,\quad \beta_i=l_i/l_j\ (i\ne j),
\]
is a diffeomorphism of the stratum onto \((0,\infty)\times\mathbb R^{p-1}\times(0,\infty)^{p-1}\), with the other uniquenesses positive by Step 0.

**Step 2: factorization of the stratum covariance.** Let \(D=\operatorname{diag}(d_1,\dots,d_p)\) with \(d_j=t\) and \(d_i=\psi_i\) for \(i\ne j\), and let
\[
A=I+\sum_{i\ne j}\beta_i\,e_ie_j^T .
\]
If \(\xi\sim N(0,D)\) and \(X=A\xi\), then \(X_j=\xi_j\) and \(X_i=\beta_iX_j+\xi_i\). Hence
\[
\Sigma=ADA^T,
\]
because \(\operatorname{Cov}(X_j)=t=l_j^2\), \(\operatorname{Cov}(X_i,X_j)=\beta_it=l_il_j\), and \(\operatorname{Cov}(X_i,X_k)=\beta_i\beta_kt+\delta_{ik}\psi_i=l_il_k+\delta_{ik}\psi_i\). Since \((e_ie_j^T)(e_ke_j^T)=0\) for \(i,k\ne j\), we have \(A^{-1}=I-\sum_{i\ne j}\beta_ie_ie_j^T\) and \(\det A=1\). Therefore
\[
\log\det\Sigma=\log t+\sum_{i\ne j}\log\psi_i,
\qquad
x^T\Sigma^{-1}x=\frac{x_j^2}{t}+\sum_{i\ne j}\frac{(x_i-\beta_ix_j)^2}{\psi_i}.
\tag{4.0}
\]

**Step 3: the stratum objective and its unique minimizer.** Taking the expectation of (4.0) under \(x\sim N(0,R)\) gives \(\operatorname{tr}(\Sigma^{-1}R)\). Hence, up to the constant \(-\log\det R-p\),
\[
F=\log t+\frac1t+\sum_{i\ne j}\left[\log\psi_i+\frac{1-2\beta_iR_{ij}+\beta_i^2}{\psi_i}\right].
\tag{4.1}
\]
The objective is separable. In \(t\), the derivative \(1/t-1/t^2\) vanishes only at \(t=1\). For each \(i\), the \(\beta_i\)-derivative \(2(\beta_i-R_{ij})/\psi_i\) vanishes only at \(\beta_i=R_{ij}\), and the \(\psi_i\)-derivative \(1/\psi_i-c_i/\psi_i^2\), with \(c_i=1-2\beta_iR_{ij}+\beta_i^2\ge 1-R_{ij}^2>0\), vanishes only at \(\psi_i=c_i\). So the stratum has a unique stationary point:
\[
t=1,\quad\beta_i=R_{ij},\quad\psi_i=1-R_{ij}^2 .
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
Reading the quadratic form (4.0) at \(t=1\), \(\beta_i=R_{ij}\), the \(j\)-th column of \(K\) is
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
- if \(\phi=0\) and \(m\ge3\), any pair \(i\ne k\) in \(A\setminus\{j\}\) has \(w_iw_k=L^2/(1-L^2)^2>0\) and residual covariance \(aL>0\).

In both cases increasing \(\psi_j\) strictly decreases \(F\), so \(\Sigma_\star\) is not a local minimum. Since \(\Sigma_\star\) was the only candidate on the stratum (Step 3), no boundary point is a local minimum.

Only \(m=2,\phi=0\) remains. Then \(A\setminus\{j\}\) has one item, and every weight from \(B\) is zero, so (4.3) vanishes. The star covariance fits the two-item retained block exactly and sets the other block to the identity. It is an endpoint representation of (1.8); its global optimality and completeness are proved in Section 8. \(\square\)

*Numerical check (not part of the proof).* On 400 random correlation matrices with \(p=3,\dots,8\) and random \(j\), the closed-form precision (4.0), the intermediate identity, (4.3) and the variance form (4.4) matched direct matrix computation to \(4\times10^{-14}\). One-sided finite differences agreed to within the expected \(O(\epsilon)\) error.

