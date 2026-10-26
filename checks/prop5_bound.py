"""Numerical check of Proposition 5 (Appendix E), the information-loss bound.

1. For 300 random two-trait populations R = Lambda Phi Lambda' + Psi_0 (p = 2..8, arbitrary
   loadings, trait correlation and uniquenesses), the covariance constructed in the proof,
   uniquenesses t Psi_0 and loading c Psi_0^{1/2} u_1, is in the one-factor model and its
   Kullback-Leibler divergence from N(0, R), computed with full matrices, equals the middle
   expression of the proposition; the right-hand bound holds as well.
2. For the exchangeable population, the minimum found by multistart ML fitting from random starts
   equals the middle expression when phi >= phi* and lies strictly below it when phi < phi*,
   where it equals F_min/2 from Theorem 1.
Numerical corroboration only; the proof is in Appendix E.
Run: python3 checks/prop5_bound.py (part of build.py --checks)
"""
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(20261026)
failures = 0


def report(name, ok):
    global failures
    failures += not ok
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)


def kl(R, S):
    """KL{N(0,R) || N(0,S)} = F(S)/2 with F = log det S - log det R + tr(S^-1 R) - p."""
    return 0.5*(np.linalg.slogdet(S)[1] - np.linalg.slogdet(R)[1] + np.trace(np.linalg.solve(S, R)) - len(R))


def bounds(Lam, Phi, Psi0):
    p = len(Psi0)
    w, V = np.linalg.eigh(Phi)
    Phih = V @ np.diag(np.sqrt(np.clip(w, 0, None))) @ V.T
    A = (Lam @ Phih)/np.sqrt(Psi0)[:, None]
    k2, k1 = np.linalg.eigvalsh(A.T @ A)
    k2 = max(k2, 0.0)
    middle = 0.5*((p - 1)*np.log1p(k2/(p - 1)) - np.log1p(k2))
    return k1, k2, middle, 0.5*(k2 - np.log1p(k2)), A


worst = 0.0
feasible = upper = True
for _ in range(300):
    p = int(rng.integers(2, 9))
    Lam = rng.uniform(-1, 1, (p, 2))*rng.uniform(0.2, 1.5)
    r = rng.uniform(-0.95, 0.95)
    Phi = np.array([[1, r], [r, 1]])*rng.uniform(0.3, 2)
    Psi0 = rng.uniform(0.1, 2, p)
    R = Lam @ Phi @ Lam.T + np.diag(Psi0)
    k1, k2, middle, right, A = bounds(Lam, Phi, Psi0)
    u1 = np.linalg.eigh(A @ A.T)[1][:, -1]
    t = 1 + k2/(p - 1)
    c2 = 1 + k1 - t
    loading = np.sqrt(max(c2, 0))*np.sqrt(Psi0)*u1
    S = np.outer(loading, loading) + np.diag(t*Psi0)
    feasible &= c2 >= -1e-12 and t >= 1
    upper &= middle <= right + 1e-14
    worst = max(worst, abs(kl(R, S) - middle))
report('constructed covariance is in the one-factor model (t >= 1, c^2 >= 0) in all 300 populations', feasible)
report(f'its KL equals the middle expression (max abs difference {worst:.1e})', worst < 1e-10)
report('middle expression <= (kappa_2 - log(1 + kappa_2))/2 in all 300 populations', upper)


def ml_min(R, starts=8):
    """Smallest F/2 over the one-factor model from random starts (loadings, log uniquenesses)."""
    p = len(R)

    def f(x):
        l, psi = x[:p], np.exp(x[p:])
        S = np.outer(l, l) + np.diag(psi)
        Si = np.linalg.inv(S)
        G = Si - Si @ R @ Si
        F = np.linalg.slogdet(S)[1] - np.linalg.slogdet(R)[1] + np.trace(Si @ R) - p
        return F, np.concatenate([2*G @ l, np.diag(G)*psi])

    best = np.inf
    for _ in range(starts):
        x0 = np.concatenate([rng.uniform(-1, 1, p), np.log(rng.uniform(0.2, 1, p))])
        res = minimize(f, x0, jac=True, method='L-BFGS-B', options={'ftol': 1e-15, 'gtol': 1e-11, 'maxiter': 5000})
        best = min(best, res.fun/2)
    return best


for m, lam in [(2, 0.7), (4, 0.5), (4, 0.7)]:
    L = lam*lam; a = 1 - L; n = m - 1; nu = 2*m - 1; b = 1 + n*L
    phis = b/(m*(2 - L))
    for phi in [0.5*phis, phis - 0.05, phis + 0.05, 0.9]:
        Lam = lam*np.kron(np.eye(2), np.ones((m, 1)))
        Phi = np.array([[1, phi], [phi, 1]])
        R = Lam @ Phi @ Lam.T + a*np.eye(2*m)
        middle = bounds(Lam, Phi, np.full(2*m, a))[2]
        found = ml_min(R)
        tag = f'm={m}, lambda={lam}, phi={phi:.3f} (phi*={phis:.3f})'
        if phi >= phis:
            report(f'{tag}: minimum {found:.10f} equals the bound {middle:.10f}', abs(found - middle) < 1e-8)
        else:
            q = m*L*phi/b
            Fa = (-n*np.log(a) - np.log(b) + n*np.log(1 - q*q))/2
            report(f'{tag}: minimum {found:.10f} = F_a/2 {Fa:.10f} < bound {middle:.10f}',
                   abs(found - Fa) < 1e-8 and found < middle - 1e-6)

print(f'{failures} failures')
raise SystemExit(1 if failures else 0)
