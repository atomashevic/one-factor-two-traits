#!/usr/bin/env python3
"""Identities displayed in the proof chain of Appendix D.

Exact SymPy checks (parts 1-6, 11) and numerical corroboration (parts 7-10, 12) of the
steps that Appendix D displays without derivation.  The checks cover these
identities and values only; they are not a proof of Lemma 1, Proposition 2 or
Theorem 3.

  1. Block-mean eigenvalues at uniquenesses (tU, tV): 1/t and (b_A b_B - h^2)/(tUV);
     T(d(v)) = H / sqrt(UV).
  2. The profile (D.1) in v equals the scale profile in d = (1/2) log(U/V).
  3. Proposition 2: dG/dt = (nu r + t^2) [Fcal(t; rho) - mu_minus] / {D (1-t^2) r}.
  4. Proposition 2: numerator of dFcal/dt; dFcal/dr; Fcal(1; rho) = mu_plus.
  5. Proposition 2: mu_minus < A_0/nu  <=>  m L_A/(1-L_A) > (m-1) L_B.
  6. Contrast-Hessian identity (S9) and the tilt constant C_d.
  7. Endpoint curvature x'(I - R_B/(1+kappa_A))x against the raw discrepancy.
  8. High-correlation limit of the share against the quartic.
  9. Share-crossing candidate against the exact value for 2 x .85 / 3 x .75.
 10. Fold locations .242 and .017 from the scalar equation of Proposition 2.
 11. Contrast quadratic form of Lemma D.5 against the raw second derivative, in exact
     rational arithmetic at stationary points with rational coordinates.
 12. The mean signal-to-noise substitution for mirrored designs (.434, .398).
"""
import math
import numpy as np
import sympy as sp

ok = True


def report(name, cond):
    global ok
    print(('PASS: ' if cond else 'FAIL: ') + name)
    ok = ok and bool(cond)


# ---- 1, 2 -------------------------------------------------------------
mA, mB, LA, LB, c, v, t = sp.symbols('m_A m_B L_A L_B c v t', positive=True)
aA, aB = 1 - LA, 1 - LB
bA, bB = 1 + (mA - 1)*LA, 1 + (mB - 1)*LB
kA, kB, cA, cB = (mA - 1)*aA, (mB - 1)*aB, mA*c, mB*c
h2 = mA*mB*c**2
U, V = bA + cB*v, bB + cA/v
H = mA*V + mB*U - U*V
e0 = mA*bB + mB*bA - bA*bB - cA*cB
report('H = e_0 + c_B k_B v + c_A k_A / v', sp.simplify(H - (e0 + cB*kB*v + cA*kA/v)) == 0)
tr = bA/(t*U) + bB/(t*V)
det = (bA*bB - h2)/(t**2*U*V)
report('block-mean eigenvalues 1/t and (b_A b_B - h^2)/(tUV)',
       sp.simplify(tr - (1/t + (bA*bB - h2)/(t*U*V))) == 0
       and sp.simplify(det - (1/t)*(bA*bB - h2)/(t*U*V)) == 0)
report('k_A V + k_B U + b_A b_B - h^2 = H', sp.simplify(kA*V + kB*U + bA*bB - h2 - H) == 0)

rng = np.random.default_rng(20261117)


def design(ma, mb, la, lb, phi):
    d = dict(mA=ma, mB=mb, LA=la**2, LB=lb**2, c=la*lb*phi)
    d['aA'], d['aB'] = 1 - d['LA'], 1 - d['LB']
    d['bA'], d['bB'] = 1 + (ma - 1)*d['LA'], 1 + (mb - 1)*d['LB']
    d['nu'] = ma + mb - 1
    return d


def prof_v(D, vv):
    cA_, cB_ = D['mA']*D['c'], D['mB']*D['c']
    U_, V_ = D['bA'] + cB_*vv, D['bB'] + cA_/vv
    H_ = D['mA']*V_ + D['mB']*U_ - U_*V_
    G = (D['nu']*math.log(H_) - (D['mB'] - 1)*math.log(U_) - (D['mA'] - 1)*math.log(V_)
         - D['nu']*math.log(D['nu']) - math.log(D['bA']*D['bB'] - D['mA']*D['mB']*D['c']**2)
         - (D['mA'] - 1)*math.log(D['aA']) - (D['mB'] - 1)*math.log(D['aB']))
    return G, U_, V_, H_


def prof_d(D, d):
    hh = D['c']*math.sqrt(D['mA']*D['mB'])
    x, y = D['bA']*math.exp(-d), D['bB']*math.exp(d)
    rm = (x + y - math.sqrt((x - y)**2 + 4*hh**2))/2
    T = (D['mA'] - 1)*D['aA']*math.exp(-d) + (D['mB'] - 1)*D['aB']*math.exp(d) + rm
    return (D['nu']*math.log(T/D['nu']) + (D['mA'] - D['mB'])*d - math.log(rm)
            - (D['mA'] - 1)*math.log(D['aA']) - (D['mB'] - 1)*math.log(D['aB']))


worst = 0.0
for _ in range(500):
    D = design(int(rng.integers(2, 7)), int(rng.integers(2, 7)), rng.uniform(.1, .95),
               rng.uniform(.1, .95), rng.uniform(.01, .99))
    vv = math.exp(rng.uniform(-3, 3))
    G, U_, V_, H_ = prof_v(D, vv)
    worst = max(worst, abs(G - prof_d(D, .5*math.log(U_/V_))))
report(f'(D.1) equals the scale profile at d = (1/2)log(U/V); max difference {worst:.1e}', worst < 1e-10)

# ---- 3, 4, 5 ----------------------------------------------------------
tt, rho, mp, mm, nu = sp.symbols('t rho mu_plus mu_minus nu', positive=True)
r = sp.sqrt(rho**2 + (1 - rho**2)*tt**2)
Dd = mp + 1 - mm*tt - r
G = nu*sp.log(Dd) - sp.log(1 - r) - (nu - 1)/2*sp.log(1 - tt**2)
A1, A0 = nu*(mp + 1) - 1, nu - 1 - mp
Fc = tt*(A1*r - A0)/(nu*r + tt**2)
lhs = sp.diff(G, tt)
rhs = (nu*r + tt**2)*(Fc - mm)/(Dd*(1 - tt**2)*r)
vals = [{tt: sp.Rational(3, 10), rho: sp.Rational(2, 5), mp: sp.Rational(7, 5), mm: sp.Rational(1, 3), nu: 7},
        {tt: sp.Rational(-1, 2), rho: sp.Rational(1, 7), mp: sp.Rational(9, 10), mm: sp.Rational(1, 5), nu: 5}]
exact = sp.simplify(lhs - rhs) == 0
report('Proposition 2: derivative of the profile in t = tanh(eta) ('
       + ('exact' if exact else '30-digit evaluation at two rational points') + ')',
       exact or all(abs(sp.N((lhs - rhs).subs(s), 30)) < 1e-25 for s in vals))
R_ = sp.Symbol('R', positive=True)
Fr = tt*(A1*R_ - A0)/(nu*R_ + tt**2)
dR = (1 - rho**2)*tt/R_
num = sp.together(sp.diff(Fr, tt) + sp.diff(Fr, R_)*dR)
target = ((A1*R_ - A0)*(nu*R_ - tt**2) + tt*dR*(A1*tt**2 + nu*A0))/(nu*R_ + tt**2)**2
report('Proposition 2: numerator of dFcal/dt', sp.simplify(num - target) == 0)
report('Proposition 2: dFcal/dr = t(A_1 t^2 + nu A_0)/(nu r + t^2)^2',
       sp.simplify(sp.diff(Fr, R_) - tt*(A1*tt**2 + nu*A0)/(nu*R_ + tt**2)**2) == 0)
report('Proposition 2: Fcal(1; rho) = mu_plus', sp.simplify(Fr.subs({tt: 1, R_: 1}) - mp) == 0)
m = sp.Symbol('m', positive=True)
bA_, bB_ = 1 + (m - 1)*LA, 1 + (m - 1)*LB
muA, muB = (m - 1)*(1 - LA)/bA_, (m - 1)*(1 - LB)/bB_
nu_ = 2*m - 1
gap = (nu_ - 1 - (muA + muB))/nu_ - (muA - muB)
target = 2*m*(m - 1)*(m*LA - (m - 1)*LB*(1 - LA))/(nu_*bA_*bB_)
report('Proposition 2: A_0/nu - mu_minus has the sign of m L_A - (m-1) L_B (1-L_A)',
       sp.simplify(gap - target) == 0)

# ---- 6 ----------------------------------------------------------------
al, be = sp.symbols('alpha beta', positive=True)
psA, psB = 1 - al**2, 1 - be**2
tau = 1 + mA*al**2/psA + mB*be**2/psB
report('(S9): tau psi_A^2 - a_A(1+alpha^2) = (L_A-alpha^2)(1+alpha^2) + (m_A-2)alpha^2 psi_A + m_B beta^2 psi_A^2/psi_B',
       sp.simplify(tau*psA**2 - aA*(1 + al**2)
                   - ((LA - al**2)*(1 + al**2) + (mA - 2)*al**2*psA + mB*be**2*psA**2/psB)) == 0)
L = sp.Symbol('L', positive=True)
b = 1 + (m - 1)*L
phis = b/(m*(2 - L))
report('C_d = (1/2)(1+phi*)/(1-phi*) (a/L) = (2m+1-L)/(2 nu L)',
       sp.simplify(sp.Rational(1, 2)*(1 + phis)/(1 - phis)*(1 - L)/L - (2*m + 1 - L)/(2*(2*m - 1)*L)) == 0)


# ---- numerical helpers -------------------------------------------------
def quartic_roots(D):
    cA_, cB_ = D['mA']*D['c'], D['mB']*D['c']
    kA_, kB_ = (D['mA'] - 1)*D['aA'], (D['mB'] - 1)*D['aB']
    e0_ = D['mA']*D['bB'] + D['mB']*D['bA'] - D['bA']*D['bB'] - cA_*cB_
    J, K = D['nu']*D['bA'] - e0_, D['nu']*D['bB'] - e0_
    f, g = D['nu'] - kB_, D['nu'] - kA_
    P = [-D['bB']*cB_*kB_, D['bB']*K - cA_*cB_*kB_ - cB_**2*f,
         D['bB']*cA_*g + cA_*K - D['bA']*cB_*f - cB_*J,
         cA_**2*g - D['bA']*J + cA_*cB_*kA_, D['bA']*cA_*kA_]
    return sorted(float(z.real) for z in np.roots(P) if abs(z.imag) < 1e-9 and z.real > 0)


def point(D, vv):
    G, U_, V_, H_ = prof_v(D, vv)
    pA, pB = H_/(D['nu']*V_), H_/(D['nu']*U_)
    return math.sqrt(1 - pA), math.sqrt(1 - pB), pA, pB, G


def share_A(la, lb, phi, alpha, beta):
    f, g = alpha/la, beta/lb
    return (f - phi*g)/((1 - phi)*(f + g))


def raw_F(R, l, psi):
    S = np.outer(l, l) + np.diag(psi)
    M = np.linalg.solve(S, R)
    return float(np.trace(M) - np.log(np.linalg.det(M)) - len(l))


# ---- 7 ----------------------------------------------------------------
ma, mb, la, lb = 4, 3, .6, .8
RA = (1 - la**2)*np.eye(ma) + la**2
RB = (1 - lb**2)*np.eye(mb) + lb**2
R0 = np.block([[RA, np.zeros((ma, mb))], [np.zeros((mb, ma)), RB]])
kap = ma*la**2/(1 - la**2)
x = np.array([.3, -.2, .5])
l0 = np.r_[np.full(ma, la), np.zeros(mb)]
psi0 = np.r_[np.full(ma, 1 - la**2), np.ones(mb)]
eps = 1e-3
second = (raw_F(R0, l0 + eps*np.r_[np.zeros(ma), x], psi0)
          - 2*raw_F(R0, l0, psi0) + raw_F(R0, l0 - eps*np.r_[np.zeros(ma), x], psi0))/eps**2
form = 2*float(x @ (np.eye(mb) - RB/(1 + kap)) @ x)
report(f'endpoint curvature: second difference {second:.6f}, formula {form:.6f}', abs(second - form) < 1e-5)

# ---- 8 ----------------------------------------------------------------
def limit_share_A(ma, la, mb, lb):
    rA, rB = la**2/(1 - la**2), lb**2/(1 - lb**2)
    jA, jB, kA_, kB_ = 2*(ma - 1)*rA, 2*(mb - 1)*rB, ma*rA, mb*rB
    den = 1 + kB_/jA + kA_/jB
    return (1 + (kA_/jB)/den - (kB_/jA)/den)/2


for (ma, la, mb, lb, ref) in [(4, .6, 4, .8, 1 - .7857), (2, .85, 3, .75, .5494)]:
    D = design(ma, mb, la, lb, .9999)
    roots = quartic_roots(D)
    al_, be_, *_ = point(D, roots[-1])
    num = share_A(la, lb, .9999, al_, be_)
    lim = limit_share_A(ma, la, mb, lb)
    report(f'limit share of A, {ma}x{la}/{mb}x{lb}: formula {lim:.4f}, quartic at phi=.9999 {num:.4f}',
           len(roots) == 1 and abs(num - lim) < 2e-4 and abs(lim - ref) < 1e-4)

# ---- 9 ----------------------------------------------------------------
ma, la, mb, lb = 2, .85, 3, .75
Qm = math.sqrt(mb*(mb - 1)/(ma*(ma - 1)))
rm = math.sqrt((ma - 1)*(mb - 1)/(ma*mb))
t0 = (Qm*lb**2 - la**2)/((Qm - 1)*la**2*lb**2)
phi0 = (1 + rm)*t0 - rm
D = design(ma, mb, la, lb, phi0)
roots = quartic_roots(D)
al_, be_, *_ = point(D, roots[-1])
report(f'share crossing 2x.85/3x.75: phi_0 = {phi0:.6f}; share of A there {share_A(la, lb, phi0, al_, be_):.6f}',
       abs(phi0 - 0.757539241684) < 1e-9 and len(roots) == 1
       and abs(share_A(la, lb, phi0, al_, be_) - .5) < 1e-8)

# ---- 10 ---------------------------------------------------------------
def n_weak_roots(m, la, lb, phi):
    """Roots of Fcal(t; rho) = -mu_minus on 0 < t < 1 (stationary points with eta < 0)."""
    LA_, LB_ = la**2, lb**2
    bA__, bB__ = 1 + (m - 1)*LA_, 1 + (m - 1)*LB_
    muA_, muB_ = (m - 1)*(1 - LA_)/bA__, (m - 1)*(1 - LB_)/bB__
    nu__ = 2*m - 1
    A1_, A0_ = nu__*(muA_ + muB_ + 1) - 1, nu__ - 1 - (muA_ + muB_)
    rho_ = m*la*lb*phi/math.sqrt(bA__*bB__)
    ts = np.linspace(1e-9, 1 - 1e-9, 400001)
    rr = np.sqrt(rho_**2 + (1 - rho_**2)*ts**2)
    y = ts*(A1_*rr - A0_)/(nu__*rr + ts**2) + (muA_ - muB_)
    return int(np.sum(np.sign(y[1:]) != np.sign(y[:-1])))


def fold(m, la, lb):
    lo, hi = 1e-6, .999
    for _ in range(40):
        mid = (lo + hi)/2
        lo, hi = (mid, hi) if n_weak_roots(m, la, lb, mid) == 2 else (lo, mid)
    return (lo + hi)/2


f1, f2 = fold(4, .68, .72), fold(4, .6, .8)
report(f'folds from the scalar equation: .68/.72 -> {f1:.4f}, .6/.8 -> {f2:.4f}',
       abs(f1 - .2415) < 2e-3 and abs(f2 - .01655) < 2e-3)
report('no weaker-block roots when condition (5.2) fails (4 x .5 / 4 x .8, phi = .001, .01, .1)',
       all(n_weak_roots(4, .5, .8, p) == 0 for p in (.001, .01, .1)))

# ---- 11: contrast quadratic form of Lemma D.5, exact rational arithmetic ----
# A stationary point with rational coordinates: choose alpha, beta and c with
# alpha*beta > c and solve the two row equations for L_A and L_B.  The form is
# then compared with the second derivative of the raw discrepancy
# tr(S^{-1}R) - log det(S^{-1}R) - p along l + eps*x*e, psi + eps*y*e, computed
# from the covariance matrices.  Exact at these points; not a symbolic identity.
def contrast_check(ma, mb, al_, be_, c_, block):
    Rq = sp.Rational
    pA, pB = 1 - al_**2, 1 - be_**2
    LA_ = al_**2 + mb*be_*pA*(al_*be_ - c_)/((ma - 1)*al_*pB)
    LB_ = be_**2 + ma*al_*pB*(al_*be_ - c_)/((mb - 1)*be_*pA)
    if not (0 < LA_ < 1 and 0 < LB_ < 1 and c_**2 < LA_*LB_):
        return None
    p_ = ma + mb
    Rm = sp.zeros(p_, p_)
    for i in range(p_):
        for j in range(p_):
            Rm[i, j] = 1 if i == j else (LA_ if i < ma and j < ma else LB_ if i >= ma and j >= ma else c_)
    l = sp.Matrix([al_]*ma + [be_]*mb)
    psi = [pA]*ma + [pB]*mb
    S = l*l.T + sp.diag(*psi)
    K = S.inv()
    grad = (K - K*Rm*K)
    stationary = all(e == 0 for e in (grad*l)) and all(grad[i, i] == 0 for i in range(p_))
    tau = 1 + ma*al_**2/pA + mb*be_**2/pB
    # an exact unit contrast would need a square root; use e = (1,-1,0,...) with e'e = 2
    # and divide the raw second derivative by 2.
    e = sp.zeros(p_, 1)
    i0 = 0 if block == 'A' else ma
    e[i0], e[i0 + 1] = 1, -1
    lo, ps, aa = (al_, pA, 1 - LA_) if block == 'A' else (be_, pB, 1 - LB_)
    good = stationary
    for (x, y) in [(Rq(1), Rq(0)), (Rq(0), Rq(1)), (Rq(2, 3), Rq(-5, 7))]:
        d1 = x*(e*l.T + l*e.T) + y*sp.diag(*[e[i] for i in range(p_)])
        d2 = 2*x**2*e*e.T
        raw = ((K*d2).trace() - (K*d1*K*d1).trace() - (K*d2*K*Rm).trace()
               + 2*(K*d1*K*d1*K*Rm).trace())/2
        kap = x/ps - lo*y/ps**2
        form = 2*x**2/ps + y**2/ps**2 - 2*aa/tau*kap**2
        good = good and sp.simplify(raw - form) == 0
    return good


Rq = sp.Rational
cases = [(3, 2, Rq(1, 2), Rq(3, 5), Rq(1, 5), 'A'), (3, 2, Rq(1, 2), Rq(3, 5), Rq(1, 5), 'B'),
         (2, 4, Rq(2, 5), Rq(1, 2), Rq(1, 10), 'A'), (4, 3, Rq(3, 5), Rq(1, 3), Rq(1, 8), 'B'),
         (2, 2, Rq(1, 2), Rq(1, 2), Rq(1, 5), 'A')]
res = [contrast_check(*cs) for cs in cases]
report(f'contrast quadratic form of Lemma D.5 equals the raw second derivative '
       f'(exact rational arithmetic, {sum(r is not None for r in res)} stationary points)',
       all(r for r in res if r is not None) and sum(r is not None for r in res) >= 4)

# ---- 12: mean signal-to-noise substitution for mirrored designs --------
def mean_c(lams):
    cb = sum(x**2/(1 - x**2) for x in lams)/len(lams)
    Le = cb/(1 + cb)
    return (1 + (len(lams) - 1)*Le)/(len(lams)*(2 - Le))


v1, v2, v3 = mean_c([.6, .65, .75, .8]), mean_c([.5, .6, .7, .8]), mean_c([.7]*4)
report(f'mean signal-to-noise substitution: {v1:.4f}, {v2:.4f}; equal loadings {v3:.4f}',
       abs(v1 - .4335) < 1e-4 and abs(v2 - .3983) < 1e-4 and abs(v3 - .40894) < 1e-5)
raise SystemExit(0 if ok else 1)
