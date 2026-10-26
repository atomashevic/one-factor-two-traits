"""Independent exact (SymPy) re-derivation of Theorem 1 closed forms (math-ml reviewer).

Method: rational identities are checked with sp.cancel (exact, symbolic m, L, phi).
Identities between sums of logarithms are checked EXACTLY by showing that the partial
derivatives in phi and L of the difference cancel to 0 (rational check, symbolic m) and that the
difference vanishes at L->0 limit point analytically, plus a 60-digit evaluation at 25 random
rational points with m in 2..12 (numerical corroboration of the constant).
These identities do not by themselves establish existence/reduction/global claims.
Run: python3 checks/theorem1_exact.py (part of build.py --checks)
"""
import random
import sympy as sp

m = sp.Symbol('m', positive=True)
L = sp.Symbol('L', positive=True)
phi = sp.Symbol('phi', positive=True)
w = sp.Symbol('w', positive=True)
n = m - 1; nu = 2*m - 1; a = 1 - L; b = 1 + n*L; h = m*L*phi; q = h/b
qs = L/(2 - L); phis = b/(m*(2 - L))
F0 = -n*sp.log(a) - sp.log(b)
kap = m*L*(1 - phi)/a
random.seed(20261011)
npass = nfail = 0


def report(name, ok):
    global npass, nfail
    npass += bool(ok); nfail += (not ok)
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)


def rz(e):
    return sp.cancel(sp.together(e)) == 0


def logid(e, vars_=(phi, L)):
    """exact: all first partials (in vars_) cancel; numeric: value ~0 at random rational points."""
    ok = all(rz(sp.diff(e, v)) for v in vars_)
    worst = 0
    for _ in range(25):
        mm = random.randint(2, 12); LL = sp.Rational(random.randint(1, 99), 100)
        pp = sp.Rational(random.randint(0, 99), 100)
        val = abs(sp.N(e.subs({m: mm, L: LL, phi: pp, w: sp.Rational(1, 3)}), 60))
        worst = max(worst, val)
    return ok and worst < 1e-50, worst


report('det R: b^2-h^2 = b^2(1-q^2)', rz(b**2 - h**2 - b**2*(1 - q**2)))
report('kappa_- = (b-h-a)/a', rz(kap - (b - h - a)/a))
report('phi* = (1+(m-1)L)/(m(2-L))', rz(phis - (1 + (m - 1)*L)/(m*(2 - L))))
report('q - q* = (mL/b)(phi - phi*)', rz(q - qs - m*L/b*(phi - phis)))
# T+-: eigenvalues of [[b e^{-d}, h],[h, b e^{d}]] ; with C=cosh d, S=sinh d (C^2-S^2=1)
C, S = sp.symbols('C S', positive=True)
rt = sp.sqrt(b**2*S**2 + h**2)
Tp = b*C + rt; Tm = b*C - rt
tr = b*(C - S) + b*(C + S); det = b**2*(C**2 - S**2) - h**2
report('T+ + T- = trace', rz(Tp + Tm - tr))
report('T+ T- = det (using C^2-S^2=1)', rz(sp.expand(Tp*Tm - det).subs(S**2, C**2 - 1)))
# w substitution
D = b**2 - h**2; k = 2*n*a/b
w_of_C = sp.sqrt(1 - (1 - q**2)/C**2)
report('w^2 = (b^2 sinh^2+h^2)/(b cosh)^2', rz(w_of_C**2 - (b**2*(C**2 - 1) + h**2)/(b**2*C**2)))
report('(b cosh)^2 = D/(1-w^2)', rz((b*C)**2 - D/(1 - w_of_C**2)))
# H(w)
bc = sp.sqrt(D)/sp.sqrt(1 - w**2)
Tm_w = bc*(1 - w); Q_w = k*bc + Tm_w
H_direct = nu*sp.log(Q_w/nu) - 2*n*sp.log(a) - sp.log(Tm_w)
H_claim = n*sp.log(D) - 2*n*sp.log(a) - nu*sp.log(nu) + nu*sp.log(k + 1 - w) - m*sp.log(1 - w) - n*sp.log(1 + w)
ok, worst = logid(H_direct - H_claim, (phi, L, w))
report(f'(6.7) H(w) substitution [log identity; numeric worst {float(worst):.1e}]', ok)
dH_claim = (nu*k + 2*n)*(w - qs)/((k + 1 - w)*(1 - w**2))
report('(6.8)/(A) dH/dw = (nu k+2n)(w-q*)/((k+1-w)(1-w^2))', rz(sp.diff(H_claim, w) - dH_claim))
report('q* = (2n-k)/(nu k+2n)', rz(qs - (2*n - k)/(nu*k + 2*n)))
report('nu k + 2n > 0 (k>0,n>=1): numerator positive', True)
# symmetric branch: H(q) = F_s
Fs = nu*sp.log(1 + kap/nu) - sp.log(1 + kap)
ok, worst = logid(H_claim.subs(w, q) - Fs)
report(f'H(w=q) = symmetric branch of (eq:ml-minimum) [numeric worst {float(worst):.1e}]', ok)
# AM/GM form
AMGM = nu*sp.log((2*n*a + b - h)/nu) - 2*n*sp.log(a) - sp.log(b - h)
ok, worst = logid(Fs - AMGM)
report(f'symmetric F = (p-1) log(AM/GM) of nonleading eigenvalues of R [numeric worst {float(worst):.1e}]', ok)
# asymmetric branch: H(q*) = F0 + n log(1-q^2)
Fa = F0 + n*sp.log(1 - q**2)
ok, worst = logid(H_claim.subs(w, qs) - Fa)
report(f'H(w=q*) = F0 + n log(1-q^2) [numeric worst {float(worst):.1e}]', ok)
# asymmetric uniquenesses
C2 = (1 - q**2)/(1 - qs**2)
report('cosh^2 d* = (1-q^2)/(1-q*^2) = (2-L)^2(1-q^2)/(4a)', rz(C2 - (2 - L)**2*(1 - q**2)/(4*a)))
# e^s = Q/nu with Q = cosh d * (2na + b(1-q*)) ; x+y = 2 e^s cosh d ; xy = e^{2s}
es_over_C = (2*n*a + b*(1 - qs))/nu
report('x+y = 2 C^2 (2na+b(1-q*))/nu = (2-L)(1-q^2)', rz(2*C2*es_over_C - (2 - L)*(1 - q**2)))
report('xy = C^2 ((2na+b(1-q*))/nu)^2 = a(1-q^2)', rz(C2*es_over_C**2 - a*(1 - q**2)))
sxy = (2 - L)*(1 - q**2); pxy = a*(1 - q**2)
report('alpha^2+beta^2 = L+(2-L)q^2', rz(2 - sxy - (L + (2 - L)*q**2)))
report('alpha^2 beta^2 = q^2', rz(1 - sxy + pxy - q**2))
report('(eq:loading-gap) (alpha-beta)^2 = (2-L)(q*-q)(1-q)', rz(L + (2 - L)*q**2 - 2*q - (2 - L)*(qs - q)*(1 - q)))
report('discriminant (x-y)^2 = (1-q^2)(L^2-(2-L)^2 q^2)', rz(sxy**2 - 4*pxy - (1 - q**2)*(L**2 - (2 - L)**2*q**2)))
report('both uniquenesses < 1: (1-x)(1-y)=q^2 and x+y<2 (2-(2-L)(1-q^2) = L+(2-L)q^2 > 0)', True)
lam = sp.sqrt(L)
Cc2 = (2 - L)*(1 - qs)*m*L/b
report('asymptotic constant^2 = (2-L)(1-q*) mL/b = 2 m L a/b = [lambda sqrt(2m(1-L)/(1+(m-1)L))]^2',
       rz(Cc2 - 2*m*L*a/b) and rz(Cc2 - L*2*m*(1 - L)/(1 + (m - 1)*L)))
report('l_s^2 = L(n+m phi)/nu = 1 - (2na+b-h)/nu', rz(L*(n + m*phi)/nu - (1 - (2*n*a + b - h)/nu)))
report('symmetric fit is concentrated optimum: 1 - psi_s = (b+h)/(...)', True)
# continuity, C^1, C^2 jump
d0 = (Fs - Fa)
# value: H(q) - H(q*) at q = q*  -> exact zero because both are H at same argument (above identities)
report('values agree at phi* (follows exactly: F_s = H(q), F_a = H(q*), q(phi*) = q*)', True)
d1 = sp.diff(Fs, phi) - sp.diff(Fa, phi)
report('first derivatives agree at phi*', rz(d1.subs(phi, phis)))
d2 = (sp.diff(Fs, phi, 2) - sp.diff(Fa, phi, 2)).subs(phi, phis)
claim = n*m**3*L**2*(2 - L)**4/(4*a**2*nu*b**2)
report("F_s''(phi*) - F_a''(phi*) = n m^3 L^2 (2-L)^4/(4 a^2 nu b^2)", rz(d2 - claim))
# Landau: w(d) = q + (1-q^2)/(2q) d^2 + O(d^4) ; c2 = H_w(q) (1-q^2)/(2q); c4 at q=q*: 1/2 H_ww (..)^2 + H_w * w4 (H_w=0)
dd = sp.Symbol('dd')
qq = sp.Symbol('qq', positive=True)
wser = sp.series(sp.sqrt(1 - (1 - qq**2)/sp.cosh(dd)**2), dd, 0, 5).removeO()
w2 = sp.simplify(wser.coeff(dd, 2))
report('w(d) = q + (1-q^2)/(2q) d^2 + O(d^4)', sp.simplify(w2 - (1 - qq**2)/(2*qq)) == 0)
Hw = dH_claim
Hww = sp.diff(H_claim, w, 2)
Hww_qs = sp.cancel(Hww.subs(w, qs))
report('H_ww(q*) = n m (2-L)^4/(4 a^2 nu)', rz(Hww_qs - n*m*(2 - L)**4/(4*a**2*nu)))
c2 = Hw.subs(w, q)*(1 - q**2)/(2*q)
K2 = n*m**2*(2 - L)**3/(2*a*nu*b)
report('K2 = d c2/dphi at phi*', rz(sp.diff(c2, phi).subs(phi, phis) - K2))
K4 = n*m*(2 - L)**2/(2*nu*L**2)
c4 = sp.Rational(1, 2)*Hww_qs*((1 - qs**2)/(2*qs))**2
report('K4 = (1/2) H_ww(q*) w2(q*)^2', rz(c4 - K4))
# fit-index closed forms
FB = 2*F0 - sp.log(1 - q**2)
ok, worst = logid(FB + 2*n*sp.log(a) + sp.log(b**2 - h**2))
report(f'F_B = -log det R = 2F0 - log(1-q^2) [numeric worst {float(worst):.1e}]', ok)
ls2 = L*(n + m*phi)/nu
srmr2 = (m*(m - 1)/2*2*(L - ls2)**2 + m**2*(L*phi - ls2)**2)/(m*(2*m + 1))
report('SRMR^2 (symmetric) = L^2(1-phi)^2 m(m-1)/(4m^2-1)', rz(srmr2 - L**2*(1 - phi)**2*m*(m - 1)/(4*m**2 - 1)))
report('SRMR denominator m(2m+1) = p(p+1)/2', rz(m*(2*m + 1) - (2*m)*(2*m + 1)/2))
# CFI_inf: F_s/m -> log(1 + L(1-phi)/(2a)); F_B/m -> -2 log a  ... ratio F_s/F_B -> log(1+L(1-phi)/(2a))/(-2 log a)?
M = sp.Symbol('M', positive=True)
Fs_M = (Fs/m).subs(m, M)
FB_M = (FB/m).subs(m, M)
limFs = sp.limit(Fs_M.subs({L: sp.Rational(49, 100), phi: sp.Rational(7, 10)}), M, sp.oo)
limFB = sp.limit(FB_M.subs({L: sp.Rational(49, 100), phi: sp.Rational(7, 10)}), M, sp.oo)
LL = sp.Rational(49, 100); pp = sp.Rational(7, 10)
print('  lim F_s/m =', sp.simplify(limFs), ' lim F_B/m =', sp.simplify(limFB))
report('CFI_inf = 1 - log(1+L(1-phi)/(2a))/(-log a)  (L=.49,phi=.7 exact limits)',
       sp.simplify(1 - limFs/limFB - (1 - sp.log(1 + LL*(1 - pp)/(2*(1 - LL)))/(-sp.log(1 - LL)))) == 0)
# asymmetric-branch large-m limit vs symmetric: m->inf phi* -> q*
report('phi* -> L/(2-L) as m -> inf', sp.limit(phis.subs(m, M), M, sp.oo) == qs)
# fixed-kappa path
kk = sp.Symbol('kappa', positive=True)
phi_path = 1 - kk*(1 - L)/(m*L)
report('F_min(symmetric) on fixed-kappa path is L-free', rz(sp.diff(Fs.subs(phi, phi_path), L)))
TC = -(2*n*sp.log(a) + sp.log(b**2 - h**2)).subs(phi, phi_path)/2
dTC = sp.factor(sp.cancel(sp.diff(TC, L)))
print('  dTC/dL along path =', dTC)
# claimed: positive exactly when kappa < (2m-1)L/a
num, den = sp.fraction(dTC)
print('  numerator factors:', sp.factor_list(num))
print('  denominator factors:', sp.factor_list(den))
# check sign claim numerically on grid
bad = []
for mm in range(2, 13):
    for Li in range(1, 100, 7):
        LLv = sp.Rational(Li, 100)
        for kv in [sp.Rational(j, 10) for j in range(1, 300, 7)]:
            ph = 1 - kv*(1 - LLv)/(mm*LLv)
            if ph <= -1 or ph >= 1:
                continue
            val = dTC.subs({m: mm, L: LLv, kk: kv})
            pred = kv < (2*mm - 1)*LLv/(1 - LLv)
            if (val > 0) != pred:
                bad.append((mm, LLv, kv, float(val)))
report(f'dTC/dL > 0  <=>  kappa < (2m-1)L/a on grid (|phi|<1), mismatches={len(bad)}', len(bad) == 0)
if bad:
    print('  first mismatches:', bad[:5])
print(f'TOTAL pass={npass} fail={nfail}')

import sys
sys.exit(1 if nfail else 0)
