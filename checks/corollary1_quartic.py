"""Numerical check of the quartic growth in Corollary 1 (Online Resource 1, Lemma 7.2).

At phi = phi*, along the symmetry-breaking path (7.1), alpha, beta = sqrt(z) +/- eps with
psi_i = 1 - l_i^2, the full ML discrepancy F is compared with the profile H(d) of Appendix A at
d(eps) = log(psi_A/psi_B)/2. Checks, at 50 digits: F(0) = H(0); F(eps) - F(0) >= H(d) - H(0);
[H(d) - H(0)]/eps^4 tends to 16 K_4 z^2/(1-z)^4 with K_4 from (7.6); and the quartic coefficient
of F exceeds that bound. Numerical corroboration only; the argument is in Online Resource 1.
Run: python3 checks/corollary1_quartic.py (part of build.py --checks)
"""
from mpmath import mp, mpf, sqrt, log, cosh, sinh, matrix, det, inverse

mp.dps = 50
failures = 0


def report(name, ok):
    global failures
    failures += not ok
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)


for m, lam in [(2, '0.7'), (4, '0.5'), (4, '0.7'), (8, '0.6'), (12, '0.3')]:
    L = mpf(lam)**2; a = 1 - L; n = m - 1; nu = 2*m - 1; b = 1 + n*L
    h = m*L*b/(m*(2 - L)); z = L/(2 - L)

    def F(eps):
        al, be = sqrt(z) + eps, sqrt(z) - eps
        pA, pB = 1 - al**2, 1 - be**2
        contrast = n*(a/pA - log(a/pA) - 1) + n*(a/pB - log(a/pB) - 1)
        M = inverse(matrix([[pA + m*al*al, m*al*be], [m*al*be, pB + m*be*be]]))*matrix([[b, h], [h, b]])
        return contrast + M[0, 0] + M[1, 1] - log(det(M)) - 2

    def H(d):
        Tm = b*cosh(d) - sqrt(b*b*sinh(d)**2 + h*h)
        return nu*log((2*n*a*cosh(d) + Tm)/nu) - 2*n*log(a) - log(Tm)

    bound = 16*(n*m*(2 - L)**2/(2*nu*L*L))*z*z/(1 - z)**4
    eps = mpf('1e-4')
    d = log((1 - (sqrt(z) + eps)**2)/(1 - (sqrt(z) - eps)**2))/2
    dF, dH = F(eps) - F(0), H(d) - H(0)
    tag = f'm={m}, lambda={lam}'
    report(f'{tag}: F(0) = H(0)', abs(F(0) - H(0)) < mpf('1e-40'))
    report(f'{tag}: F(eps) - F(0) >= H(d) - H(0)', dF >= dH)
    report(f'{tag}: profile coefficient {float(dH/eps**4):.4f} matches 16 K4 z^2/(1-z)^4 = {float(bound):.4f}',
           abs(dH/eps**4 - bound) < mpf('1e-4')*bound)
    report(f'{tag}: quartic coefficient of F {float(dF/eps**4):.4f} exceeds the bound', dF/eps**4 > bound)

print(f'{failures} failures')
raise SystemExit(1 if failures else 0)
