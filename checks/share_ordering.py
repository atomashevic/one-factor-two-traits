#!/usr/bin/env python3
"""Share ordering in the exchangeable design (Section 5.2, Appendix D).

Exact SymPy check of two identities and one inequality.  With fitted block
loadings alpha, beta, common population loading lam and trait correlation phi:
  (1) coefficient share: (w_A-w_B)/(w_A+w_B) = (1+phi)/(1-phi) * (alpha-beta)/(alpha+beta);
  (2) score share (weights Psi^{-1} l): centred share = (1+q)/(1-q) * (alpha-beta)/(alpha+beta),
      q = alpha*beta;
  (3) phi - q = phi(1-lam^2)/(1+(m-1)lam^2) > 0 when q = m lam^2 phi/(1+(m-1)lam^2).
It then reproduces the shares .786, .750 at m = 4, lam = .7, phi = .3 from Theorem 1.
The check covers these identities only; it does not verify Theorem 1.
"""
import sympy as sp
al, be, lam, phi, m = sp.symbols('alpha beta lambda varphi m', positive=True)
wA, wB = sp.symbols('w_A w_B')
sol = sp.solve([sp.Eq(al, lam*(wA+phi*wB)), sp.Eq(be, lam*(phi*wA+wB))], [wA, wB], dict=True)[0]
lhs1 = (sol[wA]-sol[wB])/(sol[wA]+sol[wB])
assert sp.simplify(lhs1 - (1+phi)/(1-phi)*(al-be)/(al+be)) == 0
print('PASS: coefficient-share identity')
sA, sB = lam*m*al/(1-al**2), lam*m*be/(1-be**2)
q = al*be
assert sp.simplify((sA-sB)/(sA+sB) - (1+q)/(1-q)*(al-be)/(al+be)) == 0
print('PASS: score-share identity')
qq = m*lam**2*phi/(1+(m-1)*lam**2)
assert sp.simplify(phi-qq - phi*(1-lam**2)/(1+(m-1)*lam**2)) == 0
print('PASS: phi - q = phi(1-lambda^2)/(1+(m-1)lambda^2)')
# Numerical reproduction from Theorem 1 at m = 4, lambda = .7, phi = .3.
M, L, P = 4, sp.Rational(49, 100), sp.Rational(3, 10)
Q = M*L*P/(1+(M-1)*L)
s2 = L+(2-L)*Q**2                      # alpha^2 + beta^2
a_ = sp.sqrt((s2+sp.sqrt(s2**2-4*Q**2))/2); b_ = Q/a_
r = (a_-b_)/(a_+b_)
coef = sp.N(sp.Rational(1, 2)+(1+P)/(1-P)*r/2, 6); score = sp.N(sp.Rational(1, 2)+(1+Q)/(1-Q)*r/2, 6)
print('coefficient share', coef, 'score share', score)
assert abs(coef-0.786) < 5e-4 and abs(score-0.750) < 5e-4 and 0.5 < score < coef
print('PASS: shares .786 and .750 reproduced; .5 < score share < coefficient share')
