#!/usr/bin/env python3
r"""Interval certificate for the unequal-size example of Section 5.2.

Design: block A has m_A = 2 items loading 17/20 = .85, block B has m_B = 3 items
loading 3/4 = .75, and the trait correlation is phi.  Online Resource 3 (S3,
Section 2.6) proves the example's claims on [1/2, 1].  This script certifies
the remaining range, with an overlap:

  (C) For every phi in (0, 51/100], the global minimiser of the exact reduced
      ML profile is the largest positive root v_3(phi) of the Appendix D
      quartic.  The profile value there is strictly smaller than at every other
      positive root, whether a minimum or a saddle.  At v_3 the fitted loadings
      satisfy alpha/lambda_A < beta/lambda_B, so S_A < 1/2: the longer block B
      is favoured.  The script prints a certified lower bound for the gap
      between the profile values.

Run:  python3 -I checks/unequal_example_certificate.py [--backend iv|exact|both] [--verbose]
Exit status: 0 if every part succeeds, 1 if any part fails.

Reduced problem (Appendix D.2; Online Resource 3 ("S3"), Sections 2.2-2.3), in
the notation of the article and Online Resource 3 (m_A, m_B, nu = m_A+m_B-1,
k_A, k_B, c_A, c_B, e_0; verify_unequal_blocks.py writes n, k, q, A, B, D, C, E, N).  For each block put L = lambda^2,
a = 1-L and b = 1+(m-1)L, and put c = lambda_A lambda_B phi.  Then
    k_A = (m_A-1)a_A,  k_B = (m_B-1)a_B,  c_A = m_A c,  c_B = m_B c,
    e_0 = m_A b_B + m_B b_A - b_A b_B - c_A c_B,
    U = b_A + c_B v,   V = b_B + c_A/v,   N = e_0 + c_B k_B v + c_A k_A/v
      = m_A V + m_B U - U V,
    G(v) = nu log N - (m_B-1) log U - (m_A-1) log V + const(phi),
    P(v) = v{v^2 V(nu V - N) - U(nu U - N)}   (a quartic),
    G'(v) = -c m_A m_B P(v) / (v U (vV)(vN)),
where v = (beta/psi_B)/(alpha/psi_A).  By S3, every positive root of P is an
admissible block-uniform stationary point with psi_A = N/(nu V) and
psi_B = N/(nu U), and the ML discrepancy there equals G(v) (S3 (S1)-(S4)).
Part N checks this against the raw discrepancy.

Rescaled coordinates.  Put xi = c v and eta = c/v, so that xi eta = c^2.  Then
U = b_A + m_B xi, V = b_B + m_A eta and N = e_0(c) + m_B k_B xi + m_A k_A eta
are affine in (xi, eta) and regular at c = 0.  Define
    W(xi, eta; c)  = N^nu / (U^(m_B-1) V^(m_A-1))   [= exp(G - const(phi))],
    Pi(xi, eta; c) = xi V T_B - eta U T_A,
    T_B = nu m_B k_B U - (m_B-1) m_B N,   T_A = nu m_A k_A V - (m_A-1) m_A N.
Part E verifies with SymPy that Pi(cv, c/v; c) = -c m_A m_B P(v)/v^2 = v G'(v) N U V.
At fixed phi, the positive roots of P are therefore the zeros of Pi on the
curve xi eta = c^2, and G(v) < G(v') if and only if W < W'.  The curve is
covered by two charts:
    chart xi  (v >= 1): coordinate xi in [c, oo), with eta = c^2/xi <= c;
    chart eta (v <= 1): coordinate eta in [c, oo), with xi = c^2/eta <= c.
Because P_4 = -c b_B m_B k_B and P_0 = c b_A m_A k_A, Cauchy's bound gives, for
every positive root,
    xi <= c + max_{j<4} |P_j| / (b_B m_B k_B)   and
    eta <= c + max_{j>0} |P_j| / (b_A m_A k_A).

Certificate for one parameter box phi in [phi_lo, phi_hi], with exact rational
endpoints (c in C = lambda_A lambda_B [phi_lo, phi_hi]):
  (i)   Isolation.  An interval X in chart xi satisfies h_c(X_lo) < 0 < h_c(X_hi)
        and h_c'(X) > 0 for all c in C, where h_c(xi) = Pi(xi, c^2/xi; c).  Then
        h_c has exactly one zero xi_3(c) = c v_3 in X.  Along the curve,
        dG/dxi has the sign of h_c, so G decreases on X up to xi_3(c) and
        increases after it.  Hence W(xi_3(c)) <= W(centre of X) for every c
        in C, and this bound is used below.
  (ii)  Preference.  On X x C, (1-psi_B)/L_B - (1-psi_A)/L_A > 0, that is,
        beta/lambda_B > alpha/lambda_A.  By S3 (S12), S_A < 1/2.
  (iii) Exclusion.  Boxes cover chart xi on [c_lo, xi_max] minus X and chart
        eta on [c_lo, eta_max].  Each box satisfies one of two conditions.
        Either 0 is not in Pi(box), so the box holds no stationary point; or
        inf W(box) > sup_c W(centre of X), so any stationary point in the box
        has a strictly larger profile value than v_3.  The second test is used
        only on boxes of width at most 2^-10 of their chart range.  Boxes treat
        (xi, eta, c) as independent intervals.  Within a chart, the relations
        xi >= c and eta = c^2/xi (or the mirror pair) only bound the other
        coordinate.  This relaxation can only widen the enclosures.
The arithmetic is either exact rational interval arithmetic (backend
`exact`, Python fractions, no rounding) or mpmath.iv with outward rounding
(backend `iv`).  Only +, -, *, / and integer powers are used, and no logarithm
or square root enters any decision.  A parameter box that fails is bisected.
The run fails if the bisection depth is exhausted.

The endpoint phi -> 0.  The first parameter box contains c = 0.  U, V, N, Pi
and W are polynomial or rational in (xi, eta, c), with denominators bounded
away from zero.  Conditions (i) to (iii) on that box therefore hold for every
c in (0, c_1] at once, so no limiting argument is needed for phi > 0.  At
phi = 0 itself the quartic degenerates.  There S3 Section 2a.1 applies: the
factor retains the block with the smaller det R, here B (exact check in
Part E).  S3 Section 2a.3 gives the convergence of the global minimiser as
phi -> 0.  As an exact consistency check, the c = 0 limits of the two minimum
roots have W ratio det R_A / det R_B.  So the profile gap tends to
2{TC(R_B) - TC(R_A)} = log(det R_A / det R_B) = .1614.

Dependencies not re-proved here:
  - Lemma 1 / Online Resource 3, Section 1: every local minimum of the unrestricted discrepancy is
    interior, positive and block-uniform, and a global minimum exists (Online Resource 3, Section 1.1).
  - S3 Section 2.3: positive roots of P correspond one to one to positive
    block-uniform stationary points, and the discrepancy there equals G.
  - S3 (S12): S_A < 1/2 if and only if alpha/lambda_A < beta/lambda_B.
  - Correct exact arithmetic in Python's fractions module, or correct directed
    rounding in mpmath.iv.  The two backends run independently.
The floating-point root guesses (mpmath) only choose X.  They do not enter any
decision.

Parts of the run:
  E  exact SymPy identities: the quartic and profile derivative of Appendix D / S3,
     the rescaled identities used here, agreement with
     proofs/unequal_blocks/verify_unequal_blocks.py, and the phi = 0 comparison;
  N  50-digit check that every positive root is a stationary point of the raw
     ML discrepancy F(Sigma) = tr(Sigma^-1 R) - log det(Sigma^-1 R) - p with
     value G(v) (corroboration only);
  C  the certificate, once per backend;
  X  negative control.  The same procedure is applied to 2 items at .90 versus
     4 at .70, which has a Maxwell point at phi_M = .3257 (review F1).  It must
     certify [.20, .21] and must fail on every subinterval of [.325, .35].
"""

import argparse
import importlib.util
import math
import platform
import sys
import time
from fractions import Fraction as Q
from pathlib import Path

import mpmath
import sympy as sp
from mpmath import iv, mp
from mpmath.libmp import to_rational

# ----------------------------------------------------------------- design
EXAMPLE = (2, 3, Q(17, 20), Q(3, 4))     # m_A, m_B, lambda_A, lambda_B
CONTROL = (2, 4, Q(9, 10), Q(7, 10))     # negative control (review F1, Maxwell point)


def use_design(m_a, m_b, lam_a, lam_b):
    """Set the module-level design constants (exact rationals)."""
    global M_A, M_B, LAM_A, LAM_B, L_A, L_B, A_A, A_B, B_A, B_B, NU, K_A, K_B
    global LL, LEAD_A, LEAD_B
    M_A, M_B, LAM_A, LAM_B = m_a, m_b, Q(lam_a), Q(lam_b)
    L_A, L_B = LAM_A ** 2, LAM_B ** 2
    A_A, A_B = 1 - L_A, 1 - L_B
    B_A, B_B = 1 + (M_A - 1) * L_A, 1 + (M_B - 1) * L_B
    NU = M_A + M_B - 1
    K_A, K_B = (M_A - 1) * A_A, (M_B - 1) * A_B
    LL = LAM_A * LAM_B                      # c = LL * phi
    LEAD_A = B_B * M_B * K_B                # |P_4| / c
    LEAD_B = B_A * M_A * K_A                # |P_0| / c


use_design(*EXAMPLE)

FAILURES = []


def fail(msg):
    FAILURES.append(msg)
    print('FAIL:', msg)


# ------------------------------------------------- backend-generic formulas
class Design:
    """Design constants converted once into a given number type."""

    def __init__(self, conv):
        self.LA, self.LB = conv(L_A), conv(L_B)
        self.bA, self.bB = conv(B_A), conv(B_B)
        self.kA, self.kB = conv(K_A), conv(K_B)
        self.e00 = conv(M_A * B_B + M_B * B_A - B_A * B_B)

    def e0(self, c):
        return self.e00 - M_A * M_B * (c * c)

    def UVN(self, xi, eta, c):
        U = self.bA + M_B * xi
        V = self.bB + M_A * eta
        N = self.e0(c) + (M_B * xi) * self.kB + (M_A * eta) * self.kA
        return U, V, N

    def TB_TA(self, U, V, N):
        TB = (NU * M_B) * self.kB * U - ((M_B - 1) * M_B) * N
        TA = (NU * M_A) * self.kA * V - ((M_A - 1) * M_A) * N
        return TB, TA

    def Pi(self, xi, eta, c):
        U, V, N = self.UVN(xi, eta, c)
        TB, TA = self.TB_TA(U, V, N)
        return xi * V * TB - eta * U * TA

    def Pi_partials(self, xi, eta, c):
        U, V, N = self.UVN(xi, eta, c)
        TB, TA = self.TB_TA(U, V, N)
        Px = (V * TB + xi * V * (M_A * M_B * M_B) * self.kB - eta * M_B * TA
              + eta * U * ((M_A - 1) * M_A * M_B) * self.kB)
        Pe = (xi * M_A * TB - xi * V * ((M_B - 1) * M_B * M_A) * self.kA - U * TA
              - eta * U * (M_B * M_A * M_A) * self.kA)
        return Px, Pe

    def W(self, xi, eta, c):
        U, V, N = self.UVN(xi, eta, c)
        return ipow(N, NU) / (ipow(U, M_B - 1) * ipow(V, M_A - 1))

    def psi(self, xi, eta, c):
        U, V, N = self.UVN(xi, eta, c)
        return N / (NU * V), N / (NU * U)          # psi_A, psi_B

    def quartic(self, c):
        """Descending coefficients P_4..P_0 of S3 (S3) / App. D, new notation."""
        cA, cB = M_A * c, M_B * c
        e = self.e0(c)
        J, K = NU * self.bA - e, NU * self.bB - e
        f, g = NU - self.kB, NU - self.kA
        bA, bB, kA, kB = self.bA, self.bB, self.kA, self.kB
        return [-bB * cB * kB,
                bB * K - cA * cB * kB - cB * cB * f,
                bB * cA * g + cA * K - bA * cB * f - cB * J,
                cA * cA * g - bA * J + cB * cA * kA,
                bA * cA * kA]


def ipow(x, n):
    """Integer power by repeated multiplication (bases here are positive)."""
    r = x
    for _ in range(n - 1):
        r = r * x
    return r


# ------------------------------------------- exact rational interval class
class RI:
    """Closed interval with exact rational endpoints (no rounding at all)."""
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        lo = Q(lo)
        hi = lo if hi is None else Q(hi)
        if lo > hi:
            raise ValueError('empty interval')
        self.lo, self.hi = lo, hi

    @staticmethod
    def _c(x):
        return x if isinstance(x, RI) else RI(x)

    def __add__(s, t):
        t = RI._c(t)
        return RI(s.lo + t.lo, s.hi + t.hi)

    __radd__ = __add__

    def __neg__(s):
        return RI(-s.hi, -s.lo)

    def __sub__(s, t):
        t = RI._c(t)
        return RI(s.lo - t.hi, s.hi - t.lo)

    def __rsub__(s, t):
        return RI._c(t) - s

    def __mul__(s, t):
        t = RI._c(t)
        p = (s.lo * t.lo, s.lo * t.hi, s.hi * t.lo, s.hi * t.hi)
        return RI(min(p), max(p))

    __rmul__ = __mul__

    def __truediv__(s, t):
        t = RI._c(t)
        if t.lo <= 0 <= t.hi:
            raise ZeroDivisionError('interval divisor contains 0')
        return s * RI(1 / t.hi, 1 / t.lo)

    def __rtruediv__(s, t):
        return RI._c(t) / s


class ExactBackend:
    name = 'exact'
    conv = staticmethod(RI)

    def design(self):
        return Design(RI)

    @staticmethod
    def ival(lo, hi):
        return RI(lo, hi)

    @staticmethod
    def lo(x):
        return x.lo

    @staticmethod
    def hi(x):
        return x.hi


class IVBackend:
    name = 'iv'

    def __init__(self, prec):
        iv.prec = prec
        self.prec = prec

    def design(self):
        iv.prec = self.prec
        return Design(self.conv)

    @staticmethod
    def conv(q):
        q = Q(q)
        return iv.mpf(q.numerator) / iv.mpf(q.denominator)   # outward rounded

    def ival(self, lo, hi):
        return iv.mpf([self.conv(lo).a, self.conv(hi).b])

    @staticmethod
    def lo(x):
        return Q(*to_rational(x._mpi_[0]))

    @staticmethod
    def hi(x):
        return Q(*to_rational(x._mpi_[1]))


# ------------------------------------------------------------ Part E (exact)
def part_E(verifier):
    print('\n== Part E: exact symbolic identities (SymPy) ==')
    v, c, xi, eta, phi = sp.symbols('v c xi eta phi', positive=True)
    mA, mB, LA, LB = sp.symbols('m_A m_B L_A L_B', positive=True)
    aA, aB = 1 - LA, 1 - LB
    bA, bB = 1 + (mA - 1) * LA, 1 + (mB - 1) * LB
    nu = mA + mB - 1
    kA, kB = (mA - 1) * aA, (mB - 1) * aB
    cA, cB = mA * c, mB * c
    e0 = mA * bB + mB * bA - bA * bB - cA * cB
    U = bA + cB * v
    V = bB + cA / v
    N = e0 + cB * kB * v + cA * kA / v
    ok = {}
    ok['N = m_A V + m_B U - U V (the quantity H of Appendix D.2)'] = \
        sp.simplify(N - (mA * V + mB * U - U * V)) == 0
    P = sp.expand(sp.cancel(v * (v ** 2 * V * (nu * V - N) - U * (nu * U - N))))
    J, K = nu * bA - e0, nu * bB - e0
    f, g = nu - kB, nu - kA
    coeff = [-bB * cB * kB, bB * K - cA * cB * kB - cB ** 2 * f,
             bB * cA * g + cA * K - bA * cB * f - cB * J,
             cA ** 2 * g - bA * J + cB * cA * kA, bA * cA * kA]
    ok['P(v) = v{v^2 V(nu V-N) - U(nu U-N)} has the coefficients of S3 (S3) / App. D'] = \
        sp.expand(P - sum(z * v ** (4 - i) for i, z in enumerate(coeff))) == 0
    G = nu * sp.log(N) - (mB - 1) * sp.log(U) - (mA - 1) * sp.log(V)
    ok["G'(v) = -c m_A m_B P(v) / (v U (vV)(vN))"] = \
        sp.cancel(sp.diff(G, v) + c * mA * mB * P / (v * U * (v * V) * (v * N))) == 0
    # rescaled coordinates xi = c v, eta = c/v
    Ur, Vr = bA + mB * xi, bB + mA * eta
    Nr = e0 + mB * kB * xi + mA * kA * eta
    TB = nu * mB * kB * Ur - (mB - 1) * mB * Nr
    TA = nu * mA * kA * Vr - (mA - 1) * mA * Nr
    Pi = xi * Vr * TB - eta * Ur * TA
    sub = {xi: c * v, eta: c / v}
    ok['U, V, N at (xi, eta) = (cv, c/v) equal U, V, N of v'] = all(
        sp.simplify(x.subs(sub) - y) == 0 for x, y in ((Ur, U), (Vr, V), (Nr, N)))
    ok['Pi(cv, c/v; c) = -c m_A m_B P(v) / v^2'] = \
        sp.cancel(Pi.subs(sub) + c * mA * mB * P / v ** 2) == 0
    Fr = nu * sp.log(Nr) - (mB - 1) * sp.log(Ur) - (mA - 1) * sp.log(Vr)
    ok['Pi = N U V (xi F_xi - eta F_eta), F = log W'] = \
        sp.cancel(Pi - Nr * Ur * Vr * (xi * sp.diff(Fr, xi) - eta * sp.diff(Fr, eta))) == 0
    # partial derivatives used for the monotonicity test
    Px = (Vr * TB + xi * Vr * mA * mB ** 2 * kB - eta * mB * TA
          + eta * Ur * (mA - 1) * mA * mB * kB)
    Pe = (xi * mA * TB - xi * Vr * (mB - 1) * mB * mA * kA - Ur * TA
          - eta * Ur * mB * mA ** 2 * kA)
    ok['coded Pi_xi, Pi_eta are the partial derivatives of Pi'] = (
        sp.expand(Px - sp.diff(Pi, xi)) == 0 and sp.expand(Pe - sp.diff(Pi, eta)) == 0)
    # the Design class used by the certificate agrees with the symbolic forms
    Ds = Design(sp.Rational)
    subs_ex = {mA: M_A, mB: M_B, LA: sp.Rational(L_A), LB: sp.Rational(L_B)}
    ok['Design.Pi, Pi_partials, W, quartic agree with the symbolic forms'] = (
        sp.expand(Ds.Pi(xi, eta, c) - Pi.subs(subs_ex)) == 0
        and sp.expand(Ds.Pi_partials(xi, eta, c)[0] - Px.subs(subs_ex)) == 0
        and sp.expand(Ds.Pi_partials(xi, eta, c)[1] - Pe.subs(subs_ex)) == 0
        and sp.cancel(Ds.W(xi, eta, c)
                      - (Nr ** nu / (Ur ** (mB - 1) * Vr ** (mA - 1))).subs(subs_ex)) == 0
        and all(sp.expand(x - y.subs(subs_ex)) == 0
                for x, y in zip(Ds.quartic(c), coeff)))
    # cross-check with the existing S3 verifier implementation
    if verifier is not None:
        theirs = verifier.quartic_coefficients(M_A, M_B, sp.Rational(17, 20),
                                               sp.Rational(3, 4), phi)
        mine = Ds.quartic(sp.Rational(LL) * phi)
        ok['quartic equals verify_unequal_blocks.quartic_coefficients (exact, in phi)'] = all(
            sp.expand(x - y) == 0 for x, y in zip(mine, theirs))
        vv = sp.Rational(7, 3)
        th = verifier.profile_value(M_A, M_B, sp.Rational(17, 20), sp.Rational(3, 4),
                                    sp.Rational(2, 5), vv, logarithm=sp.log)
        c0 = sp.Rational(LL) * sp.Rational(2, 5)
        mineG = (sp.log(Ds.W(c0 * vv, c0 / vv, c0)) - NU * sp.log(NU)
                 - sp.log(Ds.bA * Ds.bB - M_A * M_B * c0 ** 2)
                 - (M_A - 1) * sp.log(1 - Ds.LA) - (M_B - 1) * sp.log(1 - Ds.LB))
        ok['profile equals verify_unequal_blocks.profile_value (phi = 2/5, v = 7/3)'] = \
            abs(sp.N(th - mineG, 60)) < sp.Float('1e-50')
    else:
        print('  (verify_unequal_blocks.py not found: verifier cross-check skipped)')
    # phi = 0 endpoint: S3 Section 2a.1 retains the block with smaller det R
    detA = A_A ** (M_A - 1) * B_A
    detB = A_B ** (M_B - 1) * B_B
    ok['phi = 0: det R_B < det R_A, so TC(R_B) > TC(R_A): block B retained (S3 2a.1)'] = \
        detB < detA
    # c = 0 limits of the two minimum roots and of their W ratio
    e00 = M_A * B_B + M_B * B_A - B_A * B_B
    xi3_0 = ((M_B - 1) * e00 - NU * K_B * B_A) / (M_A * M_B * K_B)
    eta1_0 = ((M_A - 1) * e00 - NU * K_A * B_B) / (M_A * M_B * K_A)
    D0 = Design(Q)
    ok['c = 0: Pi(xi3_0, 0; 0) = Pi(0, eta1_0; 0) = 0'] = (
        D0.Pi(xi3_0, Q(0), Q(0)) == 0 and D0.Pi(Q(0), eta1_0, Q(0)) == 0)
    ratio = D0.W(Q(0), eta1_0, Q(0)) / D0.W(xi3_0, Q(0), Q(0))
    ok['c = 0: W(0, eta1_0) / W(xi3_0, 0) = det R_A / det R_B exactly'] = ratio == detA / detB
    for k_, v_ in ok.items():
        print(f'  [{"ok" if v_ else "FAIL"}] {k_}')
        if not v_:
            fail('Part E: ' + k_)
    print(f'  det R_A = {detA}, det R_B = {detB}; 2{{TC(R_B) - TC(R_A)}} = '
          f'log(det R_A/det R_B) = {math.log(detA / detB):.12f}')
    print(f'  c -> 0 limits: xi_3 -> {xi3_0} = {float(xi3_0):.12f} (chart xi, block B retained); '
          f'A-favouring minimum eta -> {eta1_0} = {float(eta1_0):.12f} (chart eta)')
    # structural information (not used by the certificate): discriminant zero
    Pphi = sum(z * v ** (4 - i) for i, z in enumerate(Ds.quartic(sp.Rational(LL) * phi)))
    disc = sp.factor_list(sp.discriminant(sp.expand(Pphi), v))[1]
    crit = []
    for fac, _ in disc:
        pf = sp.Poly(fac, phi)
        if pf.degree() < 1:
            continue
        for (lo, hi), _mult in pf.intervals(eps=sp.Rational(1, 10 ** 12)):
            if hi > 0 and lo < 1:
                crit.append((pf.degree(), lo, hi))
    print('  (structure only, not used by Part C) real zeros in (0,1) of the discriminant of P: '
          + ', '.join(f'degree-{d} factor, zero in [{float(lo):.12f}, {float(hi):.12f}]'
                      for d, lo, hi in crit))
    # P_4 < 0 < P_0 on 0 < phi < 1, so the number of distinct positive roots is constant
    # between consecutive discriminant zeros; one exact Sturm count on each side.
    counts = [sp.Poly(sp.expand(Pphi.subs(phi, t)), v).count_roots(0, sp.oo)
              for t in (sp.Rational(1, 4), sp.Rational(1, 2))]
    struct_ok = len(crit) == 1 and counts == [3, 1]
    print(f'  [{"ok" if struct_ok else "FAIL"}] (structure only) exactly one discriminant zero '
          f'phi_c in (0,1); Sturm counts of distinct positive roots: {counts[0]} at phi = 1/4, '
          f'{counts[1]} at phi = 1/2; hence 3 positive roots for 0 < phi < phi_c and 1 for '
          f'phi_c < phi < 1')
    if not struct_ok:
        fail('Part E: root-count structure differs from the expected one')


# ------------------------------------------------ Part N (numerical, mpmath)
def part_N():
    print('\n== Part N: reduced profile versus the raw ML discrepancy '
          '(mpmath, 50 digits; corroboration only) ==')
    mp.dps = 50
    Dm = Design(lambda q: mp.mpf(q.numerator) / q.denominator)
    lamA, lamB = mp.mpf(17) / 20, mp.mpf(3) / 4
    p = M_A + M_B
    worst = mp.mpf(0)
    for ph in (Q(1, 100), Q(3, 10), Q(12, 25), Q(1, 2)):
        phim = mp.mpf(ph.numerator) / ph.denominator
        c = mp.mpf(LL.numerator) / LL.denominator * phim
        roots = sorted(mp.re(r) for r in mp.polyroots(Dm.quartic(c), maxsteps=200, extraprec=200)
                       if abs(mp.im(r)) < mp.mpf(10) ** -30 and mp.re(r) > 0)
        loads = [lamA] * M_A + [lamB] * M_B
        Rm = mp.matrix(p, p)
        for i in range(p):
            for j in range(p):
                same = (i < M_A) == (j < M_A)
                Rm[i, j] = 1 if i == j else loads[i] * loads[j] * (1 if same else phim)
        out = []
        for vv in roots:
            xi, eta = c * vv, c / vv
            pA, pB = Dm.psi(xi, eta, c)
            al, be = mp.sqrt(1 - pA), mp.sqrt(1 - pB)
            l = mp.matrix([al] * M_A + [be] * M_B)
            psi = [pA] * M_A + [pB] * M_B
            S = l * l.T
            for i in range(p):
                S[i, i] += psi[i]
            Si = S ** -1
            F = sum((Si * Rm)[i, i] for i in range(p)) - mp.log(mp.det(Si * Rm)) - p
            Gfull = (mp.log(Dm.W(xi, eta, c)) - NU * mp.log(NU)
                     - mp.log(Dm.bA * Dm.bB - M_A * M_B * c * c)
                     - (M_A - 1) * mp.log(1 - Dm.LA) - (M_B - 1) * mp.log(1 - Dm.LB))
            Wm = Si - Si * Rm * Si                      # dF/dSigma
            score_l = 2 * (Wm * l)
            score = max([abs(score_l[i]) for i in range(p)] + [abs(Wm[i, i]) for i in range(p)])
            worst = max(worst, abs(F - Gfull), score)
            fav = 'B' if be / lamB > al / lamA else 'A'
            out.append(f'v={mp.nstr(vv, 8)} G={mp.nstr(Gfull, 12)} favours {fav}')
        print(f'  phi={float(ph):<5}: ' + ' ; '.join(out))
    ok = worst < mp.mpf(10) ** -40
    print(f'  [{"ok" if ok else "FAIL"}] at every root: |F(Sigma) - G(v)| and the sup-norm of '
          f'the ML score are <= {mp.nstr(worst, 3)} (< 1e-40)')
    if not ok:
        fail('Part N: reduced profile does not match the ML discrepancy')
    mp.dps = 15


# ------------------------------------------------- Part C (the certificate)
def float_root_xi3(phi):
    """Floating guess of the chart-xi zero of h(xi) = Pi(xi, c^2/xi; c) (bisection)."""
    mp.dps = 30
    Dm = Design(lambda q: mp.mpf(q.numerator) / q.denominator)
    c = mp.mpf(LL.numerator) / LL.denominator * mp.mpf(phi.numerator) / phi.denominator

    def h(x):
        return Dm.Pi(x, c * c / x, c)
    lo, hi = max(c, mp.mpf('0.01')), mp.mpf(50)
    if not (h(lo) < 0 < h(hi)):
        raise RuntimeError('float bracket for xi_3 failed')
    for _ in range(110):
        mid = (lo + hi) / 2
        if h(mid) < 0:
            lo = mid
        else:
            hi = mid
    r = (lo + hi) / 2
    mp.dps = 15
    return Q(str(mp.nstr(r, 25)))


class BoxFailure(Exception):
    pass


def certify_phi_box(Bk, D, phi_lo, phi_hi, max_depth):
    c_lo, c_hi = LL * phi_lo, LL * phi_hi
    C = Bk.ival(c_lo, c_hi)
    # Cauchy root bounds in the rescaled coordinates (uniform down to c = 0)
    mags = [max(abs(Bk.lo(x)), abs(Bk.hi(x))) for x in D.quartic(C)]
    xi_max = Q(math.ceil((c_hi + max(mags[1:5]) / LEAD_A) * 1024), 1024)
    eta_max = Q(math.ceil((c_hi + max(mags[0:4]) / LEAD_B) * 1024), 1024)
    # (i) isolation of xi_3 = c v_3 in chart xi
    g = [float_root_xi3(t) for t in (phi_lo, (phi_lo + phi_hi) / 2, phi_hi)]
    centre = g[1]
    w = Q(3, 2) * (max(g) - min(g)) + Q(1, 10 ** 9)
    for _attempt in range(8):
        X_lo, X_hi = centre - w, centre + w
        h_lo = D.Pi(Bk.ival(X_lo, X_lo), C * C / Bk.ival(X_lo, X_lo), C)
        h_hi = D.Pi(Bk.ival(X_hi, X_hi), C * C / Bk.ival(X_hi, X_hi), C)
        X = Bk.ival(X_lo, X_hi)
        eX = C * C / X                                   # encloses eta = c^2/xi
        Px, Pe = D.Pi_partials(X, eX, C)
        dh = Px - (eX / X) * Pe                          # encloses h_c'(xi)
        if Bk.hi(h_lo) < 0 and Bk.lo(h_hi) > 0 and Bk.lo(dh) > 0 and X_lo > c_hi:
            break
        w *= 2
    else:
        raise BoxFailure('isolation of xi_3')
    # G decreases on X up to xi_3(c) and increases after it, so
    # W(xi_3(c)) <= W(centre, c^2/centre; c) for every c in C.
    Pc = Bk.ival(centre, centre)
    W_up = Bk.hi(D.W(Pc, C * C / Pc, C))
    # (ii) preference for block B at xi_3: (1-psi_B)/L_B > (1-psi_A)/L_A, i.e.
    # beta/lambda_B > alpha/lambda_A.  Admissibility 0 < psi < 1 of every positive
    # root is S3 Section 2.3 (in the limit c = 0, psi_A = 1: the A loading vanishes).
    pA, pB = D.psi(X, eX, C)
    pref = (1 - pB) / D.LB - (1 - pA) / D.LA
    if not Bk.lo(pref) > 0:
        raise BoxFailure('preference for block B at xi_3')
    # (iii) exclusion on the rest of both charts
    best_ratio = None
    nboxes = 0
    stats = {'no_root': 0, 'worse': 0, 'worse_above_X': 0}
    pieces = [('xi', c_lo, X_lo), ('xi', X_hi, xi_max), ('eta', c_lo, eta_max)]
    for piece, (chart, lo0, hi0) in enumerate(pieces):
        if hi0 <= lo0:
            continue
        n0 = 16
        w_min = (hi0 - lo0) / 2 ** 10           # profile test only on boxes this small
        stack = [(lo0 + (hi0 - lo0) * i / n0, lo0 + (hi0 - lo0) * (i + 1) / n0, 0)
                 for i in range(n0)]
        while stack:
            lo, hi, depth = stack.pop()
            nboxes += 1
            # other coordinate: c^2/t with t in [lo, hi], c in C, and <= c (chart condition)
            other_lo = c_lo * c_lo / hi
            other_hi = c_hi if lo == 0 else min(c_hi, c_hi * c_hi / lo)
            T = Bk.ival(lo, hi)
            O = Bk.ival(other_lo, other_hi)
            x_, e_ = (T, O) if chart == 'xi' else (O, T)
            pv = D.Pi(x_, e_, C)
            if Bk.lo(pv) > 0 or Bk.hi(pv) < 0:
                stats['no_root'] += 1
                continue
            if hi - lo > w_min and depth < max_depth:
                mid = (lo + hi) / 2
                stack.append((lo, mid, depth + 1))
                stack.append((mid, hi, depth + 1))
                continue
            w_lo = Bk.lo(D.W(x_, e_, C))
            if w_lo > W_up:
                stats['worse'] += 1
                if piece == 1:
                    stats['worse_above_X'] += 1
                ratio = w_lo / W_up
                best_ratio = ratio if best_ratio is None else min(best_ratio, ratio)
                continue
            if depth >= max_depth:
                raise BoxFailure(f'exclusion in chart {chart} near '
                                 f'[{float(lo):.6g}, {float(hi):.6g}]')
            mid = (lo + hi) / 2
            stack.append((lo, mid, depth + 1))
            stack.append((mid, hi, depth + 1))
    return {'X': (X_lo, X_hi), 'ratio': best_ratio, 'boxes': nboxes, 'stats': stats,
            'xi_max': xi_max, 'eta_max': eta_max, 'W_up': W_up}


def log_lower(ratio):
    """Lower bound for log(ratio), ratio a positive Fraction (mpmath.iv; reporting only)."""
    old = iv.prec
    iv.prec = 120
    x = iv.log(iv.mpf(ratio.numerator) / iv.mpf(ratio.denominator))
    lo = Q(*to_rational(x._mpi_[0]))
    iv.prec = old
    return lo


def part_C(Bk, phi_max, width, max_depth, max_phi_depth, verbose):
    print(f'\n== Part C: interval certificate on 0 < phi <= {phi_max} '
          f'(backend {Bk.name}) ==')
    t0 = time.time()
    D = Bk.design()
    n_init = int(phi_max / width)
    todo = [(phi_max * i / n_init, phi_max * (i + 1) / n_init, 0) for i in range(n_init)]
    todo.reverse()
    done = []
    while todo:
        lo, hi, d = todo.pop()
        try:
            res = certify_phi_box(Bk, D, lo, hi, max_depth)
            done.append((lo, hi, res))
        except BoxFailure as e:
            if d >= max_phi_depth:
                fail(f'[{Bk.name}] phi subinterval [{lo}, {hi}] not certified: {e}')
                done.append((lo, hi, None))
                continue
            mid = (lo + hi) / 2
            todo.append((mid, hi, d + 1))
            todo.append((lo, mid, d + 1))
    done.sort(key=lambda z: z[0])
    ok_boxes = [z for z in done if z[2] is not None]
    # coverage check: consecutive exact rational endpoints from 0 to phi_max
    cover = done[0][0] == 0 and done[-1][1] == phi_max and all(
        done[i][1] == done[i + 1][0] for i in range(len(done) - 1))
    if not cover:
        fail(f'[{Bk.name}] phi subintervals do not tile [0, {phi_max}]')
    if verbose:
        print('  phi subinterval   X (encloses xi_3 = c v_3)          v-boxes  '
              'certified gap G(other) - G(v_3) >=')
        for lo, hi, r in done:
            if r is None:
                print(f'  [{lo}, {hi}]  NOT CERTIFIED')
                continue
            gtxt = (f'{float(log_lower(r["ratio"])):.6f}' if r['ratio'] is not None
                    else 'no other stationary point')
            print(f'  [{float(lo):.4f}, {float(hi):.4f}]  [{float(r["X"][0]):.10f}, '
                  f'{float(r["X"][1]):.10f}]  {r["boxes"]:5d}  {gtxt}')
    ratios = [z[2]['ratio'] for z in ok_boxes if z[2]['ratio'] is not None]
    gap = log_lower(min(ratios)) if ratios else None
    widths = sorted({z[1] - z[0] for z in done})
    n_no = sum(z[2]['stats']['no_root'] for z in ok_boxes)
    n_worse = sum(z[2]['stats']['worse'] for z in ok_boxes)
    print(f'  phi subintervals: {len(done)} (widths {", ".join(str(w) for w in widths)}); '
          f'certified: {len(ok_boxes)}; exact tiling of [0, {phi_max}]: {cover}')
    print(f'  v-boxes examined in certified subintervals: '
          f'{sum(z[2]["boxes"] for z in ok_boxes)}; decided as "no stationary point": {n_no}; '
          f'as "profile strictly larger than at v_3": {n_worse}; all others were bisected')
    largest = all(z[2]['stats']['worse_above_X'] == 0 for z in ok_boxes)
    print(f'  [{"ok" if largest else "FAIL"}] no stationary point lies above X in any subinterval, '
          f'so v_3 is the largest positive root of P')
    if not largest:
        fail(f'[{Bk.name}] v_3 not certified to be the largest positive root')
    if done and done[0][2] is not None:
        print(f'  the first subinterval [0, {done[0][1]}] contains c = 0; its certificate holds '
              f'for every phi in (0, {done[0][1]}] simultaneously')
    if gap is not None:
        print(f'  minimum certified profile gap: G(v) - G(v_3) >= {float(gap):.6f} for every '
              f'other positive stationary point v, at every phi in (0, {phi_max}]')
    # final run of subintervals in which no other stationary point is possible
    k = len(ok_boxes)
    while k > 0 and ok_boxes[k - 1][2]['ratio'] is None:
        k -= 1
    if k < len(ok_boxes) and len(ok_boxes) == len(done):
        u0 = ok_boxes[k][0]
        print(f'  v_3 is certified to be the only positive stationary point for phi in '
              f'[{u0}, {phi_max}] = [{float(u0):.6f}, {float(phi_max)}]')
    print('  phi band       subintervals  certified gap (lower bound)')
    band = Q(1, 20)
    k = 0
    while band * k < phi_max:
        b_lo, b_hi = band * k, min(band * (k + 1), phi_max)
        sel = [z for z in ok_boxes if b_lo <= z[0] < b_hi]
        rs = [z[2]['ratio'] for z in sel if z[2]['ratio'] is not None]
        gtxt = f'{float(log_lower(min(rs))):.6f}' if rs else 'no other stationary point'
        print(f'  ({float(b_lo):.2f}, {float(b_hi):.2f}]  {len(sel):12d}  {gtxt:>27}')
        k += 1
    print(f'  time {time.time() - t0:.1f} s')
    return {'n': len(done), 'gap': gap, 'ok': len(ok_boxes) == len(done) and cover}


# ------------------------------------------- Part X (negative control)
def part_X(Bk, max_depth):
    """Run the same procedure on a design where its claim is false somewhere.

    For 2 items at .90 versus 4 at .70, the minimum favouring the longer block
    is global below the Maxwell point phi_M = .3257 and the other minimum is
    global above it (review F1).  A sound procedure must fail on every
    subinterval reaching phi_M or beyond.  A useful one should still certify
    subintervals well below phi_M.
    """
    print(f'\n== Part X: negative control, 2 items at .90 versus 4 at .70 '
          f'(backend {Bk.name}) ==')
    phi_M = Q(325654740638, 10 ** 12)     # review F1 (floating point), used only to sort boxes
    use_design(*CONTROL)
    try:
        D = Bk.design()
        rows = []
        # four subintervals in [.20, .21] (claim true there) and ten in [.325, .35]
        for i in list(range(80, 84)) + list(range(130, 140)):
            lo, hi = Q(i, 400), Q(i + 1, 400)
            try:
                certify_phi_box(Bk, D, lo, hi, max_depth)
                rows.append((lo, hi, True, ''))
            except BoxFailure as e:
                rows.append((lo, hi, False, str(e)))
    finally:
        use_design(*EXAMPLE)
    for lo, hi, ok, why in rows:
        print(f'  [{float(lo):.4f}, {float(hi):.4f}]  '
              f'{"certified" if ok else "not certified (" + why + ")"}')
    sound = all(not ok for lo, hi, ok, _ in rows if hi > phi_M)
    useful = all(ok for lo, hi, ok, _ in rows if hi <= Q(21, 100))
    print(f'  [{"ok" if sound else "FAIL"}] no subinterval reaching phi_M = .3257 or beyond '
          f'is certified')
    print(f'  [{"ok" if useful else "FAIL"}] the subintervals in [.20, .21], well below phi_M, '
          f'are certified')
    if not (sound and useful):
        fail('Part X: negative control did not behave as required')


# --------------------------------------------------------------------- main
def load_verifier():
    path = Path(__file__).resolve().parents[1] / 'proofs' / 'unequal_blocks' / 'verify_unequal_blocks.py'
    if not path.exists():
        return None
    try:
        spec = importlib.util.spec_from_file_location('verify_unequal_blocks', path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod
    except Exception as e:                      # missing numpy/scipy etc.
        print(f'  (could not import verifier: {e!r}; cross-check skipped)')
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--backend', choices=['iv', 'exact', 'both'], default='both')
    ap.add_argument('--phi-max', default='51/100', help='upper end (rational), default 51/100')
    ap.add_argument('--width', default='1/400', help='initial phi subinterval width')
    ap.add_argument('--verbose', action='store_true', help='print every phi subinterval')
    ap.add_argument('--max-depth', type=int, default=40, help='v-box bisection depth')
    ap.add_argument('--max-phi-depth', type=int, default=8, help='phi bisection depth')
    ap.add_argument('--iv-prec', type=int, default=113, help='mpmath.iv precision (bits)')
    args = ap.parse_args()
    t0 = time.time()
    print('Interval certificate: 2 items at .85 (block A) versus 3 items at .75 (block B)')
    print(f'Python {platform.python_version()} | mpmath {mpmath.__version__} '
          f'(backend {mpmath.libmp.BACKEND}) | sympy {sp.__version__} | {platform.platform()}')
    print(f'options: backend={args.backend} phi_max={args.phi_max} width={args.width} '
          f'max_depth={args.max_depth} max_phi_depth={args.max_phi_depth} iv_prec={args.iv_prec}')
    verifier = load_verifier()
    mp.dps = 15
    part_E(verifier)
    part_N()
    phi_max, width = Q(args.phi_max), Q(args.width)
    backends = {'iv': [IVBackend(args.iv_prec)], 'exact': [ExactBackend()],
                'both': [IVBackend(args.iv_prec), ExactBackend()]}[args.backend]
    results = [part_C(Bk, phi_max, width, args.max_depth, args.max_phi_depth, args.verbose)
               for Bk in backends]
    part_X(backends[-1], args.max_depth)
    print('\n== Summary ==')
    for Bk, r in zip(backends, results):
        gtxt = f'{float(r["gap"]):.6f}' if r['gap'] is not None else 'n/a'
        print(f'  backend {Bk.name:5s}: {r["n"]} phi subintervals, all certified: {r["ok"]}, '
              f'min certified gap {gtxt}')
    print(f'  total time {time.time() - t0:.1f} s')
    if FAILURES:
        print(f'RESULT: NOT CERTIFIED ({len(FAILURES)} failure(s))')
        return 1
    print(f'RESULT: CERTIFIED. For every phi in (0, {phi_max}] the global minimiser of the reduced '
          'profile is the largest positive root v_3, it favours block B (alpha/.85 < beta/.75), and '
          'every other positive stationary root has a strictly larger profile value.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
