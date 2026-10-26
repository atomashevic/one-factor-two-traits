#!/usr/bin/env python3
"""Independent verification of one-factor ML fits beyond exchangeability.

Run: python3 verify_unequal_blocks.py --output verification_results_unequal.json
Dependencies: numpy, scipy, sympy, mpmath.
Default: 52 designs, 10 random full-parameter starts each, seed 20261008;
80-digit folds and scalar solves, 60-digit mirrored thresholds; exact SymPy
identities and Sturm certificates. The reference program is never imported or
executed. Supplied JSON numbers are used only in the final comparison stage.

Random optimization supports implementation correctness, not mathematical proof.
Exact symbolic identities are distinguished from numerical residual checks.
All verification routines and baseline numbers are included in this one file.
"""


# Independent component: verification_core.py

import numpy as np

from numpy.polynomial import Polynomial

from scipy.optimize import minimize

def make_population(n, k, la, lb, phi):
    la = np.broadcast_to(np.asarray(la, dtype=float), (n,))
    lb = np.broadcast_to(np.asarray(lb, dtype=float), (k,))
    R = np.empty((n + k, n + k))
    R[:n, :n] = np.multiply.outer(la, la)
    R[n:, n:] = np.multiply.outer(lb, lb)
    R[:n, n:] = phi * np.multiply.outer(la, lb)
    R[n:, :n] = R[:n, n:].T
    np.fill_diagonal(R, 1.0)
    return R

class RankOneLikelihood:
    def __init__(self, R):
        self.R = np.asarray(R, dtype=float)
        self.p = len(R)
        sign, self.logdet = np.linalg.slogdet(R)
        if sign != 1:
            raise ValueError('Population matrix must be positive definite')

    def value_score(self, parameters):
        l, eta = np.split(np.asarray(parameters), 2)
        reciprocal = np.exp(-eta)
        u = reciprocal * l
        gamma = float(l @ u)
        t = 1.0 + gamma
        Ru = self.R @ u
        z = float(u @ Ru)
        value = (np.diag(self.R) @ reciprocal - z / t + eta.sum()
                 + np.log1p(gamma) - self.logdet - self.p)
        coefficient = (t + z) / (t * t)
        loading_score = 2.0 * (coefficient * u - reciprocal * Ru / t)
        uniqueness_score = (1.0 - np.diag(self.R) * reciprocal
                            + 2.0 * u * Ru / t - coefficient * l * u)
        return float(value), np.r_[loading_score, uniqueness_score]

    def dense_value(self, parameters):
        l, eta = np.split(parameters, 2)
        S = np.diag(np.exp(eta)) + np.multiply.outer(l, l)
        return (float(np.trace(np.linalg.solve(S, self.R)))
                + np.linalg.slogdet(S)[1] - self.logdet - self.p)

    def hessian(self, parameters):
        """Analytic directional derivative of the full likelihood gradient."""
        l, eta = np.split(parameters, 2)
        psi = np.exp(eta)
        reciprocal = 1.0 / psi
        u = reciprocal * l
        K = np.diag(reciprocal) - np.outer(u, u) / (1.0 + l @ u)
        E = K @ self.R @ K
        D = K - E
        H = np.empty((2*self.p, 2*self.p))
        for col in range(2*self.p):
            dl = np.zeros(self.p)
            de = np.zeros(self.p)
            if col < self.p:
                dl[col] = 1.0
            else:
                de[col-self.p] = 1.0
            dS = np.outer(l, dl) + np.outer(dl, l) + np.diag(psi * de)
            dD = -K @ dS @ K + K @ dS @ E + E @ dS @ K
            H[:, col] = np.r_[2.0*(dD @ l + D @ dl),
                              psi*np.diag(dD) + psi*de*np.diag(D)]
        return (H+H.T)/2.0

    def optimize(self, initial):
        result = minimize(self.value_score, initial, method='L-BFGS-B', jac=True,
                          bounds=[(None, None)]*self.p + [(-24.0, 4.0)]*self.p,
                          options={'maxiter':4000, 'maxls':50,
                                   'ftol':2e-15, 'gtol':1e-10})
        # A few analytic Newton corrections remove premature function-tolerance
        # stops.  Do not cross an indefinite Hessian or accept an uphill step.
        z = result.x.copy()
        for _ in range(8):
            f, g = self.value_score(z)
            if np.linalg.norm(g, np.inf) < 2e-11:
                break
            H = self.hessian(z)
            eig = np.linalg.eigvalsh(H)
            if eig[0] < 1e-8:
                break
            step = np.linalg.solve(H, g)
            if np.linalg.norm(step) > 0.25:
                break
            trial = z-step
            if np.any(trial[self.p:] < -24) or np.any(trial[self.p:] > 4):
                break
            f2, g2 = self.value_score(trial)
            if f2 <= f+2e-13 and np.linalg.norm(g2) < np.linalg.norm(g):
                z = trial
            else:
                break
        f, g = self.value_score(z)
        return {'theta':z, 'F':f, 'score':np.linalg.norm(g,np.inf),
                'iterations':int(result.nit), 'optimizer_success':bool(result.success)}

def stationary_quartic(n, k, la, lb, phi):
    """Ascending coefficients, no squaring or resultant roots are introduced."""
    L, M = la*la, lb*lb
    a, b = 1+(n-1)*L, 1+(k-1)*M
    A, B, c, q = n-a, k-b, la*lb*phi, n+k-1
    C, D = k*c, n*c
    E = n*b+k*a-a*b-C*D
    v = Polynomial([0.0,1.0])
    U = Polynomial([a,C])
    Vbar = Polynomial([D,b])
    Nbar = Polynomial([D*A,E,C*B])
    Jbar = q*v*U-Nbar
    Kbar = q*Vbar-Nbar
    return v*Vbar*Kbar-U*Jbar

def stationary_candidates(n, k, la, lb, phi):
    """Positive stationary points reconstructed independently from a quartic."""
    p = n+k
    a, b = 1+(n-1)*la*la, 1+(k-1)*lb*lb
    A, B, c, q = n-a, k-b, la*lb*phi, p-1
    C, D = k*c, n*c
    E = n*b+k*a-a*b-C*D
    objective = RankOneLikelihood(make_population(n,k,la,lb,phi))
    out=[]
    for root in stationary_quartic(n,k,la,lb,phi).roots():
        if abs(root.imag) > 1e-7 or root.real <= 0:
            continue
        v = float(root.real)
        U, V = a+C*v, b+D/v
        N = E+C*B*v+D*A/v
        x,y = N/(q*V),N/(q*U)
        if not (0 < x < 1 and 0 < y < 1):
            continue
        alpha,beta = np.sqrt(1-x),np.sqrt(1-y)
        theta=np.r_[np.repeat(alpha,n),np.repeat(beta,k),
                    np.repeat(np.log(x),n),np.repeat(np.log(y),k)]
        f,g=objective.value_score(theta)
        evals=np.linalg.eigvalsh(objective.hessian(theta))
        fa,fb=alpha/la,beta/lb
        # These positive stationary solutions have positive composite weights.
        share_A=(fa-phi*fb)/((1-phi)*(fa+fb))
        out.append({'v':v,'x':x,'y':y,'alpha':alpha,'beta':beta,
                    'theta':theta,'F':f,'share_A':share_A,
                    'score':float(np.linalg.norm(g,np.inf)),
                    'hessian_eigenvalues':evals,
                    'local_minimum':bool(evals[0]>-1e-8)})
    return sorted(out,key=lambda r:r['v'])

def profile_point(n,k,la,lb,phi,d):
    """Direct two by two eigendecomposition, independent of quartic algebra."""
    aA,aB=1-la*la,1-lb*lb
    bA,bB=1+(n-1)*la*la,1+(k-1)*lb*lb
    h=la*lb*phi*np.sqrt(n*k)
    diag=np.array([bA*np.exp(-d),bB*np.exp(d)])
    eigenvalues,eigenvectors=np.linalg.eigh(np.array([[diag[0],h],[h,diag[1]]]))
    T=(n-1)*aA*np.exp(-d)+(k-1)*aB*np.exp(d)+eigenvalues[0]
    t=T/(n+k-1)
    x,y=t*np.exp(d),t*np.exp(-d)
    theta1=eigenvalues[1]/t
    e=eigenvectors[:,1]
    if e.sum()<0:
        e=-e
    alpha=np.sqrt(x*(theta1-1)/n)*e[0]
    beta=np.sqrt(y*(theta1-1)/k)*e[1]
    G=((n+k-1)*np.log(t)+(n-k)*d-np.log(eigenvalues[0])
       -(n-1)*np.log(aA)-(k-1)*np.log(aB))
    return {'F':float(G),'alpha':float(alpha),'beta':float(beta),
            'x':float(x),'y':float(y),'s':float(np.log(t))}


# Independent component: scalar_math.py

from math import sqrt, log

def quartic_coefficients(n, k, loading_a, loading_b, phi):
    """Coefficients in descending order; accepts ordinary or mpmath scalars."""
    la, lb = loading_a**2, loading_b**2
    a, b = 1 + (n - 1)*la, 1 + (k - 1)*lb
    aa, bb = n - a, k - b
    c = loading_a*loading_b*phi
    cc, dd = k*c, n*c
    q = n + k - 1
    ee = n*b + k*a - a*b - cc*dd
    jj, kk = q*a - ee, q*b - ee
    ff, gg = q - bb, q - aa
    return [
        -b*cc*bb,
        b*kk - dd*cc*bb - cc*cc*ff,
        b*dd*gg + dd*kk - a*cc*ff - cc*jj,
        dd*dd*gg - a*jj + cc*dd*aa,
        a*dd*aa,
    ]

def reconstruct(n, k, loading_a, loading_b, phi, v):
    """Return (x,y,tau) at the profiled scale for any positive v."""
    ba = 1 + (n - 1)*loading_a**2
    bb = 1 + (k - 1)*loading_b**2
    c = loading_a*loading_b*phi
    u, w = ba + k*c*v, bb + n*c/v
    nn = n*w + k*u - u*w
    q = n + k - 1
    return nn/(q*w), nn/(q*u), q*u*w/nn

def profile_value(n, k, loading_a, loading_b, phi, v, logarithm=log):
    la, lb = loading_a**2, loading_b**2
    ba, bb = 1 + (n - 1)*la, 1 + (k - 1)*lb
    c = loading_a*loading_b*phi
    u, w = ba + k*c*v, bb + n*c/v
    nn = n*w + k*u - u*w
    q = n + k - 1
    determinant = ba*bb - n*k*c*c
    return (q*logarithm(nn) - (k - 1)*logarithm(u)
            - (n - 1)*logarithm(w) - q*logarithm(q)
            - logarithm(determinant)
            - (n - 1)*logarithm(1 - la)
            - (k - 1)*logarithm(1 - lb))

def share_a(n, k, loading_a, loading_b, phi, v, square_root=sqrt):
    x, y, _ = reconstruct(n, k, loading_a, loading_b, phi, v)
    f, g = square_root(1-x)/loading_a, square_root(1-y)/loading_b
    return (f - phi*g)/((1-phi)*(f+g))

def share_limit(n, k, loading_a, loading_b):
    ra, rb = loading_a**2/(1-loading_a**2), loading_b**2/(1-loading_b**2)
    ja, jb = 2*(n-1)*ra, 2*(k-1)*rb
    ka, kb = n*ra, k*rb
    denominator = 1 + kb/ja + ka/jb
    u, v = (kb/ja)/denominator, (ka/jb)/denominator
    return (1 + v - u)/2

def share_crossing(n, k, loading_a, loading_b, square_root=sqrt):
    """Algebraic candidate, admissible only if 0<t<1 and 0<phi<1."""
    one = loading_a*0 + 1
    ratio = square_root(one*k*(k-1)/(n*(n-1)))
    if ratio == 1:
        raise ValueError('Equal sizes have no isolated normalized-share crossing')
    la, lb = loading_a**2, loading_b**2
    t = (ratio*lb - la)/((ratio-1)*la*lb)
    r = square_root(one*(n-1)*(k-1)/(n*k))
    return t, (1+r)*t-r

def weaker_endpoint_stable(m, loading_a, loading_b):
    """Equal sizes and loading_a<loading_b. Strict iff a positive fold exists."""
    la, lb = loading_a**2, loading_b**2
    return la > (m-1)*lb/(m+(m-1)*lb)

def scalar_symbolic_checks():
    """Exact identities checked with SymPy, using independent symbolic formulas."""
    import sympy as s
    v,n,k,a,b,c = s.symbols('v n k a b c', nonzero=True)
    q=n+k-1
    A,B,C,D=n-a,k-b,k*c,n*c
    U,V=a+C*v,b+D/v
    N=s.expand(n*V+k*U-U*V)
    E=n*b+k*a-a*b-C*D
    J,K=q*a-E,q*b-E
    f,g=q-B,q-A
    coefficients=[-b*C*B,b*K-D*C*B-C*C*f,
                  b*D*g+D*K-a*C*f-C*J,D*D*g-a*J+C*D*A,a*D*A]
    polynomial=sum(z*v**(4-i) for i,z in enumerate(coefficients))
    explicit=s.cancel(v*(v*v*V*(q*V-N)-U*(q*U-N)))
    assert s.cancel(explicit-polynomial)==0
    derivative=q*s.diff(N,v)/N-(k-1)*s.diff(U,v)/U-(n-1)*s.diff(V,v)/V
    multiplier=-c*n*k/(v*U*(v*V)*(v*N))
    assert s.cancel(derivative-multiplier*polynomial)==0
    x,alpha,tau,aa,bbeta,y=s.symbols('x alpha tau aa bbeta y', nonzero=True)
    H=s.Matrix([[2/x-2*aa/(tau*x*x),2*aa*alpha/(tau*x**3)],
                [2*aa*alpha/(tau*x**3),1/x**2-2*aa*alpha**2/(tau*x**4)]])
    det_expected=2/x**3*(1-aa*(x+2*alpha**2)/(tau*x*x))
    assert s.cancel(H.det()-det_expected)==0
    certificate=unequal_example_certificate()
    return {'quartic_identity':True,'derivative_identity':True,'contrast_determinant':True,
            'unequal_example_discriminant_identity':certificate['discriminant_identity'],
            'unequal_example_no_discriminant_roots_half_to_one':certificate['discriminant_roots_half_to_one']==0,
            'unequal_example_one_positive_root_at_half':certificate['positive_quartic_roots_at_half']==1,
            'unequal_example_share_derivative_identity':certificate['share_crossing_transverse']}

def unequal_example_certificate():
    """Exact algebraic/Sturm certificate for 2*.85 versus 3*.75.

    No floating point numbers enter any assertion.
    """
    import sympy as s
    v,phi=s.symbols('v phi')
    n,k=2,3
    loading_a,loading_b=s.Rational(17,20),s.Rational(3,4)
    coefs=quartic_coefficients(n,k,loading_a,loading_b,phi)
    polynomial=sum(z*v**(4-i) for i,z in enumerate(coefs))
    H=(s.Integer(80519601508683987527503641)*phi**10
       +s.Integer(2682578090235904365337528572)*phi**8
       -s.Integer(7230857150333348166462593850)*phi**6
       +s.Integer(5608822138312392897886943748)*phi**4
       -s.Integer(1711849644321614301603840471)*phi**2
       +s.Integer(177529911547540532703125000))
    expected=(s.Rational(4913,70368744177664000000000000000000)
              *(459*phi**2-689)**2*H)
    assert s.expand(s.discriminant(polynomial,v)-expected)==0
    assert s.Poly(H,phi).count_roots(s.Rational(1,2),1)==0
    assert s.Poly(459*phi**2-689,phi).count_roots(s.Rational(1,2),1)==0
    assert s.Poly(polynomial.subs(phi,s.Rational(1,2)),v).count_roots(0,s.oo)==1
    t0=s.Rational(3088,2601)-s.Rational(512,2601)*s.sqrt(3)
    phi0=s.Rational(2576,2601)-s.Rational(1049,7803)*s.sqrt(3)
    t2,p2=share_crossing(n,k,loading_a,loading_b,s.sqrt)
    assert s.simplify(t0-t2)==0 and s.simplify(phi0-p2)==0
    assert phi0>s.Rational(1,2) and phi0<1
    a,b=loading_a**2,loading_b**2
    x,y=1-a*t0,1-b*t0
    r=phi0-t0
    A11=-2*(n-1)*a*(1-a)*t0/x**2+k*b*(-r-t0)/y
    A12=k*b*((r-t0)/y+2*b*t0*r/y**2)
    A21=n*a*((r-t0)/x+2*a*t0*r/x**2)
    A22=-2*(k-1)*b*(1-b)*t0/y**2+n*a*(-r-t0)/x
    J=s.Matrix([[A11,A12],[A21,A22]]).applyfunc(s.simplify)
    forcing=s.Matrix([k*b/y,n*a/x])
    uv=J.inv()*(-forcing)
    derivative=s.simplify((1+phi0)/(4*(1-phi0))*(uv[0]-uv[1]))
    expected_derivative=(315514305-180683667*s.sqrt(3))/5739392
    assert s.simplify(derivative-expected_derivative)==0
    assert s.Integer(315514305)**2>3*s.Integer(180683667)**2
    return {'discriminant_identity':True,'discriminant_roots_half_to_one':0,
            'positive_quartic_roots_at_half':1,'crossing_expression':str(phi0),
            'share_derivative_expression':str(expected_derivative),
            'share_crossing_transverse':True}


# Independent component: cusp_math.py

import json

from pathlib import Path

import mpmath as mp

import numpy as np

import sympy as sp

DPS = 80

ROOT_TOL = mp.mpf("1e-65")

VERIFY_TOL = mp.mpf("1e-55")

def equal_coefficients(m, lambda_a, lambda_b, phi):
    """Coefficients, in descending powers, of the exact scalar quartic."""
    m, la, lb, phi = map(mp.mpf, (m, lambda_a, lambda_b, phi))
    q = 2 * m - 1
    ba = 1 + (m - 1) * la * la
    bb = 1 + (m - 1) * lb * lb
    aa, ab = m - ba, m - bb
    h = m * la * lb * phi
    e = m * (ba + bb) - ba * bb - h * h
    return [
        -bb * h * ab,
        bb * (q * bb - e) - q * h * h,
        h * (3 * m - 2) * (bb - ba),
        q * h * h - ba * (q * ba - e),
        ba * h * aa,
    ]

def poly_value(coefficients, v, order=0):
    co = list(coefficients)
    for _ in range(order):
        n = len(co) - 1
        co = [(n - j) * x for j, x in enumerate(co[:-1])]
    return mp.polyval(co, v)

def profile_g(d, m, lambda_a, lambda_b, phi):
    """Exact profile in x=exp(s+d), y=exp(s-d), after minimizing s."""
    m, la, lb, phi = map(mp.mpf, (m, lambda_a, lambda_b, phi))
    n, q = m - 1, 2 * m - 1
    aa, ab = 1 - la * la, 1 - lb * lb
    ba, bb = 1 + n * la * la, 1 + n * lb * lb
    h = m * la * lb * phi
    determinant = ba * bb - h * h
    trace = ba * mp.exp(-d) + bb * mp.exp(d)
    discriminant = (ba * mp.exp(-d) - bb * mp.exp(d)) ** 2 + 4 * h * h
    r = 2 * determinant / (trace + mp.sqrt(discriminant))
    qsum = n * aa * mp.exp(-d) + n * ab * mp.exp(d) + r
    return q * mp.log(qsum / q) - n * mp.log(aa) - n * mp.log(ab) - mp.log(r)

def candidate_from_root(m, lambda_a, lambda_b, phi, v):
    """Recover all fitted quantities from a positive quartic root."""
    m, la, lb, phi, v = map(mp.mpf, (m, lambda_a, lambda_b, phi, v))
    q = 2 * m - 1
    ba, bb = 1 + (m - 1) * la * la, 1 + (m - 1) * lb * lb
    aa, ab = m - ba, m - bb
    h = m * la * lb * phi
    e = m * (ba + bb) - ba * bb - h * h
    u, w = ba + h * v, bb + h / v
    numerator = e + h * ab * v + h * aa / v
    x, y = numerator / (q * w), numerator / (q * u)
    if not (0 < x < 1 and 0 < y < 1):
        raise ValueError("Nonphysical quartic root")
    alpha, beta = mp.sqrt(1 - x), mp.sqrt(1 - y)
    d = mp.log(x / y) / 2
    signed_tilt = (1 + phi) / (2 * (1 - phi)) * (beta / lb - alpha / la) / (beta / lb + alpha / la)
    gamma = m * alpha * alpha / x + m * beta * beta / y
    score_a = ba * alpha / x + m * la * lb * phi * beta / y - (1 + gamma) * alpha
    score_b = bb * beta / y + m * la * lb * phi * alpha / x - (1 + gamma) * beta
    curvature = mp.diff(lambda z: profile_g(z, m, la, lb, phi), d, 2)
    return {
        "v": v,
        "x": x,
        "y": y,
        "alpha": alpha,
        "beta": beta,
        "d": d,
        "F": profile_g(d, m, la, lb, phi),
        "signed_B_tilt": signed_tilt,
        "share": mp.mpf(".5") + abs(signed_tilt),
        "curvature_G_dd": curvature,
        "score_residual": max(abs(score_a), abs(score_b)),
        "profile_score_residual": abs(mp.diff(lambda z: profile_g(z, m, la, lb, phi), d)),
        "quartic_residual": abs(poly_value(equal_coefficients(m, la, lb, phi), v)),
    }

def positive_stationary_candidates(m, lambda_a, lambda_b, phi):
    """Refine every positive real quartic root to current mpmath precision.

    Double-precision polynomial roots are used only as starting guesses.
    Their final residuals and the independent profile scores are checked.
    A double root itself should instead be computed with solve_fold.
    """
    coeff = equal_coefficients(m, lambda_a, lambda_b, phi)
    guesses = np.roots([float(c) for c in coeff])
    roots = []
    for z in guesses:
        if z.real <= 0 or abs(z.imag) > 1e-7 * max(1.0, abs(z.real)):
            continue
        root = mp.findroot(
            lambda v: poly_value(coeff, v),
            mp.mpf(str(z.real)),
            df=lambda v: poly_value(coeff, v, 1),
            tol=ROOT_TOL,
            maxsteps=100,
        )
        if root > 0 and not any(abs(root - old) < mp.mpf("1e-45") for old in roots):
            roots.append(root)
    rows = [candidate_from_root(m, lambda_a, lambda_b, phi, v) for v in sorted(roots)]
    for row in rows:
        if max(row["score_residual"], row["profile_score_residual"]) > VERIFY_TOL:
            raise ArithmeticError("Quartic and exact-profile scores disagree")
    return rows

def cusp_coefficients(lam, m):
    lam, m = map(mp.mpf, (lam, m))
    l = lam * lam
    a, n, q = 1 - l, m - 1, 2 * m - 1
    b = 1 + n * l
    phi_star = b / (m * (2 - l))
    k1 = 2 * m * n * lam * (2 - l) / (a * b)
    k2 = n * m * m * (2 - l) ** 3 / (a * q * b)
    k4 = 2 * n * m * (2 - l) ** 2 / (q * l * l)
    d_coeff = (k1 / k4) ** (mp.mpf(1) / 3)
    fold_coeff = 3 * k4 / k2 * (k1 / (2 * k4)) ** (mp.mpf(2) / 3)
    tilt_factor = (2 * m + 1 - l) / (2 * q * l)
    return {
        "phi_star": phi_star,
        "k1": k1,
        "k2": k2,
        "k4": k4,
        "d_coefficient": d_coeff,
        "fold_shift_coefficient": fold_coeff,
        "tilt_conversion_from_abs_d": tilt_factor,
        "tilt_coefficient": d_coeff * tilt_factor,
        "fold_over_d_squared_limit": 3 * k4 / (mp.power(2, mp.mpf(2) / 3) * k2),
        "fold_over_tilt_squared_limit": 3 * k4 / (mp.power(2, mp.mpf(2) / 3) * k2 * tilt_factor ** 2),
    }

def solve_fold(m, lambda_a, lambda_b, start=None):
    """Solve the exact P(v,phi)=P_v(v,phi)=0 condition with mpmath."""
    m, la, lb = map(mp.mpf, (m, lambda_a, lambda_b))
    if not la < lb:
        raise ValueError("Use lambda_a < lambda_b")
    if m * la * la / (1 - la * la) <= (m - 1) * lb * lb:
        raise ValueError("The weaker-block local minimum has no positive fold")
    lam, eps = (la + lb) / 2, (lb - la) / 2
    constants = cusp_coefficients(lam, m)
    if start is None:
        d = -(constants["k1"] * eps / (2 * constants["k4"])) ** (mp.mpf(1) / 3)
        v0 = mp.exp(2 * d / (lam * lam))
        phi0 = max(mp.mpf(".001"), constants["phi_star"] - 3 * constants["k4"] * d * d / constants["k2"])
    else:
        v0, phi0 = map(mp.mpf, start)
    def equations(v, phi):
        co = equal_coefficients(m, la, lb, phi)
        return poly_value(co, v), poly_value(co, v, 1)
    v, phi = mp.findroot(equations, (v0, phi0), tol=ROOT_TOL, maxsteps=150)
    if not (0 < v < 1 and 0 < phi < 1):
        raise ArithmeticError("Did not obtain the physical weaker-block fold")
    row = candidate_from_root(m, la, lb, phi, v)
    co = equal_coefficients(m, la, lb, phi)
    row.update({
        "epsilon": eps,
        "lambda_A": la,
        "lambda_B": lb,
        "phi_fold": phi,
        "P_v_residual": abs(poly_value(co, v, 1)),
        "P_vv": poly_value(co, v, 2),
        "P_phi": mp.diff(lambda z: poly_value(equal_coefficients(m, la, lb, z), v), phi),
        "G_ddd": mp.diff(lambda z: profile_g(z, m, la, lb, phi), row["d"], 3),
        "G_d_phi": mp.diff(lambda z: mp.diff(lambda d: profile_g(d, m, la, lb, z), row["d"]), phi),
    })
    if max(row["quartic_residual"], row["P_v_residual"], abs(row["curvature_G_dd"]), row["profile_score_residual"]) > VERIFY_TOL:
        raise ArithmeticError("Fold failed the independent exact-profile verification")
    return row

def critical_epsilon(lam, m):
    """Return the extinction boundary, or None if it is outside the domain."""
    lam, m = map(mp.mpf, (lam, m))
    upper = min(lam, 1 - lam)
    def f(e):
        la, lb = lam - e, lam + e
        return m * la * la - (m - 1) * lb * lb * (1 - la * la)
    if f(upper) > 0:
        return None
    lo, hi = mp.mpf(0), upper
    for _ in range(300):
        mid = (lo + hi) / 2
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

def cusp_symbolic_checks():
    """Exact identities checked, without numeric substitution."""
    l, m, h, lam = sp.symbols("L m h lambda", positive=True)
    n, q, a = m - 1, 2 * m - 1, 1 - l
    b, c = 1 + n * l, 2 * n * a + 1 + n * l
    hs = l * b / (2 - l)
    r0, z0 = b - h, c - h
    r2, z2 = b - b ** 2 / h, c - b ** 2 / h
    r4, z4 = b - 4 * b ** 2 / h + 3 * b ** 4 / h ** 3, c - 4 * b ** 2 / h + 3 * b ** 4 / h ** 3
    g2 = q * z2 / z0 - r2 / r0
    g4 = q * (z4 / z0 - 3 * (z2 / z0) ** 2) - (r4 / r0 - 3 * (r2 / r0) ** 2)
    gde = -4 * m * n * lam * (n * l + h) / (h * (c - h))
    k1 = 2 * m * n * lam * (2 - l) / (a * b)
    k2 = n * m ** 2 * (2 - l) ** 3 / (a * q * b)
    k4 = 2 * n * m * (2 - l) ** 2 / (q * l ** 2)
    phi_star = b / (m * (2 - l))
    expressions = {
        "C1_mixed_derivative_G_d_epsilon_equals_minus_k1": gde.subs(h, hs) + k1,
        "C2_G_dd_vanishes_at_threshold": g2.subs(h, hs),
        "C3_G_dd_phi_equals_k2": (m * l * sp.diff(g2, h)).subs(h, hs) - k2,
        "C4_G_dddd_equals_six_k4": g4.subs(h, hs) - 6 * k4,
        "C5_tilt_conversion_identity": a * (1 + phi_star) / (2 * l * (1 - phi_star)) - (2 * m + 1 - l) / (2 * q * l),
    }
    ba, bb, v = sp.symbols("b_A b_B v")
    aa, ab = m - ba, m - bb
    e = m * (ba + bb) - ba * bb - h ** 2
    u, w = ba + h * v, bb + h / v
    num = e + h * ab * v + h * aa / v
    poly_def = v * (v ** 2 * w * (q * w - num) - u * (q * u - num))
    coefficients = [-bb * h * ab, bb * (q * bb - e) - q * h ** 2,
                    h * (3 * m - 2) * (bb - ba), q * h ** 2 - ba * (q * ba - e), ba * h * aa]
    poly_explicit = sum(z * v ** (4 - j) for j, z in enumerate(coefficients))
    expressions["C6_equal_size_quartic_expansion"] = poly_def - poly_explicit
    expressions["C7_exchangeable_cusp_triple_root"] = poly_explicit.subs({ba: b, bb: b, h: hs}) + n * a * l * b ** 2 / (2 - l) * (v - 1) ** 3 * (v + 1)
    phi = sp.Symbol("phi")
    la, lb = sp.Rational(17, 25), sp.Rational(18, 25)
    p_mild = poly_explicit.subs({m: 4, ba: 1 + 3 * la ** 2, bb: 1 + 3 * lb ** 2, h: 4 * la * lb * phi}, simultaneous=True)
    z = phi ** 2
    mild_q = (513630418748951851474432698956382732288 * z ** 4
              - 1076409817540179744160183194931606587648 * z ** 3
              + 464339406756728873014939509768353273328 * z ** 2
              - 84180363108670733656948180878978886248 * z
              + 3540366223188038489387282515648442363)
    expressions["C8_mild_design_discriminant_factorization"] = sp.discriminant(p_mild, v) - sp.Rational(322486272, 211758236813575084767080625169910490512847900390625) * (374544 * phi ** 2 - 595681) ** 2 * mild_q
    return {key: bool(sp.cancel(sp.factor(value)) == 0) for key, value in expressions.items()}

def encode_mp(obj):
    if isinstance(obj, mp.mpf):
        return mp.nstr(obj, 65)
    if isinstance(obj, dict):
        return {key: encode_mp(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [encode_mp(value) for value in obj]
    return obj

def run_cusp_checks(write_path=None):
    mp.mp.dps = DPS
    lam, m = mp.mpf(".7"), 4
    checks = cusp_symbolic_checks()
    constants = cusp_coefficients(lam, m)
    folds = []
    start = None
    for es in (".001", ".0025", ".005", ".01", ".02", ".04", ".06", ".08", ".1", ".11", ".12"):
        eps = mp.mpf(es)
        row = solve_fold(m, lam - eps, lam + eps, start=start)
        start = row["v"], row["phi_fold"]
        folds.append(row)
        checks["fold_double_root_eps_" + es] = max(row["quartic_residual"], row["P_v_residual"], abs(row["curvature_G_dd"])) < VERIFY_TOL
    tilt_rows = []
    for es in (".0005", ".001", ".002", ".004", ".008", ".016"):
        eps = mp.mpf(es)
        rows = positive_stationary_candidates(m, lam - eps, lam + eps, constants["phi_star"])
        if len(rows) != 1:
            raise ArithmeticError("Expected one real positive stationary root at phi_star")
        row = rows[0]
        row["epsilon"] = eps
        row["tilt"] = row["share"] - mp.mpf(".5")
        row["tilt_leading_approximation"] = constants["tilt_coefficient"] * eps ** (mp.mpf(1) / 3)
        tilt_rows.append(row)
        checks["cusp_tilt_score_eps_" + es] = row["score_residual"] < VERIFY_TOL and row["curvature_G_dd"] > 0
    scaling = []
    for es in ("1e-3", "1e-6", "1e-9", "1e-12"):
        eps = mp.mpf(es)
        row = positive_stationary_candidates(m, lam - eps, lam + eps, constants["phi_star"])[0]
        fold = solve_fold(m, lam - eps, lam + eps)
        tilt = row["share"] - mp.mpf(".5")
        scaling.append({
            "epsilon": eps,
            "d": row["d"],
            "tilt": tilt,
            "phi_fold": fold["phi_fold"],
            "d_over_leading_term": row["d"] / (constants["d_coefficient"] * eps ** (mp.mpf(1) / 3)),
            "tilt_over_leading_term": tilt / (constants["tilt_coefficient"] * eps ** (mp.mpf(1) / 3)),
            "fold_shift_over_leading_term": (constants["phi_star"] - fold["phi_fold"]) / (constants["fold_shift_coefficient"] * eps ** (mp.mpf(2) / 3)),
            "fold_shift_over_d_squared": (constants["phi_star"] - fold["phi_fold"]) / row["d"] ** 2,
            "fold_shift_over_tilt_squared": (constants["phi_star"] - fold["phi_fold"]) / tilt ** 2,
        })
    checks["cusp_scaling_ratio_at_1e_minus_12"] = all(abs(scaling[-1][key] - 1) < mp.mpf("1e-6") for key in (
        "d_over_leading_term", "tilt_over_leading_term", "fold_shift_over_leading_term"))
    critical = critical_epsilon(lam, m)
    extinction_residual = abs(m * (lam - critical) ** 2 - (m - 1) * (lam + critical) ** 2 * (1 - (lam - critical) ** 2))
    checks["weaker_minimum_extinction_equation"] = extinction_residual < VERIFY_TOL
    sample_designs = []
    for la, lb in ((".68", ".72"), (".6", ".8")):
        for phi in (mp.mpf(".05"), constants["phi_star"], mp.mpf(".9")):
            rows = positive_stationary_candidates(m, la, lb, phi)
            sample_designs.append({"lambda_A": la, "lambda_B": lb, "phi": phi, "candidates": rows})
    results = {
        "provenance": "Independent derivation from R1 and R2; reference implementation not reused",
        "precision_decimal_digits": DPS,
        "mpmath_version": mp.__version__,
        "sympy_version": sp.__version__,
        "root_solver_tolerance": ROOT_TOL,
        "independent_profile_score_tolerance": VERIFY_TOL,
        "constants_lambda_0_7_m_4": constants,
        "symbolic_identity_checks": {key: value for key, value in checks.items() if key.startswith("C")},
        "folds_lambda_0_7_plus_minus_epsilon_m_4": folds,
        "tilts_at_symmetric_threshold": tilt_rows,
        "asymptotic_scaling_checks": scaling,
        "weak_minimum_extinction": {
            "epsilon_critical": critical,
            "lambda_A_critical": lam - critical,
            "lambda_B_critical": lam + critical,
            "equation_residual": extinction_residual,
        },
        "sample_stationary_designs": sample_designs,
        "checks": checks,
        "all_pass": all(checks.values()),
    }
    if write_path:
        Path(write_path).write_text(json.dumps(encode_mp(results), indent=2) + "\n")
    return results


# Independent component: mirrored_numeric.py

import mpmath as mp

def mp_population(la,lb,phi):
    p=len(la)+len(lb)
    n=len(la)
    out=mp.eye(p)
    full=list(la)+list(lb)
    for i in range(p):
        for j in range(i):
            out[i,j]=out[j,i]=full[i]*full[j]*(1 if (i<n)==(j<n) else phi)
    return out

def mp_full_hessian(loadings,psi,R):
    p=len(loadings)
    l=mp.matrix(loadings)
    u=mp.matrix([l[i]/psi[i] for i in range(p)])
    K=mp.diag([1/t for t in psi])-u*u.T/(1+(l.T*u)[0])
    E=K*R*K
    D=K-E
    H=mp.zeros(2*p)
    for col in range(2*p):
        dl=mp.zeros(p,1)
        de=mp.zeros(p,1)
        if col<p: dl[col]=1
        else: de[col-p]=1
        dS=l*dl.T+dl*l.T+mp.diag([psi[i]*de[i] for i in range(p)])
        dD=-K*dS*K+K*dS*E+E*dS*K
        gl=2*(dD*l+D*dl)
        for i in range(p):
            H[i,col]=gl[i]
            H[p+i,col]=psi[i]*dD[i,i]+psi[i]*de[i]*D[i,i]
    return (H+H.T)/2

def mirrored_scores(loadings,phi,lam):
    den=[1-t*t for t in loadings]
    tau=1+2*sum(t*t/d for t,d in zip(loadings,den))
    S=sum(a*t/d for a,t,d in zip(lam,loadings,den))
    return [tau*t-(1-a*a)*t/d-(1+phi)*a*S
            for a,t,d in zip(lam,loadings,den)]

def antisymmetric_hessian(loadings,phi,lam):
    n=len(lam)
    p=2*n
    psi=[1-t*t for t in loadings]*2
    H=mp_full_hessian(list(loadings)*2,psi,mp_population(lam,lam,phi))
    basis=mp.zeros(2*p,2*n)
    for j in range(n):
        basis[j,j]=basis[p+j,n+j]=1/mp.sqrt(2)
        basis[n+j,j]=basis[p+n+j,n+j]=-1/mp.sqrt(2)
    return basis.T*H*basis

def mirrored_threshold(lam,guess,precision=60):
    with mp.workdps(precision):
        lam=[mp.mpf(str(a)) for a in lam]
        def equations(*unknown):
            l=list(unknown[:-1]);phi=unknown[-1]
            HA=antisymmetric_hessian(l,phi,lam)
            return tuple(mirrored_scores(l,phi,lam))+ (mp.det(HA),)
        start=[mp.mpf('.8')*a for a in lam]+[mp.mpf(str(guess))]
        solution=mp.findroot(equations,tuple(start),tol=mp.mpf('1e-48'),maxsteps=50)
        l=[solution[i] for i in range(len(lam))];phi=solution[len(lam)]
        residual=max(abs(z) for z in equations(*solution))
        eig=mp.eigsy(antisymmetric_hessian(l,phi,lam),eigvals_only=True)
        if not (0<phi<1 and abs(eig[0])<mp.mpf('1e-44') and eig[1]>0):
            raise ArithmeticError('The solved point is not a simple stability threshold')
        cbar=sum(a*a/(1-a*a) for a in lam)/len(lam)
        L=cbar/(1+cbar)
        mean_prediction=(1+(len(lam)-1)*L)/(len(lam)*(2-L))
        return {'loadings':[mp.nstr(t,55) for t in lam],
                'threshold':mp.nstr(phi,55),
                'symmetric_fitted_loadings':[mp.nstr(t,55) for t in l],
                'antisymmetric_hessian_eigenvalues':[mp.nstr(e,55) for e in eig],
                'equation_residual':mp.nstr(residual,5),
                'mean_c_prediction':mp.nstr(mean_prediction,55),
                'precision_dps':precision,'tolerance':'1e-48'}


# Independent component: run_core_checks.py

import json,time

import numpy as np

def plain(z):
    if isinstance(z,np.ndarray): return z.tolist()
    if isinstance(z,np.generic): return z.item()
    if isinstance(z,dict): return {k:plain(v) for k,v in z.items()}
    if isinstance(z,(list,tuple)): return [plain(v) for v in z]
    return z

def run_random_checks(seed=20261008,random_count=40,starts=10):
    rng=np.random.default_rng(seed)
    designs=[(2,2,.7,.7,.001),(2,2,.7,.7,.999),
             (2,3,.85,.75,.001),(3,2,.75,.85,.999),
             (2,4,.15,.95,.01),(4,2,.95,.15,.99),
             (3,3,.98,.96,.005),(3,3,.08,.12,.95),
             (4,4,.68,.72,.05),(4,4,.6,.8,.005),
             (2,5,.4,.8,.02),(5,2,.8,.4,.98)]
    for _ in range(random_count):
        n,k=rng.integers(2,7,size=2)
        la,lb=rng.uniform(.18,.96,size=2)
        phi=rng.uniform(.002,.998)
        designs.append((int(n),int(k),float(la),float(lb),float(phi)))
    checks=[]
    for index,design in enumerate(designs):
        n,k,la,lb,phi=design
        objective=RankOneLikelihood(make_population(*design))
        candidates=stationary_candidates(*design)
        minima=[c for c in candidates if c['local_minimum']]
        target=min(c['F'] for c in candidates)
        fits=[]
        for j in range(starts):
            loads=rng.uniform(-.85,.85,n+k)
            if j%3==0:
                loads=np.abs(loads)
            eta=rng.uniform(np.log(.15),np.log(1.3),n+k)
            fits.append(objective.optimize(np.r_[loads,eta]))
        best=min(fits,key=lambda r:r['F'])
        settled=[]
        for fit in fits:
            theta=fit['theta'].copy()
            if theta[:n+k].sum()<0:
                theta[:n+k]*=-1
            l=theta[:n+k]
            psi=np.exp(theta[n+k:])
            mineig=float(np.linalg.eigvalsh(objective.hessian(theta))[0])
            spread=float(max(np.ptp(l[:n]),np.ptp(l[n:]),
                             np.ptp(psi[:n]),np.ptp(psi[n:])))
            well_resolved=(fit['score']<5e-8 and mineig>1e-6)
            settled.append({'F':fit['F'],'score_inf':fit['score'],
                            'min_full_hessian_eigenvalue':mineig,
                            'within_block_range':spread,
                            'well_resolved':well_resolved,
                            'optimizer_success':fit['optimizer_success']})
        resolved=[r for r in settled if r['well_resolved']]
        worst_spread=max([r['within_block_range'] for r in resolved],default=0.0)
        minimum_error=abs(best['F']-target)
        all_not_below=all(r['F']>=target-2e-8 for r in fits)
        well_resolved_stationary=all(min(abs(r['F']-c['F']) for c in minima)<2e-8
                                     for r in resolved)
        passed=(minimum_error<2e-8 and all_not_below
                and worst_spread<2e-5 and well_resolved_stationary)
        check={'check':f'random_full_optimizer_{index:02d}',
               'status':'PASS' if passed else 'FAIL','classification':'Numerical only',
               'parameters':dict(zip(['m_A','m_B','lambda_A','lambda_B','phi'],design)),
               'tolerance':{'F_absolute':2e-8,'within_block_range':2e-5,
                            'score_inf_for_spread_check':5e-8,
                            'minimum_curvature_for_spread_check':1e-6},
               'actual':{'best_F':best['F'],'quartic_global_F':target,
                         'absolute_F_error':minimum_error,
                         'max_well_resolved_within_block_range':worst_spread,
                         'positive_stationary_count':len(candidates),
                         'positive_minimum_count':len(minima),
                         'resolved_count':len(resolved),'starts':starts,
                         'runs':settled}}
        checks.append(check)
        print(f"{check['status']} {check['check']}: F error {minimum_error:.3g}, "
              f"spread {worst_spread:.3g}, resolved {len(resolved)}/{starts}",flush=True)
    return checks


# Independent component: verification_extra.py

import mpmath as mp

import numpy as np

import sympy as sp

def general_mp_candidates(n,k,la,lb,phi):
    n,k,la,lb,phi=map(mp.mpf,(n,k,la,lb,phi))
    coefficients=quartic_coefficients(n,k,la,lb,phi)
    guesses=np.roots([float(z) for z in coefficients])
    rows=[]
    roots=[]
    for z in guesses:
        if z.real<=0 or abs(z.imag)>1e-7:
            continue
        root=mp.findroot(lambda v:mp.polyval(coefficients,v),mp.mpf(str(z.real)),
                         tol=mp.mpf('1e-65'),maxsteps=100)
        if root<=0 or any(abs(root-r)<mp.mpf('1e-40') for r in roots):
            continue
        roots.append(root)
        x,y,tau=reconstruct(n,k,la,lb,phi,root)
        alpha,beta=mp.sqrt(1-x),mp.sqrt(1-y)
        ba,bb=1+(n-1)*la*la,1+(k-1)*lb*lb
        residual=max(abs(ba*alpha/x+la*lb*phi*k*beta/y-tau*alpha),
                     abs(bb*beta/y+la*lb*phi*n*alpha/x-tau*beta))
        if residual>mp.mpf('1e-55'):
            raise ArithmeticError('General scalar root failed original scores')
        value=profile_value(n,k,la,lb,phi,root,logarithm=mp.log)
        curvature=mp.diff(lambda v:profile_value(n,k,la,lb,phi,v,logarithm=mp.log),root,2)
        share=share_a(n,k,la,lb,phi,root,square_root=mp.sqrt)
        rows.append({'v':root,'x':x,'y':y,'alpha':alpha,'beta':beta,'F':value,
                     'share_A':share,'curvature_G_vv':curvature,
                     'score_residual':residual})
    return sorted(rows,key=lambda r:r['v'])

def extension_symbolic_checks():
    r=sp.Rational
    R=sp.Matrix([[1,r(1,10),r(1,2),r(1,2)],
                 [r(1,10),1,r(1,2),r(1,2)],
                 [r(1,2),r(1,2),1,r(1,10)],
                 [r(1,2),r(1,2),r(1,10),1]])
    l=sp.Matrix([1,r(1,10),r(1,2),r(1,2)])
    S=l*l.T+sp.diag(0,r(99,100),r(3,4),r(3,4))
    K=S.inv()
    D=K-K*R*K
    expected={r(9,10):2,r(21,10):1,r(1,10):1}
    u,n,k,L,v=sp.symbols('u n k L v')
    x,y=1-u,1-v
    tau=1+n*u/x+k*v/y
    identity=tau*x*x-(1-L)*(1+u)-((L-u)*(1+u)+(n-2)*u*x+k*v*x*x/y)
    return {
        'E1_interchangeability_counterexample_spectrum':R.eigenvals()==expected,
        'E2_interchangeability_counterexample_inward_derivative':D[0,0]==r(2,165),
        'E3_interchangeability_counterexample_discrepancy':sp.trace(K*R)==4 and S.det()/R.det()==r(275,84),
        'E4_contrast_positive_identity':sp.cancel(identity)==0,
    }

def extra_checks(seed=20261008):
    checks=[]
    data={}
    def record(name,ok,value,tolerance,kind='Numerical only'):
        item={'check':name,'status':'PASS' if ok else 'FAIL',
              'classification':kind,'tolerance':tolerance,'actual':value}
        checks.append(item)
        print(item['status']+' '+name,flush=True)

    rng=np.random.default_rng(seed+1)
    objective=RankOneLikelihood(make_population(3,4,.63,.81,.37))
    z=np.r_[rng.uniform(-.4,.7,7),rng.uniform(-1.2,.3,7)]
    direction=rng.normal(size=14);direction/=np.linalg.norm(direction)
    h=2e-5
    f,g=objective.value_score(z)
    finite_score=(objective.value_score(z+h*direction)[0]-objective.value_score(z-h*direction)[0])/(2*h)
    finite_hessian=(objective.value_score(z+h*direction)[1]-objective.value_score(z-h*direction)[1])/(2*h)
    verror=abs(f-objective.dense_value(z))
    gerror=abs(finite_score-g@direction)
    herror=np.max(np.abs(finite_hessian-objective.hessian(z)@direction))
    record('dense_objective_and_analytic_derivatives',verror<1e-11 and gerror<1e-7 and herror<1e-7,
           {'value_error':verror,'directional_gradient_error':gerror,
            'directional_hessian_error':herror},
           {'dense_value':1e-11,'finite_difference_score':1e-7,'finite_difference_hessian':1e-7})

    boundary_rates=[]
    for n,k,la,lb,phi in [(2,2,.7,.7,.001),(2,5,.85,.3,.999),(3,2,.4,.9,.01)]:
        R=make_population(n,k,la,lb,phi)
        for j in [0,n]:
            load=R[:,j].copy();load[j]=1
            psi=1-load*load;psi[j]=0
            S=np.outer(load,load)+np.diag(psi)
            K=np.linalg.inv(S)
            rate=(K-K@R@K)[j,j]
            ix=[i for i in range(n+k) if i!=j]
            w=R[ix,j]/(1-R[ix,j]**2)
            conditional=R[np.ix_(ix,ix)]-np.outer(R[ix,j],R[ix,j])
            np.fill_diagonal(conditional,0)
            formula=-float(w@conditional@w)
            boundary_rates.append({'design':[n,k,la,lb,phi],'index':j,
                                   'direct_rate':rate,'formula_rate':formula,
                                   'absolute_error':abs(rate-formula)})
    record('Heywood_inward_derivatives_in_two_item_cases',
           all(r['direct_rate']<0 and r['absolute_error']<1e-10 for r in boundary_rates),
           boundary_rates,{'rate_error':1e-10,'strict_sign_threshold':0})

    mpa,mpb=mp.mpf('.85'),mp.mpf('.75')
    cost_A=-mp.log(1-mpa**4)
    cost_B=-2*mp.log(1-mpb**2)-mp.log(1+2*mpb**2)
    equal_LB=mp.findroot(lambda L:-2*mp.log(1-L)-mp.log(1+2*L)-cost_A,mp.mpf('.51'))
    data['endpoint_example']={'TC_A':cost_A/2,'TC_B':cost_B/2,
                              'kappa_A':2*mpa**2/(1-mpa**2),
                              'kappa_B':3*mpb**2/(1-mpb**2),
                              'equal_TC_lambda_B_for_lambda_A_0_85':mp.sqrt(equal_LB)}
    endpoint=RankOneLikelihood(make_population(2,3,float(mpa),float(mp.sqrt(equal_LB)),0))
    family=[]
    L=float(mpa**2)
    for t in [L,(L+np.sqrt(L))/2,np.sqrt(L),(1+np.sqrt(L))/2,1.0]:
        load=np.array([t,L/t,0,0,0])
        psi=1-load*load
        S=np.outer(load,load)+np.diag(psi)
        value=np.trace(np.linalg.solve(S,endpoint.R))+np.linalg.slogdet(S)[1]-endpoint.logdet-5
        family.append({'retained_loading_1':t,'F':float(value)})
    load=np.r_[0,0,np.repeat(float(mp.sqrt(equal_LB)),3)]
    S=np.outer(load,load)+np.diag(1-load*load)
    other=np.trace(np.linalg.solve(S,endpoint.R))+np.linalg.slogdet(S)[1]-endpoint.logdet-5
    data['equal_TC_endpoint_representations']=family
    record('unequal_size_equal_TC_endpoint_and_two_item_family',
           max(abs(r['F']-float(cost_A)) for r in family)<1e-12 and abs(other-float(cost_A))<1e-12,
           {'retained_A_family':family,'retained_B_F':float(other),
            'common_cost':cost_A},1e-12)

    designs=[]
    for la,lb in [('.68','.72'),('.6','.8')]:
        for phi in ['.05','.40894','.9']:
            rows=general_mp_candidates(4,4,la,lb,phi)
            designs.append({'m_A':4,'m_B':4,'lambda_A':la,'lambda_B':lb,
                            'phi':phi,'stationary_candidates':rows})
    size_rows=[]
    for phi in ['.01','.5','.75','.775','.9']:
        rows=general_mp_candidates(2,3,'.85','.75',phi)
        best=min(rows,key=lambda row:row['F'])
        size_rows.append({'phi':phi,'candidate_count':len(rows),'best':best})
    n,k=mp.mpf(2),mp.mpf(3)
    Q=mp.sqrt(k*(k-1)/(n*(n-1)))
    r=mp.sqrt((n-1)*(k-1)/(n*k))
    t=(Q*mpb**2-mpa**2)/((Q-1)*mpa**2*mpb**2)
    crossing=(1+r)*t-r
    crossrows=general_mp_candidates(n,k,mpa,mpb,crossing)
    limit=share_limit(n,k,mpa,mpb)
    near=general_mp_candidates(n,k,mpa,mpb,1-mp.mpf('1e-12'))[0]
    data.update({'equal_size_reference_designs':designs,'unequal_size_reference_design':size_rows,
                 'share_crossing':{'phi':crossing,'normalized_loading_squared':t,
                                   'candidates':crossrows,'limit_share_A':limit,
                                   'share_A_at_phi_1_minus_1e_minus_12':near['share_A']}})
    record('unequal_size_exact_crossing_original_scores',len(crossrows)==1 and
           abs(crossrows[0]['share_A']-mp.mpf('.5'))<mp.mpf('1e-55') and
           crossrows[0]['curvature_G_vv']>0,
           data['share_crossing'],{'share': '1e-55','score':'1e-55'})
    record('unequal_size_high_phi_share_limit',abs(near['share_A']-limit)<mp.mpf('1e-10'),
           {'computed_share':near['share_A'],'limiting_share':limit,
            'absolute_difference':abs(near['share_A']-limit)},'1e-10')
    local=stationary_candidates(4,4,.68,.72,.05)
    data['mild_phi_0_05_full_hessians']=[{
        'F':row['F'],'alpha':row['alpha'],'beta':row['beta'],
        'eigenvalues_in_loading_log_uniqueness_coordinates':row['hessian_eigenvalues']}
        for row in local]
    record('mild_phi_0_05_full_hessian_min_saddle_min',
           len(local)==3 and local[0]['hessian_eigenvalues'][0]>0 and
           np.count_nonzero(local[1]['hessian_eigenvalues']<0)==1 and
           local[2]['hessian_eigenvalues'][0]>0,
           data['mild_phi_0_05_full_hessians'],{'eigenvalue_zero_cutoff':1e-8})
    return checks,data


# Independent component: contrast_crosscheck.py

import json

import mpmath as mp

mp.mp.dps=65

def contrast_stationary(n,k,la,lb,phi):
    la,lb,phi=map(mp.mpf,(la,lb,phi)); ba=1+(n-1)*la*la;bb=1+(k-1)*lb*lb;c=la*lb*phi
    def score(al,be):
        x,y=1-al*al,1-be*be;tau=1+n*al*al/x+k*be*be/y
        return ba*al/x+c*k*be/y-tau*al,bb*be/y+c*n*al/x-tau*be
    for guess in [(la,lb),(phi*la,lb),(la,phi*lb),(la/2,lb/2)]:
        try:
            al,be=mp.findroot(score,guess,tol=mp.mpf('1e-55'),maxsteps=150)
            if 0<al<1 and 0<be<1 and min(al,be)>mp.mpf('1e-30'):
                return al,be
        except (ValueError,ZeroDivisionError):
            pass
    raise ValueError('No physical contrast_stationary point')

def contrast_direct_case(n,k,la,lb,phi):
    la,lb,phi=map(mp.mpf,(la,lb,phi));al,be=contrast_stationary(n,k,la,lb,phi)
    p=n+k;x,y=1-al*al,1-be*be;tau=1+n*al*al/x+k*be*be/y
    R=mp.matrix(p)
    for i in range(p):
        for j in range(p):
            R[i,j]=1 if i==j else (la*la if i<n and j<n else lb*lb if i>=n and j>=n else la*lb*phi)
    load=mp.matrix([al]*n+[be]*k);psi=[x]*n+[y]*k
    results=[]
    for block,offset,uni,loading,a in [('A',0,x,al,1-la*la),('B',n,y,be,1-lb*lb)]:
        contrast=mp.matrix(p,1);contrast[offset]=1/mp.sqrt(2);contrast[offset+1]=-1/mp.sqrt(2)
        def full_f(t,u):
            lv=load+t*contrast
            diagonal=mp.diag([psi[i]*mp.exp(u*contrast[i]) for i in range(p)])
            sigma=lv*lv.T+diagonal
            product=sigma**-1*R
            return sum(product[i,i] for i in range(p))+mp.log(mp.det(sigma))
        direct=mp.matrix([[mp.diff(lambda t:full_f(t,0),0,2),mp.diff(lambda u:mp.diff(lambda t:full_f(t,u),0),0)],
                          [mp.diff(lambda u:mp.diff(lambda t:full_f(t,u),0),0),mp.diff(lambda u:full_f(0,u),0,2)]])
        formula=mp.matrix([[2/uni-2*a/(tau*uni**2),2*a*loading/(tau*uni**2)],
                           [2*a*loading/(tau*uni**2),1-2*a*loading**2/(tau*uni**2)]])
        error=max(abs(direct[i,j]-formula[i,j]) for i in range(2) for j in range(2))
        assert error<mp.mpf('1e-48')
        results.append({'block':block,'error':mp.nstr(error,12),'direct':[[mp.nstr(direct[i,j],30) for j in range(2)] for i in range(2)],'determinant':mp.nstr(mp.det(formula),30)})
    return {'design':[n,k,str(la),str(lb),str(phi)],'alpha':mp.nstr(al,40),'beta':mp.nstr(be,40),'checks':results}


REFERENCE_BASELINE = {'lemma1_max_within_block_spread_60_designs': 1.2569678529243862e-08,
 'mirrored_threshold_[np.float64(0.7), np.float64(0.7), np.float64(0.7), np.float64(0.7)]': 0.40894039813188643,
 'mean_c_prediction_[np.float64(0.7), np.float64(0.7), np.float64(0.7), np.float64(0.7)]': 0.40894039735099336,
 'mirrored_threshold_[np.float64(0.6), np.float64(0.65), np.float64(0.75), np.float64(0.8)]': 0.4410748997855084,
 'mean_c_prediction_[np.float64(0.6), np.float64(0.65), np.float64(0.75), np.float64(0.8)]': 0.43354668262529106,
 'mirrored_threshold_[np.float64(0.5), np.float64(0.6), np.float64(0.7), np.float64(0.8)]': 0.41314149273933515,
 'mean_c_prediction_[np.float64(0.5), np.float64(0.6), np.float64(0.7), np.float64(0.8)]': 0.39833573259365895,
 'share_mild_.68_.72_phi0.40894': 0.774667148187578,
 'share_mild_.68_.72_phi0.9': 0.5725938545279485,
 'share_strong_.60_.80_phi0.40894': 0.9242738809944869,
 'share_strong_.60_.80_phi0.9': 0.8076929286279912,
 'mild_phi.05_F_local_vs_global': [1.2480239020621688, 0.9878350891063103],
 'fold_eps0.001': 0.385,
 'fold_eps0.0025': 0.365,
 'fold_eps0.005': 0.34,
 'fold_eps0.01': 0.301,
 'fold_eps0.02': 0.2415,
 'fold_eps0.04': 0.1555,
 'tilt_at_phistar_eps0.0005': 0.0836000411549962,
 'tilt_at_phistar_eps0.001': 0.1051130282685081,
 'tilt_at_phistar_eps0.002': 0.131997784319825,
 'tilt_at_phistar_eps0.004': 0.16542325766354926,
 'tilt_at_phistar_eps0.008': 0.20661948097659244,
 'tilt_at_phistar_eps0.016': 0.25662055049883925,
 'sizes_2x.85_3x.75_shareA_phi0.01': 0.0023000340539457926,
 'sizes_2x.85_3x.75_shareA_phi0.5': 0.18348636099358043,
 'sizes_2x.85_3x.75_shareA_phi0.75': 0.4964914399282228,
 'sizes_2x.85_3x.75_shareA_phi0.775': 0.5071170516130915,
 'sizes_2x.85_3x.75_shareA_phi0.9': 0.5366662618021859}


def audit_reference_numbers(reference,random_checks,cusp_data,mirror_data,extra_data):
    actual={}
    original_spread_key='lemma1_max_within_block_spread_60_designs'
    actual[original_spread_key]=max(r['actual']['max_well_resolved_within_block_range'] for r in random_checks)
    mirrorkeys=[key for key in reference if key.startswith('mirrored_threshold_')]
    meankeys=[key for key in reference if key.startswith('mean_c_prediction_')]
    for key,row in zip(mirrorkeys,mirror_data): actual[key]=float(row['threshold'])
    for key,row in zip(meankeys,mirror_data): actual[key]=float(row['mean_c_prediction'])
    for row in extra_data['equal_size_reference_designs']:
        label='mild_.68_.72' if row['lambda_A']=='.68' else 'strong_.60_.80'
        if row['phi'] in ['.40894','.9']:
            best=min(row['stationary_candidates'],key=lambda x:x['F'])
            actual[f"share_{label}_phi{float(row['phi'])}"]=float(max(best['share_A'],1-best['share_A']))
        elif row['phi']=='.05' and row['lambda_A']=='.68':
            rr=row['stationary_candidates']
            actual['mild_phi.05_F_local_vs_global']=[float(rr[0]['F']),float(rr[-1]['F'])]
    for row in cusp_data['folds_lambda_0_7_plus_minus_epsilon_m_4']:
        key='fold_eps'+str(float(row['epsilon']))
        if key in reference: actual[key]=float(row['phi_fold'])
    for row in cusp_data['tilts_at_symmetric_threshold']:
        actual['tilt_at_phistar_eps'+str(float(row['epsilon']))]=float(row['tilt'])
    for row in extra_data['unequal_size_reference_design']:
        actual['sizes_2x.85_3x.75_shareA_phi'+str(float(row['phi']))]=float(row['best']['share_A'])
    checks=[]
    for key,expected in reference.items():
        if key not in actual: raise KeyError('Missing independent reference comparison: '+key)
        observed=actual[key]
        error=float(np.max(np.abs(np.asarray(observed)-np.asarray(expected))))
        record={'check':'reference_'+key,'status':'PASS' if error<=1e-3 else 'FAIL',
                'classification':'Numerical only','tolerance':1e-3,
                'actual':observed,'reference':expected,'absolute_difference':error}
        if key==original_spread_key:
            record['interpretation']='Different random design experiments; both spreads are reported, not an exact replication of the reference sample.'
        checks.append(record)
        print(record['status']+' '+record['check']+': absolute difference '+format(error,'.4g'),flush=True)
    return checks


def json_ready(value):
    if isinstance(value,mp.mpf): return mp.nstr(value,65)
    if isinstance(value,np.ndarray): return [json_ready(x) for x in value.tolist()]
    if isinstance(value,np.generic): return value.item()
    if isinstance(value,dict): return {str(k):json_ready(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [json_ready(x) for x in value]
    if isinstance(value,sp.Basic): return str(value)
    return value


def main():
    import argparse,platform,time,sys
    import scipy
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed',type=int,default=20261008)
    parser.add_argument('--starts',type=int,default=10)
    parser.add_argument('--random-designs',type=int,default=40)
    parser.add_argument('--output',type=Path,default=Path('verification_results_unequal.json'))
    parser.add_argument('--reference-json',type=Path,default=None)
    args=parser.parse_args()
    started=time.monotonic()
    mp.mp.dps=80
    all_checks=[]
    random_results=run_random_checks(seed=args.seed,random_count=args.random_designs,starts=args.starts)
    all_checks.extend(random_results)
    exact=dict(scalar_symbolic_checks(),**extension_symbolic_checks())
    for name,passed in exact.items():
        record={'check':name,'status':'PASS' if passed else 'FAIL',
                'classification':'Proved with computer algebra','tolerance':'exact',
                'actual':bool(passed)}
        all_checks.append(record)
        print(record['status']+' '+name,flush=True)
    cusp_data=run_cusp_checks()
    for name,passed in cusp_data['checks'].items():
        exact_check=name.startswith('C') and name[1:2].isdigit()
        tolerance='exact' if exact_check else ('1e-6' if name.startswith('cusp_scaling') else '1e-55')
        record={'check':name,'status':'PASS' if passed else 'FAIL',
                'classification':'Proved with computer algebra' if exact_check else 'Numerical only',
                'tolerance':tolerance,'actual':bool(passed)}
        all_checks.append(record)
        print(record['status']+' '+name,flush=True)
    mirrored=[]
    for index,(loads,guess) in enumerate([([.7]*4,.409),([.6,.65,.75,.8],.441),([.5,.6,.7,.8],.413)]):
        row=mirrored_threshold(loads,guess,precision=60)
        mirrored.append(row)
        passed=mp.mpf(row['equation_residual'])<mp.mpf('1e-48')
        record={'check':f'mirrored_threshold_{index}_high_precision',
                'status':'PASS' if passed else 'FAIL','classification':'Numerical only',
                'tolerance':'1e-48','actual':row}
        all_checks.append(record)
        print(record['status']+' '+record['check'],flush=True)
    extra_results,extra_data=extra_checks(seed=args.seed)
    all_checks.extend(extra_results)
    direct_contrasts=[]
    with mp.workdps(65):
        for index,design in enumerate([(2,2,'.3','.8','.002'),(2,3,'.85','.75','.4'),
                                       (3,2,'.65','.9','.999'),(4,4,'.68','.72','.05'),
                                       (2,5,'.98','.21','.4')]):
            row=contrast_direct_case(*design)
            direct_contrasts.append(row)
            record={'check':f'direct_full_covariance_contrast_hessian_{index}',
                    'status':'PASS','classification':'Numerical only',
                    'tolerance':'1e-48','precision_decimal_digits':65,'actual':row}
            all_checks.append(record)
            print(record['status']+' '+record['check'],flush=True)
    extra_data['direct_full_covariance_contrast_checks']=direct_contrasts
    reference=REFERENCE_BASELINE if args.reference_json is None else json.loads(args.reference_json.read_text())
    reference_checks=audit_reference_numbers(reference,random_results,cusp_data,mirrored,extra_data)
    all_checks.extend(reference_checks)
    strongfold=next(row['phi_fold'] for row in cusp_data['folds_lambda_0_7_plus_minus_epsilon_m_4']
                    if abs(row['epsilon']-mp.mpf('.1'))<mp.mpf('1e-50'))
    prose_refinements=[
        {'claim':'Unequal-size equal-share crossing','paper_approximation':'.755',
         'correct_value':extra_data['share_crossing']['phi'],
         'absolute_difference':abs(extra_data['share_crossing']['phi']-mp.mpf('.755')),
         'reporting_threshold':'1e-3'},
        {'claim':'Fold for .6 versus .8','paper_approximation':'.02',
         'correct_value':strongfold,'absolute_difference':abs(strongfold-mp.mpf('.02')),
         'reporting_threshold':'1e-3'}]
    elapsed=time.monotonic()-started
    failures=[item['check'] for item in all_checks if item['status']!='PASS']
    payload={
        'title':'Independent verification: unequal blocks in a one-factor ML fit',
        'provenance':{'reference_program_reused':False,'reference_program_executed':False,
                      'reference_json_used_only_for_comparison':True,
                      'numerical_evidence_is_not_proof':True},
        'environment':{'python':platform.python_version(),'numpy':np.__version__,
                       'scipy':scipy.__version__,'sympy':sp.__version__,'mpmath':mp.__version__},
        'settings':{'seed':args.seed,'random_designs':args.random_designs,
                    'fixed_edge_designs':12,'starts_per_design':args.starts,
                    'full_optimizer_coordinates':'unrestricted loadings, log uniquenesses',
                    'log_uniqueness_bounds':[-24,4],
                    'high_precision_scalar_dps':80,'mirrored_threshold_dps':60},
        'summary':{'checks':len(all_checks),'passed':len(all_checks)-len(failures),
                   'failed':len(failures),'failed_checks':failures,
                   'full_optimizations':(12+args.random_designs)*args.starts,
                   'maximum_random_objective_error':max(r['actual']['absolute_F_error'] for r in random_results),
                   'maximum_random_within_block_range':max(r['actual']['max_well_resolved_within_block_range'] for r in random_results),
                   'reference_comparisons':len(reference_checks),
                   'reference_disagreements_above_1e_minus_3':[r['check'] for r in reference_checks if r['status']=='FAIL'],
                   'elapsed_seconds':elapsed},
        'checks':all_checks,'cusp_and_fold_data':cusp_data,
        'mirrored_threshold_data':mirrored,'additional_data':extra_data,
        'exact_unequal_example_certificate':unequal_example_certificate(),
        'paper_prose_refinements_above_1e_minus_3':prose_refinements,
        'scope_limits':['Random starts are not a proof of completeness or global optimality.',
                        'The scalar quartic classifies positive stationary points and all local minima; other saddles are not enumerated.',
                        'No general bound or improved approximation for mirrored heterogeneous thresholds is certified.']}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(json_ready(payload),indent=2)+'\n')
    print(f"TOTAL {len(all_checks)-len(failures)} PASS, {len(failures)} FAIL; {elapsed:.2f} seconds",flush=True)
    print('WROTE '+str(args.output.resolve()),flush=True)
    return 1 if failures else 0


if __name__=='__main__':
    raise SystemExit(main())
