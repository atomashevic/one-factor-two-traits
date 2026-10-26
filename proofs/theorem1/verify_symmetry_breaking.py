#!/usr/bin/env python3
"""Independent numerical VERIFICATION of the two-trait one-factor theorems.

This script is NOT a proof and does not use interval arithmetic. The accompanying
proof supplies the exact mathematics. Numerical optimization can miss minima.

Requirements: Python >= 3.10, numpy, scipy. No sympy or mpmath is needed.

Examples:
    python verify_symmetry_breaking.py --output verification_results.json
    python verify_symmetry_breaking.py --full-grid --starts 25 --jobs 4 \
        --seed 20261008 --output verification_results.json

The default runs 128 random algebraic/covariance checks, 24 independently sampled
populations with 8 random full-parameter ML starts each, 12 populations of full
ULS optimization, all named examples, the boundary/flat-family checks, and
high-precision finite-difference checks of the critical coefficients.

--full-grid selects an INDEPENDENT 144-population grid:
    m = 2, 3, 4;
    lambda = .4, .5, .6, .7, .8, .9;
    phi = 0, .05, .15, .30, .45, .60, .80, .95.
The task did not specify its original 144 populations individually. This grid
therefore must not be described as an exact rerun of an undocumented original.

ML optimization uses all 2p parameters (unrestricted l and log(Psi)), with no
block tying, diagonal-reproduction constraint, analytic starting values, or
analytic candidate in the pool of numerical starts. Half the starts are random
common-sign vectors and half have random mixed signs. log(Psi) is bounded to
[-20, 5], so the optimizer itself explores positive uniquenesses only; exact
zero-uniqueness cases are checked separately and the bounds are reported.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from decimal import Decimal, localcontext
import itertools
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

# Avoid multiplying BLAS thread counts when --jobs is used.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy import __version__ as scipy_version
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize

LOG_PSI_BOUNDS = (-20.0, 5.0)
GRID_M = [2, 3, 4]
GRID_LOADING = [.4, .5, .6, .7, .8, .9]
GRID_PHI = [0.0, .05, .15, .30, .45, .60, .80, .95]


def population(m: int, lam: float, phi: float) -> np.ndarray:
    L = lam * lam
    R = np.full((2 * m, 2 * m), L * phi)
    R[:m, :m] = L
    R[m:, m:] = L
    np.fill_diagonal(R, 1.0)
    return R


def covariance(l: np.ndarray, psi: np.ndarray) -> np.ndarray:
    return np.outer(l, l) + np.diag(psi)


def discrepancy(R: np.ndarray, Sigma: np.ndarray) -> float:
    """Direct SPD covariance calculation, independent of the closed forms."""
    chol = cho_factor(Sigma, lower=True, check_finite=False)
    inv_r = cho_solve(chol, R, check_finite=False)
    logdet_s = 2 * np.log(np.diag(chol[0])).sum()
    sign, logdet_r = np.linalg.slogdet(R)
    if sign <= 0:
        raise ValueError("R must be positive definite")
    return float(np.trace(inv_r) + logdet_s - logdet_r - len(R))


def ml_formula(m: int, lam: float, phi: float) -> dict:
    L = lam * lam
    a, n, N = 1 - L, m - 1, 2 * m - 1
    b = 1 + n * L
    t = abs(phi)
    critical = b / (m * (2 - L))
    q = m * L * t / b
    qstar = L / (2 - L)
    F0 = -n * math.log(a) - math.log(b)
    km = m * L * (1 - t) / a
    Fs = N * math.log1p(km / N) - math.log1p(km)
    if t < critical:
        S = L + (2 - L) * q * q
        # The factored radicand is more stable close to the threshold.
        gap = math.sqrt(max(0.0, (2 - L) * (qstar - q) * (1 - q)))
        common = math.sqrt(S + 2 * q)
        alpha, beta = (common + gap) / 2, (common - gap) / 2
        F = F0 + n * math.log1p(-q * q)
        branch = "asymmetric"
    else:
        alpha = beta = math.sqrt(L * (n + m * t) / N)
        F, branch = Fs, "symmetric"
    if phi < 0:
        beta = -beta
    l = np.r_[np.full(m, alpha), np.full(m, beta)]
    psi = 1 - l * l
    baseline = 2 * F0 - math.log1p(-q * q)
    return dict(m=m, lam=lam, phi=phi, threshold=critical, q=q, qstar=qstar,
                alpha=alpha, beta=beta, l=l, psi=psi, F=F, Fs=Fs,
                F0=F0, baseline=baseline, CFI=1 - F / baseline, branch=branch)


def uls_formula(m: int, lam: float, phi: float) -> dict:
    L, t = lam * lam, abs(phi)
    if t < 1 / (2 * m):
        root = math.sqrt(max(0.0, 1 - (2 * m * t)**2))
        alpha = math.sqrt(L * (1 + root) / 2)
        # Product formula avoids cancellation for small phi.
        beta = m * L * t / alpha
        val = m * (m - 1) * L * L / 2 - m*m*(m-1)*L*L*t*t
    else:
        alpha = beta = math.sqrt(L * (m - 1 + m * t) / (2 * m - 1))
        val = m*m*(m-1)*L*L*(1-t)**2/(2*m-1)
    if phi < 0:
        beta = -beta
    return dict(l=np.r_[np.full(m, alpha), np.full(m, beta)], F=val,
                alpha=alpha, beta=beta, threshold=1/(2*m))


def uls_value_gradient(l: np.ndarray, R: np.ndarray):
    E = np.outer(l, l) - R
    np.fill_diagonal(E, 0.0)
    return float(.5 * np.sum(E * E)), 2 * E @ l


def uls_hessian(l: np.ndarray, R: np.ndarray) -> np.ndarray:
    H = 4 * np.outer(l, l) - 2 * R
    np.fill_diagonal(H, 2 * (np.dot(l, l) - l * l))
    return H


class Checks:
    def __init__(self):
        self.metrics = {}
        self.failures = []
        self.count = 0

    def error(self, name, err, tolerance, context=None):
        err = float(abs(err))
        self.count += 1
        old = self.metrics.get(name)
        if old is None or err > old["maximum_error"]:
            self.metrics[name] = dict(maximum_error=err, tolerance=float(tolerance),
                                      worst_case=context)
        if not math.isfinite(err) or err > tolerance:
            self.failures.append(dict(check=name, error=err,
                                      tolerance=float(tolerance), context=context))

    def close(self, name, observed, expected, tolerance=1e-10, context=None):
        # Absolute error for quantities <= 1, relative error otherwise.
        err = abs(float(observed) - float(expected)) / max(1.0, abs(float(expected)))
        self.error(name, err, tolerance, context)

    def truth(self, name, condition, context=None):
        self.error(name, 0.0 if condition else 1.0, 0.0, context)

    def result(self):
        return dict(number_of_assertions=self.count, metrics=self.metrics,
                    passed=not self.failures, failures=self.failures)


def random_identity_checks(rng, count, checks):
    """Cross-check exact identities using direct matrices and eigenvalues."""
    cases = []
    for index in range(count):
        m = int(rng.integers(2, 13))
        lam = float(rng.uniform(.08, .98))
        phi = float(rng.uniform(0, .98))
        context = dict(index=index, m=m, lam=lam, phi=phi)
        cases.append((m, lam, phi))
        L, a, n, N = lam*lam, 1-lam*lam, m-1, 2*m-1
        b, h = 1+n*L, m*L*phi
        R = population(m, lam, phi)
        fit = ml_formula(m, lam, phi)
        l, psi = fit["l"], fit["psi"]
        direct = discrepancy(R, covariance(l, psi))
        checks.close("closed_minimum_vs_direct_covariance", direct, fit["F"],
                     3e-10, context)
        checks.close("population_determinant", np.linalg.slogdet(R)[1],
                     2*n*math.log(a)+math.log(b*b-h*h), 3e-11, context)
        u = l / psi
        gamma = float(l @ u)
        residual = R @ u - (1 + gamma) * l
        checks.error("Lawley_score_at_closed_solution",
                     np.max(np.abs(residual))/max(1, np.max(np.abs(R@u))),
                     3e-10, context)
        checks.error("diagonal_reproduction_at_closed_solution",
                     np.max(np.abs(np.diag(covariance(l, psi))-1)), 2e-13, context)
        kp, km = m*L*(1+phi)/a, m*L*(1-phi)/a
        lhs = N*(kp-km)-km*(2+kp+km)
        rhs = 2*m*L/a*((N+1+m*L/a)*phi-(1+m*L/a))
        checks.close("threshold_equivalence_algebra", lhs, rhs, 2e-10, context)
        checks.truth("threshold_equivalence_sign",
                     (lhs >= 0) == (phi >= fit["threshold"]), context)

        # Arbitrary positive block uniquenesses; not constrained to a solution.
        x, y = np.exp(rng.uniform(-2.5, 2.5, 2))
        ps = np.r_[np.full(m, x), np.full(m, y)]
        W = R / np.sqrt(np.outer(ps, ps))
        evals, evecs = np.linalg.eigh(W)
        M = np.array([[b/x, h/math.sqrt(x*y)], [h/math.sqrt(x*y), b/y]])
        reduced = np.sort(np.r_[np.full(n, a/x), np.full(n, a/y),
                                np.linalg.eigvalsh(M)])
        checks.error("block_eigenvalue_formula",
                     np.max(np.abs(evals-reduced))/max(1, evals[-1]), 3e-12, context)
        checks.close("leading_eigenvalue_from_mean_block", evals[-1],
                     np.linalg.eigvalsh(M)[-1], 3e-12, context)
        white_l = evecs[:, -1] * math.sqrt(max(0, evals[-1]-1))
        fixed_l = np.sqrt(ps) * white_l
        gs = evals - np.log(evals) - 1
        concentrated = float(gs[:-1].sum() if evals[-1] >= 1 else gs.sum())
        checks.close("fixed_Psi_concentration", discrepancy(R,covariance(fixed_l,ps)),
                     concentrated, 5e-10, context)

        d = float(rng.uniform(-2, 2))
        rr = math.sqrt(b*b*math.sinh(d)**2+h*h)
        tm, tp = b*math.cosh(d)-rr, b*math.cosh(d)+rr
        Q = 2*n*a*math.cosh(d)+tm
        scale = Q/N
        checks.truth("profile_scale_is_active", tp/scale > 1, context)
        H = N*math.log(Q/N)-2*n*math.log(a)-math.log(tm)
        w = rr/(b*math.cosh(d))
        k = 2*n*a/b
        H_w = (n*math.log(b*b-h*h)-2*n*math.log(a)-N*math.log(N)
               +N*math.log(k+1-w)-m*math.log(1-w)-n*math.log(1+w))
        checks.close("profile_coordinate_identity", H, H_w, 2e-10, context)
        analytic = (N*k+2*n)*(w-fit["qstar"])/((k+1-w)*(1-w*w))
        eps = 1e-24
        wc = complex(w, eps)
        numerical = (N*np.log(k+1-wc)-m*np.log(1-wc)-n*np.log(1+wc)).imag/eps
        checks.close("profile_derivative_complex_step", numerical, analytic,
                     3e-10, context)

        # Root-elimination identity used by the within-block reduction lemma.
        tx, ty = rng.uniform(.02, .98, 2)
        px, py = 1-tx*tx, 1-ty*ty
        T = a*(1+tx*ty)/(px*py)
        gamma_rest = T-1-tx*tx/px-ty*ty/py
        checks.close("mixed_root_elimination_identity", gamma_rest*px*py,
                     (1+tx*ty)*(tx*ty-L), 3e-12, context)

        # Reflection on an arbitrary fitted covariance, and on the minimizer.
        signs = np.r_[np.ones(m), -np.ones(m)]
        lr = rng.normal(0, .5, 2*m)
        pr = np.exp(rng.uniform(-1, 1, 2*m))
        Rneg = population(m, lam, -phi)
        checks.close("reflection_arbitrary_covariance",
                     discrepancy(Rneg,covariance(lr,pr)),
                     discrepancy(R,covariance(signs*lr,pr)), 3e-11, context)
        reflected = ml_formula(m, lam, -phi)
        checks.close("reflection_optimal_covariance",
                     discrepancy(Rneg,covariance(reflected["l"],reflected["psi"])),
                     fit["F"], 3e-10, context)

        # C2, with residual variance explicitly included before standardizing.
        primary, secondary = rng.normal(0, .6, 2)
        rho, tau2 = float(rng.uniform(-.95,.95)), float(rng.uniform(.2,2))
        trait = np.array([[1, rho], [rho, 1]])
        loads = np.vstack([np.tile([primary,secondary],(m,1)),
                           np.tile([secondary,primary],(m,1))])
        signal = loads @ trait @ loads.T
        v = primary**2+secondary**2+2*primary*secondary*rho
        wsignal = 2*primary*secondary+rho*(primary**2+secondary**2)
        raw = signal + tau2*np.eye(2*m)
        normalized = raw/(tau2+v)
        effective = population(m, math.sqrt(v/(tau2+v)), wsignal/v)
        checks.error("cross_loading_standardization_general_tau",
                     np.max(np.abs(normalized-effective)), 3e-12, context)
        raw_unit_error = signal + np.eye(2*m)
        unit_error_eff = population(m,math.sqrt(v/(1+v)),wsignal/v)
        checks.error("cross_loading_formula_with_unit_residual",
                     np.max(np.abs(raw_unit_error/(1+v)-unit_error_eff)),
                     3e-12, context)
        multiplier = math.sqrt(.6/v)
        unit_item_signal = signal*multiplier**2
        already_unit = unit_item_signal+.4*np.eye(2*m)
        already_unit_eff = population(m,math.sqrt(.6),wsignal/v)
        checks.error("cross_loading_formula_with_unit_item_variance",
                     np.max(np.abs(already_unit-already_unit_eff)), 3e-12, context)

        uf = uls_formula(m,lam,phi)
        checks.close("ULS_closed_minimum_vs_direct", uls_value_gradient(uf["l"],R)[0],
                     uf["F"], 3e-11, context)
        z = L*(n+m*phi)/N
        eig = np.linalg.eigvalsh(uls_hessian(np.full(2*m,math.sqrt(z)),R))
        expected = np.sort(np.r_[np.full(2*n,2*((2*m-3)*z+L)),
                        4*L*(n+m*phi), 4*L*n*(2*m*phi-1)/N])
        checks.error("ULS_Hessian_eigenvalues",np.max(np.abs(eig-expected)),
                     3e-10, context)
        rand_l = rng.uniform(0,1,2*m)
        xblock, yblock = rand_l[:m], rand_l[m:]
        SA, SB = xblock.sum(), yblock.sum()
        QA, QB = float(xblock@xblock), float(yblock@yblock)
        K, Cm = m*(m-1), m*(m-1)/2
        aa, bb = (SA*SA-QA)/K, (SB*SB-QB)/K
        relaxed = (Cm*((L-aa)**2+(L-bb)**2)+QA*QB
                   -2*L*phi*math.sqrt((QA+K*aa)*(QB+K*bb))+m*m*L*L*phi*phi)
        direct_g = uls_value_gradient(rand_l,R)[0]
        checks.error("ULS_Jensen_lower_bound_violation",max(0,relaxed-direct_g),
                     5e-11, context)
    return cases


def named_examples(checks):
    rows = []
    critical = ml_formula(2,.7,0)["threshold"]
    specs = [(2,.7,0,.274568,.412924), (2,.7,.2,.257113,.310047),
             (2,.7,.473,None,None), (2,.7,critical-.02,.172604,.172889),
             (2,.7,.513,None,None), (4,.8,.3,1.778314,1.972),
             (4,.8,.6,1.040351,None), (4,.8,0,1.9934,None)]
    for m,lam,phi,expected,expected_s in specs:
        f = ml_formula(m,lam,phi)
        direct = discrepancy(population(m,lam,phi),covariance(f["l"],f["psi"]))
        row = {key:f[key] for key in ["m","lam","phi","threshold","alpha","beta",
                                      "F","Fs","CFI","branch"]}
        row["direct_covariance_F"] = direct
        row["quoted_rounded_minimum"] = expected
        row["quoted_rounded_symmetric"] = expected_s
        if m == 2 and phi == .473:
            row["note"] = "Literal .473; task's quoted values instead use phi_star-.02."
        if m == 2 and abs(phi-(critical-.02)) < 1e-15:
            row["note"] = "Unrounded phi_star-.02 reproduces task's .172604/.172889."
        checks.close("named_direct_F",direct,f["F"],3e-12,row)
        if expected is not None:
            tolerance = 5e-5 if expected == 1.9934 else 5e-7
            checks.error("quoted_rounded_minimum",abs(f["F"]-expected),tolerance,row)
        if expected_s is not None:
            tolerance = 5e-4 if expected_s == 1.972 else 5e-7
            checks.error("quoted_rounded_symmetric",abs(f["Fs"]-expected_s),tolerance,row)
        rows.append(row)
    return rows


def boundary_checks(rng, checks):
    rows = []
    specs = [(2,.7,0),(3,.7,0),(4,.8,0),(2,.7,.2),(4,.8,.3)]
    specs += [(int(rng.integers(2,9)),float(rng.uniform(.2,.95)),
               float(rng.uniform(.02,.95))) for _ in range(12)]
    for m,lam,phi in specs:
        R = population(m,lam,phi)
        l = R[:,0].copy()
        psi = 1-l*l
        S = covariance(l,psi)
        F = discrepancy(R,S)
        star_expected = np.log(psi[1:]).sum()-np.linalg.slogdet(R)[1]
        inv = np.linalg.inv(S)
        derivative = float((inv-inv@R@inv)[0,0])
        weights = R[1:,0]/(1-R[1:,0]**2)
        delta = R[1:,1:]-S[1:,1:]
        derivative_sum = float(-weights@delta@weights)
        step = 1e-5*float(np.linalg.eigvalsh(S)[0])
        E = np.zeros_like(S); E[0,0] = step
        fd = (-3*F+4*discrepancy(R,S+E)-discrepancy(R,S+2*E))/(2*step)
        row = dict(m=m,lam=lam,phi=phi,star_F=F,
                   inward_derivative=derivative,finite_difference=fd,step=step)
        checks.close("boundary_star_conditional_value",F,star_expected,3e-11,row)
        checks.close("boundary_inward_derivative_identity",derivative,
                     derivative_sum,3e-10,row)
        checks.close("boundary_inward_derivative_finite_difference",fd,
                     derivative_sum,4e-5,row)
        if m == 2 and phi == 0:
            checks.error("boundary_derivative_exception",abs(derivative),2e-11,row)
        else:
            checks.truth("boundary_derivative_strictly_negative",derivative<0,row)
        rows.append(row)
    family = []
    for lam in [.2,.7,.95]:
        L = lam*lam
        R = population(2,lam,0)
        for t in np.linspace(L,1,9):
            l = np.array([t,L/t,0,0])
            psi = 1-l*l
            psi[np.abs(psi)<1e-14] = 0
            S = covariance(l,psi)
            target = np.eye(4); target[:2,:2] = R[:2,:2]
            value = discrepancy(R,S)
            expected = -math.log1p(-L*L)
            row = dict(lam=lam,t=float(t),F=value,
                       minimum_uniqueness=float(psi.min()),
                       minimum_covariance_eigenvalue=float(np.linalg.eigvalsh(S)[0]))
            checks.error("m2_phi0_family_covariance",np.max(np.abs(S-target)),2e-12,row)
            checks.close("m2_phi0_family_value",value,expected,3e-11,row)
            checks.truth("m2_phi0_family_positive_definite",row["minimum_covariance_eigenvalue"]>0,row)
            family.append(row)
    return dict(star_checks=rows,exceptional_family_including_endpoints=family)


def critical_checks(checks):
    rows = []
    for mint, lamstr in [(2,".3"),(2,".7"),(4,".8"),(5,".2"),(10,".95")]:
        with localcontext() as ctx:
            ctx.prec = 80
            D = Decimal
            m, L = D(mint), D(lamstr)**2
            a,n,N = 1-L,m-1,2*m-1
            b = 1+n*L
            qstar = L/(2-L)
            critical = b/(m*(2-L))
            h = b*qstar
            def profile(d):
                co,si = (d.exp()+(-d).exp())/2,(d.exp()-(-d).exp())/2
                tm = b*co-(b*b*si*si+h*h).sqrt()
                Q = 2*n*a*co+tm
                return N*(Q/N).ln()-2*n*a.ln()-tm.ln()
            d = D("1e-8")
            K4numerical = (profile(d)-profile(D(0)))/(d**4)
            K4 = n*m*(2-L)**2/(2*N*L*L)
            k,A = 2*n*a/b,N*(2*n*a/b)+2*n
            shift = D("1e-18")
            slope = m*L/b
            q = qstar+slope*shift
            K2numerical = A*(q-qstar)/(2*q*(k+1-q))/shift
            K2 = n*m*m*(2-L)**3/(2*a*N*b)
            c = m*L/a
            Fssecond = -N*c*c/(N+c*(1-critical))**2+c*c/(1+c*(1-critical))**2
            Fasecond = -2*n*slope*slope*(1+qstar*qstar)/(1-qstar*qstar)**2
            J = n*m**3*L*L*(2-L)**4/(4*a*a*N*b*b)
            row = dict(m=mint,lam=float(lamstr),threshold=float(critical),
                       K4_exact=float(K4),K4_high_precision=float(K4numerical),
                       K2_exact=float(K2),K2_high_precision=float(K2numerical),
                       phi_second_derivative_jump=float(J),
                       phi_second_derivative_jump_direct=float(Fssecond-Fasecond))
            checks.close("Landau_quartic_high_precision",K4numerical,K4,2e-12,row)
            checks.close("Landau_quadratic_slope",K2numerical,K2,2e-12,row)
            checks.close("Fmin_C1_not_C2_second_derivative_jump",Fssecond-Fasecond,J,2e-12,row)
            # Verify matching first derivatives independently.
            Fsp = -N*c/(N+c*(1-critical))+c/(1+c*(1-critical))
            Fap = -2*n*slope*qstar/(1-qstar*qstar)
            checks.close("Fmin_matching_first_derivatives",Fsp,Fap,2e-12,row)
            phi = float(critical)*(1-1e-6)
            fit = ml_formula(mint,float(lamstr),phi)
            ratio = abs(fit["alpha"]-fit["beta"])/math.sqrt(float(critical)-phi)
            C = math.sqrt(float(2*m*a*L/b))
            row.update(C_exact=C,C_near_threshold=ratio)
            checks.close("loading_pitchfork_coefficient",ratio,C,5e-6,row)
            # Smooth antisymmetric diagonal-reproducing path, including phi=0.
            Erows = []
            for phi in [D(0),critical*D(".7"),critical*D("1.1")]:
                z = L*(n+m*phi)/N
                r,U,V = 1-z,b+m*L*phi,b-m*L*phi
                t,eps = z.sqrt(),D("1e-12")
                def path(e):
                    x,y = 1-(t+e)**2,1-(t-e)**2
                    pp,mm,off = U-e*e,r+N*e*e,2*n*t*e
                    det = pp*mm-off*off
                    return n*(a*(1/x+1/y)+(x*y).ln())+(U*mm+V*pp)/det+det.ln()
                Enumeric = (path(eps)-path(D(0)))/(eps*eps)
                Eexact = 4*n*m*(1+z)*((2-L)*z-L)/(r**3*(1+N*z))
                erow = dict(phi=float(phi),E_exact=float(Eexact),E_numeric=float(Enumeric))
                checks.close("smooth_antisymmetric_ML_path",Enumeric,Eexact,3e-12,erow)
                checks.truth("smooth_path_stability_sign",(Eexact<0)==(phi<critical),erow)
                Erows.append(erow)
            row["smooth_saddle_path"] = Erows
            rows.append(row)
    # Requested ULS loading/length checks explicitly surround the threshold.
    uls_rows = []
    for m,lam in itertools.product([2,3,5],[.5,.8]):
        threshold = 1/(2*m)
        for phi in [threshold-.01,threshold,threshold+.01]:
            L=lam*lam; N=2*m-1
            z = L*(m-1+m*phi)/N
            H = uls_hessian(np.full(2*m,math.sqrt(z)),population(m,lam,phi))
            contrast = np.r_[np.ones(m),-np.ones(m)]/math.sqrt(2*m)
            observed = float(contrast@H@contrast)
            expected = 4*L*(m-1)*(2*m*phi-1)/N
            row = dict(m=m,lam=lam,phi=phi,threshold=threshold,
                       contrast_curvature=observed,expected=expected)
            checks.close("named_ULS_stability_threshold",observed,expected,2e-12,row)
            uls_rows.append(row)
    return dict(ML_pitchfork=rows,ULS_thresholds=uls_rows)


def ml_objective_gradient(params, R, logdet_r):
    p = len(R)
    l, psi = params[:p], np.exp(params[p:])
    S = covariance(l,psi)
    try:
        chol = cho_factor(S,lower=True,check_finite=False)
        inv = cho_solve(chol,np.eye(p),check_finite=False)
        value = np.sum(inv*R)+2*np.log(np.diag(chol[0])).sum()-logdet_r-p
        score = inv-inv@R@inv
        grad = np.r_[2*score@l,psi*np.diag(score)]
        if not np.isfinite(value) or not np.all(np.isfinite(grad)):
            raise FloatingPointError
        return float(value),grad
    except (np.linalg.LinAlgError,ValueError,FloatingPointError):
        return 1e100,np.zeros(2*p)


def optimize_population(task):
    m,lam,phi,starts,seed = task
    rng = np.random.default_rng(seed)
    R = population(m,lam,phi); p=2*m
    logdet_r = np.linalg.slogdet(R)[1]
    best = None
    start_records = []
    for j in range(starts):
        if j%2 == 0:
            l0 = rng.uniform(.03,.95,p)*rng.choice([-1,1])
        else:
            l0 = rng.normal(0,.55,p)
        eta0 = rng.uniform(-1.6,.5,p)
        x0 = np.r_[l0,eta0]
        opt = minimize(ml_objective_gradient,x0,args=(R,logdet_r),jac=True,
                       method="L-BFGS-B",bounds=[(None,None)]*p+[LOG_PSI_BOUNDS]*p,
                       options=dict(maxiter=2200,maxls=45,ftol=2e-14,gtol=2e-9))
        start_records.append(dict(F=float(opt.fun),success=bool(opt.success),
                                  iterations=int(opt.nit),evaluations=int(opt.nfev),
                                  gradient_inf=float(np.max(np.abs(opt.jac)))))
        if best is None or opt.fun < best.fun:
            best=opt
    l,psi = best.x[:p].copy(),np.exp(best.x[p:])
    if l.sum()<0: l=-l
    f=ml_formula(m,lam,phi)
    S=covariance(l,psi)
    target=covariance(f["l"],f["psi"])
    swap=np.r_[np.arange(m,2*m),np.arange(m)]
    covariance_error=min(np.max(np.abs(S-target)),np.max(np.abs(S-target[np.ix_(swap,swap)])))
    u=l/psi; gamma=float(l@u)
    score=R@u-(1+gamma)*l
    exceptional=m==2 and phi==0
    return dict(m=m,lam=lam,phi=phi,starts=starts,seed=int(seed),
                numerical_F=float(best.fun),closed_F=f["F"],
                absolute_F_error=abs(float(best.fun)-f["F"]),
                normalized_F_error=abs(float(best.fun)-f["F"])/max(1,abs(f["F"])),
                covariance_error=float(covariance_error),
                diagonal_reproduction_error=float(np.max(np.abs(np.diag(S)-1))),
                stationarity_normalized_error=float(np.max(np.abs(score))/max(1,np.max(np.abs(R@u)))),
                gradient_inf=float(np.max(np.abs(best.jac))),
                loading_block_nonuniformity=float(max(np.ptp(l[:m]),np.ptp(l[m:]))),
                psi_block_nonuniformity=float(max(np.ptp(psi[:m]),np.ptp(psi[m:]))),
                sign_violation=float(max(0,-np.min(l))) if phi>0 else 0.0,
                exceptional_flat_family=exceptional,
                lower_logpsi_bound_distance=float(np.min(best.x[p:])-LOG_PSI_BOUNDS[0]),
                upper_logpsi_bound_distance=float(LOG_PSI_BOUNDS[1]-np.max(best.x[p:])),
                l=l.tolist(),psi=psi.tolist(),start_results=start_records)


def optimization_checks(populations,starts,seed,jobs,checks):
    tasks=[(m,lam,phi,starts,seed+104729*(i+1)) for i,(m,lam,phi) in enumerate(populations)]
    rows=[]
    def add(row):
        context={key:row[key] for key in ["m","lam","phi","starts"]}
        checks.error("full_ML_optimization_objective",row["normalized_F_error"],3e-7,context)
        checks.error("full_ML_optimization_covariance",row["covariance_error"],8e-4,context)
        checks.error("full_ML_diagonal_reproduction",row["diagonal_reproduction_error"],8e-5,context)
        checks.error("full_ML_stationarity",row["stationarity_normalized_error"],8e-5,context)
        checks.error("full_ML_common_sign",row["sign_violation"],5e-5,context)
        if not row["exceptional_flat_family"]:
            checks.error("full_ML_block_uniformity",row["loading_block_nonuniformity"],8e-4,context)
            checks.error("full_ML_uniqueness_block_uniformity",row["psi_block_nonuniformity"],8e-4,context)
        checks.truth("full_ML_upper_Psi_bound_inactive",row["upper_logpsi_bound_distance"]>1,context)
        # A lower-bound approach is allowed for the known m=2 endpoint family.
        if not row["exceptional_flat_family"]:
            checks.truth("full_ML_lower_Psi_bound_inactive",row["lower_logpsi_bound_distance"]>1,context)
        rows.append(row)
        if len(rows)%4 == 0 or len(rows)==len(tasks):
            print(f"ML: {len(rows)}/{len(tasks)} populations; "
                  f"largest F error {max(r['absolute_F_error'] for r in rows):.3g}",flush=True)
    if jobs>1:
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            for row in pool.map(optimize_population,tasks): add(row)
    else:
        for task in tasks: add(optimize_population(task))
    return dict(populations=len(rows),starts_per_population=starts,
                total_optimizations=len(rows)*starts,
                all_starts_random=True,analytic_start_used=False,
                parameter_dimension="2p = 4m",log_uniqueness_bounds=list(LOG_PSI_BOUNDS),
                minimum_allowed_uniqueness=math.exp(LOG_PSI_BOUNDS[0]),
                maximum_objective_error=max(r["absolute_F_error"] for r in rows),
                maximum_covariance_error=max(r["covariance_error"] for r in rows),rows=rows)


def uls_optimization_checks(populations,starts,seed,checks):
    rng=np.random.default_rng(seed)
    rows=[]
    for m,lam,phi in populations:
        R=population(m,lam,phi); p=2*m
        best=None
        for j in range(starts):
            l0=rng.uniform(.02,.95,p) if j%2==0 else rng.uniform(-1,1,p)
            opt=minimize(uls_value_gradient,l0,args=(R,),jac=True,method="L-BFGS-B",
                         bounds=[(-1,1)]*p,
                         options=dict(maxiter=1800,maxls=40,ftol=2e-14,gtol=2e-9))
            if best is None or opt.fun<best.fun:best=opt
        expected=uls_formula(m,lam,phi)
        row=dict(m=m,lam=lam,phi=phi,starts=starts,numerical_G=float(best.fun),
                 closed_G=expected["F"],absolute_error=abs(best.fun-expected["F"]),
                 l=best.x.tolist())
        checks.close("full_ULS_optimization",best.fun,expected["F"],3e-8,row)
        rows.append(row)
    return dict(populations=len(rows),starts_per_population=starts,
                all_starts_random=True,rows=rows)


def gradient_checks(rng,checks):
    for m in [2,3,4]:
        R=population(m,.73,.31); p=2*m
        x=np.r_[rng.normal(0,.4,p),rng.uniform(-1,.5,p)]
        value,grad=ml_objective_gradient(x,R,np.linalg.slogdet(R)[1])
        d=rng.normal(size=2*p);d/=np.linalg.norm(d)
        eps=1e-5
        fp=ml_objective_gradient(x+eps*d,R,np.linalg.slogdet(R)[1])[0]
        fm=ml_objective_gradient(x-eps*d,R,np.linalg.slogdet(R)[1])[0]
        checks.close("ML_optimizer_gradient_finite_difference",(fp-fm)/(2*eps),grad@d,3e-8,dict(m=m))
        l=x[:p]; dd=d[:p]
        _,gu=uls_value_gradient(l,R)
        fp=uls_value_gradient(l+eps*dd,R)[0];fm=uls_value_gradient(l-eps*dd,R)[0]
        checks.close("ULS_optimizer_gradient_finite_difference",(fp-fm)/(2*eps),gu@dd,3e-8,dict(m=m))


def parse_args():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--full-grid",action="store_true")
    parser.add_argument("--starts",type=int,default=8,help="Random ML starts per population (default 8).")
    parser.add_argument("--starts25",action="store_const",const=25,dest="starts",
                        help="Shortcut for --starts 25.")
    parser.add_argument("--seed",type=int,default=20261008)
    parser.add_argument("--random-checks",type=int,default=128,help="At least 100 random identity populations.")
    parser.add_argument("--random-populations",type=int,default=24)
    parser.add_argument("--uls-populations",type=int,default=12)
    parser.add_argument("--uls-starts",type=int,default=8)
    parser.add_argument("--jobs",type=int,default=1)
    parser.add_argument("--output",type=Path,default=Path("verification_results.json"))
    args=parser.parse_args()
    if args.random_checks<100:parser.error("--random-checks must be at least 100")
    if min(args.starts,args.random_populations,args.uls_populations,args.uls_starts,args.jobs)<1:
        parser.error("All counts must be positive")
    return args


def main():
    args=parse_args();started=time.time()
    rng=np.random.default_rng(args.seed)
    checks=Checks()
    print("Numerical verification only; this script is not a proof.",flush=True)
    gradient_checks(rng,checks)
    named=named_examples(checks)
    identity_cases=random_identity_checks(rng,args.random_checks,checks)
    boundary=boundary_checks(rng,checks)
    critical=critical_checks(checks)
    if args.full_grid:
        populations=list(itertools.product(GRID_M,GRID_LOADING,GRID_PHI))
        selection=dict(kind="independent_144_population_grid",m=GRID_M,lam=GRID_LOADING,
                       phi=GRID_PHI,note="Original task grid not specified individually; this grid is independent.")
    else:
        populations=[]
        for i in range(args.random_populations):
            m=int(rng.integers(2,5));lam=float(rng.uniform(.35,.92))
            phi=0.0 if i%8==0 else float(rng.uniform(.01,.95))
            populations.append((m,lam,phi))
        selection=dict(kind="random_populations",count=len(populations),
                       m="2..4",lam="uniform .35..92",phi="uniform .01..95; every eighth phi=0")
    full_ml=optimization_checks(populations,args.starts,args.seed,args.jobs,checks)
    uls_pops=[(2,.5,0),(2,.8,.1),(3,.5,1/6),(5,.8,.08)]
    # The remaining ULS populations are independent of the analytic identity pool.
    while len(uls_pops)<args.uls_populations:
        uls_pops.append((int(rng.integers(2,7)),float(rng.uniform(.35,.92)),float(rng.uniform(.01,.95))))
    full_uls=uls_optimization_checks(uls_pops[:args.uls_populations],args.uls_starts,args.seed+900001,checks)
    results=dict(label="NUMERICAL VERIFICATION, NOT PROOF",passed=not checks.failures,
        elapsed_seconds=time.time()-started,seed=args.seed,
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy_version),
        limitations=["Floating-point and finite-difference checks are not certified proofs.",
                     "Random local optimization can miss a better minimum.",
                     "ML optimizer bounds log(Psi) to [-20,5]; exact boundary points are checked separately.",
                     "No analytic minimizer or block-uniform parameterization initializes the unrestricted optimizer.",
                     "The full-grid option specifies an independent 144-population grid; the original individual grid points were not supplied."],
        configuration=dict(random_identity_populations=args.random_checks,
                           population_selection=selection,jobs=args.jobs),
        named_sanity_checks=named,boundary_and_exception_checks=boundary,
        critical_and_stability_checks=critical,full_ML_optimization=full_ml,
        full_ULS_optimization=full_uls,checks=checks.result())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(f"{'PASS' if results['passed'] else 'FAIL'}: {checks.count} checks; "
          f"{len(checks.failures)} failures; {results['elapsed_seconds']:.1f} s",flush=True)
    print(f"Results: {args.output.resolve()}",flush=True)
    if checks.failures:
        for fail in checks.failures[:12]:print(json.dumps(fail),flush=True)
    return 0 if results["passed"] else 1


if __name__=="__main__":
    sys.exit(main())
