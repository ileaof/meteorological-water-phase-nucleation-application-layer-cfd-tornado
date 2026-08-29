"""Composite two-level Poisson -- the correct coarse-fine interface stencil (1-D).

This is the *crux* of the AMR projection (`docs/amr_design.md` Milestone 2): a
composite solve where a fine patch replaces part of a coarse grid and the
coarse-fine interface is handled so the result is **2nd-order accurate** and the
interface flux is **single-valued (conservative)**.

The interface flux uses a **2nd-order ghost**: a quadratic through the two nearest
fine cells and the adjacent coarse cell, evaluated at the ghost point that is
symmetric to the fine boundary cell about the interface -- so ``(ghost - fine)/hf``
is the gradient *at* the interface to 2nd order.  The **same** flux value is used
for the coarse cell's face and the fine cells' faces, so the interface is
conservative by construction.  The composite operator is assembled sparsely and
solved directly (a reference solve; the multigrid V-cycle in
:mod:`storm_dynamics.poisson_mg` is the scalable smoother/solver for each level).

Verified (``test_composite_poisson_1d_second_order``): on a manufactured solution
``phi = sin(2 pi x)`` the error falls by ~4x per 2x refinement (error ∝ h²).

The 1-D kernel nails the normal-direction interface stencil; :func:`solve_2d` and
:func:`solve_3d` add the tangential coarse interpolation (linear in 2-D, bilinear
in 3-D) and are verified 2nd-order including the patch edges/corners.
:func:`project_divergence_2d` then wires it into a **two-level MAC projection**: a
face-flux velocity is made discretely divergence-free to machine precision
*including at the coarse-fine interface* (the composite Poisson in its anelastic
role), because the divergence, gradient and Laplacian share the single-valued
interface flux (``L = div . grad``).
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


def _ghost_weights(r: int):
    """Lagrange weights for the 2nd-order interface ghost (nodes at -3/2, -1/2 fine
    and +r/2 coarse, in hf units; evaluated at the ghost point +1/2)."""
    s0, s1, s2, xg = -1.5, -0.5, r / 2.0, 0.5
    a0 = (xg - s1) * (xg - s2) / ((s0 - s1) * (s0 - s2))
    a1 = (xg - s0) * (xg - s2) / ((s1 - s0) * (s1 - s2))
    a2 = (xg - s0) * (xg - s1) / ((s2 - s0) * (s2 - s1))
    return a0, a1, a2


def solve_1d(f_coarse, f_fine, nc, r, ci0, ci1, anchor_value=None):
    """Solve ``lap(phi) = f`` on a composite 1-D grid: periodic coarse grid of ``nc``
    cells with a fine patch (refinement ``r``) over coarse cells ``[ci0, ci1)``.

    ``f_coarse`` (nc,) and ``f_fine`` (r*(ci1-ci0),) are the RHS on each level.
    ``anchor_value`` pins the null space (periodic Poisson is singular): coarse
    cell 0 is set to it (default 0).  Returns ``(phi_coarse, phi_fine)`` with the
    covered coarse cells left as NaN.
    """
    hc = 1.0 / nc
    hf = hc / r
    nf = r * (ci1 - ci0)
    cids = {i: k for k, i in enumerate(i for i in range(nc) if not (ci0 <= i < ci1))}
    fbase = len(cids)
    N = fbase + nf
    fid = lambda k: fbase + k
    gcol = lambda kind, idx: fid(idx) if kind == "f" else cids[idx]

    a0, a1, a2 = _ghost_weights(r)
    Fright = {("f", nf - 2): a0 / hf, ("f", nf - 1): (a1 - 1.0) / hf, ("c", ci1): a2 / hf}
    Fleft = {("f", 1): -a0 / hf, ("f", 0): (1.0 - a1) / hf, ("c", ci0 - 1): -a2 / hf}

    rows, cols, data, rhs = [], [], [], np.zeros(N)
    add = lambda rr, cc, vv: (rows.append(rr), cols.append(cc), data.append(vv))

    for i in cids:                                    # coarse: (F_{i+1/2}-F_{i-1/2})/hc
        row = cids[i]
        if i == ci0 - 1:
            for (k, idx), w in Fleft.items():
                add(row, gcol(k, idx), w / hc)
        else:
            add(row, cids[(i + 1) % nc], 1.0 / hc ** 2); add(row, row, -1.0 / hc ** 2)
        if i == ci1:
            for (k, idx), w in Fright.items():
                add(row, gcol(k, idx), -w / hc)
        else:
            add(row, cids[(i - 1) % nc], 1.0 / hc ** 2); add(row, row, -1.0 / hc ** 2)
        rhs[row] = f_coarse[i]

    for k in range(nf):                               # fine: (F_{k+1/2}-F_{k-1/2})/hf
        row = fid(k)
        if k == nf - 1:
            for (kk, idx), w in Fright.items():
                add(row, gcol(kk, idx), w / hf)
        else:
            add(row, fid(k + 1), 1.0 / hf ** 2); add(row, row, -1.0 / hf ** 2)
        if k == 0:
            for (kk, idx), w in Fleft.items():
                add(row, gcol(kk, idx), -w / hf)
        else:
            add(row, fid(k - 1), 1.0 / hf ** 2); add(row, row, -1.0 / hf ** 2)
        rhs[row] = f_fine[k]

    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N)).tolil()
    A[cids[0], :] = 0; A[cids[0], cids[0]] = 1.0
    rhs[cids[0]] = 0.0 if anchor_value is None else anchor_value
    sol = spla.spsolve(A.tocsr(), rhs)

    phi_c = np.full(nc, np.nan)
    for i in cids:
        phi_c[i] = sol[cids[i]]
    phi_f = np.array([sol[fid(k)] for k in range(nf)])
    return phi_c, phi_f


def manufactured_error(nc, r=2):
    """max error of the composite solve vs phi=sin(2pi x); returns (err, hc)."""
    hc = 1.0 / nc
    ci0, ci1 = nc // 3, 2 * nc // 3
    nf = r * (ci1 - ci0)
    xc = (np.arange(nc) + 0.5) * hc
    xf = ci0 * hc + (np.arange(nf) + 0.5) * (hc / r)
    exact_c = np.sin(2 * np.pi * xc)
    exact_f = np.sin(2 * np.pi * xf)
    fc = -4 * np.pi ** 2 * exact_c
    ff = -4 * np.pi ** 2 * exact_f
    pc, pf = solve_1d(fc, ff, nc, r, ci0, ci1, anchor_value=exact_c[0])
    ec = np.nanmax(np.abs(pc - exact_c))
    ef = np.abs(pf - exact_f).max()
    return float(max(ec, ef)), hc


def solve_2d(f_coarse, f_fine, nc, r, ci0, ci1, cj0, cj1, anchor_value=0.0):
    """Composite 2-D Poisson: periodic coarse grid (nc x nc) with a fine patch
    (refinement ``r``) over coarse cells ``[ci0,ci1) x [cj0,cj1)``.

    Normal direction: the same 2nd-order ghost as :func:`solve_1d`. Tangential:
    linear interpolation of the coarse solution to the fine cell position. The
    interface flux is oriented ``d(phi)/d(+axis)`` and is **single-valued**
    (conservative). Assembled sparsely, solved directly.  Returns
    ``(phi_coarse, phi_fine)`` (covered coarse cells are NaN).
    """
    hc = 1.0 / nc; hf = hc / r
    nfx, nfy = r * (ci1 - ci0), r * (cj1 - cj0)
    covered = lambda I, J: (ci0 <= I < ci1) and (cj0 <= J < cj1)
    cids, g = {}, 0
    for I in range(nc):
        for J in range(nc):
            if not covered(I, J):
                cids[(I, J)] = g; g += 1
    fbase = g
    fid = lambda a, b: fbase + a * nfy + b
    N = fbase + nfx * nfy
    a0, a1, a2 = _ghost_weights(r)
    gcol = lambda k: (fid(k[1], k[2]) if k[0] == "f" else cids[(k[1], k[2])])

    def tang(colfix_is_x, fixed, b):            # linear coarse interp along the edge
        base = cj0 if colfix_is_x else ci0
        Jc = base + b // r
        t = ((b % r) + 0.5) / r - 0.5
        Jn = Jc + (1 if t > 0 else -1)
        w0, w1 = 1 - abs(t), abs(t)
        return ({(fixed, Jc): w0, (fixed, Jn % nc): w1} if colfix_is_x
                else {(Jc, fixed): w0, (Jn % nc, fixed): w1})

    def ghost(edge, b):
        if edge == "R":
            d = {("f", nfx - 2, b): a0, ("f", nfx - 1, b): a1}; cc = tang(True, ci1, b)
        elif edge == "L":
            d = {("f", 1, b): a0, ("f", 0, b): a1}; cc = tang(True, ci0 - 1, b)
        elif edge == "T":
            d = {("f", b, nfy - 2): a0, ("f", b, nfy - 1): a1}; cc = tang(False, cj1, b)
        else:
            d = {("f", b, 1): a0, ("f", b, 0): a1}; cc = tang(False, cj0 - 1, b)
        for (I, J), w in cc.items():
            d[("c", I, J)] = d.get(("c", I, J), 0.0) + a2 * w
        return d

    def flux(edge, b):                          # d(phi)/d(+axis), single-valued
        fb = {"R": (nfx - 1, b), "L": (0, b), "T": (b, nfy - 1), "B": (b, 0)}[edge]
        key = ("f", fb[0], fb[1])
        if edge in ("R", "T"):
            d = {k: v / hf for k, v in ghost(edge, b).items()}
            d[key] = d.get(key, 0.0) - 1.0 / hf
        else:
            d = {k: -v / hf for k, v in ghost(edge, b).items()}
            d[key] = d.get(key, 0.0) + 1.0 / hf
        return d

    rows, cols, data, rhs = [], [], [], np.zeros(N)
    add = lambda rr, cc, vv: (rows.append(rr), cols.append(cc), data.append(vv))

    for a in range(nfx):                        # fine cells
        for b in range(nfy):
            row = fid(a, b)
            if a == nfx - 1:
                for k, w in flux("R", b).items():
                    add(row, gcol(k), w / hf)
            else:
                add(row, fid(a + 1, b), 1 / hf ** 2); add(row, row, -1 / hf ** 2)
            if a == 0:
                for k, w in flux("L", b).items():
                    add(row, gcol(k), -w / hf)
            else:
                add(row, fid(a - 1, b), 1 / hf ** 2); add(row, row, -1 / hf ** 2)
            if b == nfy - 1:
                for k, w in flux("T", a).items():
                    add(row, gcol(k), w / hf)
            else:
                add(row, fid(a, b + 1), 1 / hf ** 2); add(row, row, -1 / hf ** 2)
            if b == 0:
                for k, w in flux("B", a).items():
                    add(row, gcol(k), -w / hf)
            else:
                add(row, fid(a, b - 1), 1 / hf ** 2); add(row, row, -1 / hf ** 2)
            rhs[row] = f_fine[a, b]

    for (I, J), row in cids.items():            # coarse cells
        rhs[row] = f_coarse[I, J]
        for (dI, dJ, edge, along) in ((1, 0, "L", J), (-1, 0, "R", J),
                                      (0, 1, "B", I), (0, -1, "T", I)):
            In, Jn = (I + dI) % nc, (J + dJ) % nc
            if covered(In, Jn):
                base = (along - (cj0 if edge in ("L", "R") else ci0)) * r
                s2 = -(1.0 if edge in ("R", "T") else -1.0)
                for bb in range(base, base + r):
                    for k, w in flux(edge, bb).items():
                        add(row, gcol(k), s2 * (w / r) / hc)
            else:
                add(row, cids[(In, Jn)], 1 / hc ** 2); add(row, row, -1 / hc ** 2)

    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N)).tolil()
    A[cids[(0, 0)], :] = 0; A[cids[(0, 0)], cids[(0, 0)]] = 1.0; rhs[cids[(0, 0)]] = anchor_value
    sol = spla.spsolve(A.tocsr(), rhs)
    phi_c = np.full((nc, nc), np.nan)
    for (I, J) in cids:
        phi_c[I, J] = sol[cids[(I, J)]]
    phi_f = np.array([[sol[fid(a, b)] for b in range(nfy)] for a in range(nfx)])
    return phi_c, phi_f


def manufactured_error_2d(nc, r=2):
    """max error of the 2-D composite solve vs phi=sin(2pi x)sin(2pi y); (err, hc)."""
    hc = 1.0 / nc
    ci0, ci1, cj0, cj1 = nc // 3, 2 * nc // 3, nc // 3, 2 * nc // 3
    nfx, nfy = r * (ci1 - ci0), r * (cj1 - cj0)
    xc = (np.arange(nc) + 0.5) * hc
    Xc, Yc = np.meshgrid(xc, xc, indexing="ij")
    ex_c = np.sin(2 * np.pi * Xc) * np.sin(2 * np.pi * Yc)
    xf = ci0 * hc + (np.arange(nfx) + 0.5) * (hc / r)
    yf = cj0 * hc + (np.arange(nfy) + 0.5) * (hc / r)
    Xf, Yf = np.meshgrid(xf, yf, indexing="ij")
    ex_f = np.sin(2 * np.pi * Xf) * np.sin(2 * np.pi * Yf)
    pc, pf = solve_2d(-8 * np.pi ** 2 * ex_c, -8 * np.pi ** 2 * ex_f,
                      nc, r, ci0, ci1, cj0, cj1, anchor_value=ex_c[0, 0])
    return float(max(np.nanmax(np.abs(pc - ex_c)), np.abs(pf - ex_f).max())), hc


def solve_3d(f_coarse, f_fine, nc, r, lo, hi, anchor_value=0.0):
    """Composite 3-D Poisson: periodic coarse grid (nc^3) with a fine box patch
    (refinement ``r``) over coarse cells ``[lo:hi)`` (``lo=(ci0,cj0,ck0)``,
    ``hi=(ci1,cj1,ck1)``).

    Same construction as :func:`solve_2d` in each of the six faces: a 2nd-order
    quadratic ghost in the normal direction, a **bilinear** coarse interpolation in
    the two tangential directions, and a single-valued interface flux oriented
    ``d(phi)/d(+axis)`` (conservative).  Assembled sparsely, solved directly.
    Returns ``(phi_coarse, phi_fine)`` (covered coarse cells are NaN).
    """
    (ci0, cj0, ck0), (ci1, cj1, ck1) = lo, hi
    hc = 1.0 / nc; hf = hc / r
    nfx, nfy, nfz = r * (ci1 - ci0), r * (cj1 - cj0), r * (ck1 - ck0)
    cov = lambda I, J, K: ci0 <= I < ci1 and cj0 <= J < cj1 and ck0 <= K < ck1
    cids, g = {}, 0
    for I in range(nc):
        for J in range(nc):
            for K in range(nc):
                if not cov(I, J, K):
                    cids[(I, J, K)] = g; g += 1
    fb = g
    fid = lambda a, b, c: fb + (a * nfy + b) * nfz + c
    N = fb + nfx * nfy * nfz
    a0, a1, a2 = _ghost_weights(r)
    gcol = lambda k: (fid(k[1], k[2], k[3]) if k[0] == "f" else cids[(k[1], k[2], k[3])])

    def tangw(t):                          # 1-D linear coarse weights around cell t//r
        Jc = t // r; s = ((t % r) + 0.5) / r - 0.5
        Jn = Jc + (1 if s > 0 else -1)
        return [(Jc, 1 - abs(s)), (Jn, abs(s))]

    def coarse_at(axis, fixedI, u, v, b0, b1):   # bilinear over the 2 tangential axes
        out = {}
        for (Ju, wU) in [(b0 + jj, w) for jj, w in tangw(u)]:
            for (Jv, wV) in [(b1 + kk, w) for kk, w in tangw(v)]:
                key = ((fixedI, Ju % nc, Jv % nc) if axis == 0 else
                       (Ju % nc, fixedI, Jv % nc) if axis == 1 else
                       (Ju % nc, Jv % nc, fixedI))
                out[key] = out.get(key, 0.0) + wU * wV
        return out

    def ghost(face, u, v):
        if face == "Xp":
            d = {("f", nfx - 2, u, v): a0, ("f", nfx - 1, u, v): a1}; cc = coarse_at(0, ci1, u, v, cj0, ck0)
        elif face == "Xm":
            d = {("f", 1, u, v): a0, ("f", 0, u, v): a1}; cc = coarse_at(0, ci0 - 1, u, v, cj0, ck0)
        elif face == "Yp":
            d = {("f", u, nfy - 2, v): a0, ("f", u, nfy - 1, v): a1}; cc = coarse_at(1, cj1, u, v, ci0, ck0)
        elif face == "Ym":
            d = {("f", u, 1, v): a0, ("f", u, 0, v): a1}; cc = coarse_at(1, cj0 - 1, u, v, ci0, ck0)
        elif face == "Zp":
            d = {("f", u, v, nfz - 2): a0, ("f", u, v, nfz - 1): a1}; cc = coarse_at(2, ck1, u, v, ci0, cj0)
        else:
            d = {("f", u, v, 1): a0, ("f", u, v, 0): a1}; cc = coarse_at(2, ck0 - 1, u, v, ci0, cj0)
        for key, w in cc.items():
            d[("c",) + key] = d.get(("c",) + key, 0.0) + a2 * w
        return d

    def flux(face, u, v):                  # d(phi)/d(+axis), single-valued
        bc = {"Xp": (nfx - 1, u, v), "Xm": (0, u, v), "Yp": (u, nfy - 1, v),
              "Ym": (u, 0, v), "Zp": (u, v, nfz - 1), "Zm": (u, v, 0)}[face]
        key = ("f",) + bc
        if face in ("Xp", "Yp", "Zp"):
            d = {k: w / hf for k, w in ghost(face, u, v).items()}; d[key] = d.get(key, 0.0) - 1.0 / hf
        else:
            d = {k: -w / hf for k, w in ghost(face, u, v).items()}; d[key] = d.get(key, 0.0) + 1.0 / hf
        return d

    rows, cols, data, rhs = [], [], [], np.zeros(N)
    add = lambda rr, cc, vv: (rows.append(rr), cols.append(cc), data.append(vv))

    for a in range(nfx):                   # fine cells (6-point + interface fluxes)
        for b in range(nfy):
            for c in range(nfz):
                row = fid(a, b, c)
                for (idx, n, fp, fm, uu, vv) in (
                        (a, nfx, "Xp", "Xm", b, c), (b, nfy, "Yp", "Ym", a, c),
                        (c, nfz, "Zp", "Zm", a, b)):
                    if idx == n - 1:
                        for k, w in flux(fp, uu, vv).items():
                            add(row, gcol(k), w / hf)
                    else:
                        nb = fid(a + (fp == "Xp"), b + (fp == "Yp"), c + (fp == "Zp"))
                        add(row, nb, 1 / hf ** 2); add(row, row, -1 / hf ** 2)
                    if idx == 0:
                        for k, w in flux(fm, uu, vv).items():
                            add(row, gcol(k), -w / hf)
                    else:
                        nb = fid(a - (fm == "Xm"), b - (fm == "Ym"), c - (fm == "Zm"))
                        add(row, nb, 1 / hf ** 2); add(row, row, -1 / hf ** 2)
                rhs[row] = f_fine[a, b, c]

    dirs = [(1, 0, 0, "Xm", 1, 2, cj0, ck0), (-1, 0, 0, "Xp", 1, 2, cj0, ck0),
            (0, 1, 0, "Ym", 0, 2, ci0, ck0), (0, -1, 0, "Yp", 0, 2, ci0, ck0),
            (0, 0, 1, "Zm", 0, 1, ci0, cj0), (0, 0, -1, "Zp", 0, 1, ci0, cj0)]
    for (I, J, K), row in cids.items():    # coarse cells (interface faces refluxed)
        rhs[row] = f_coarse[I, J, K]
        pos = (I, J, K)
        for (dI, dJ, dK, face, ta, tb, b0, b1) in dirs:
            nb = ((I + dI) % nc, (J + dJ) % nc, (K + dK) % nc)
            if cov(*nb):
                s2 = -(1.0 if face in ("Xp", "Yp", "Zp") else -1.0)
                ua, ub = (pos[ta] - b0) * r, (pos[tb] - b1) * r
                for uu in range(ua, ua + r):
                    for vv in range(ub, ub + r):
                        for k, w in flux(face, uu, vv).items():
                            add(row, gcol(k), s2 * (w / r ** 2) / hc)
            else:
                add(row, cids[nb], 1 / hc ** 2); add(row, row, -1 / hc ** 2)

    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N)).tolil()
    anc = cids[(0, 0, 0)]
    A[anc, :] = 0; A[anc, anc] = 1.0; rhs[anc] = anchor_value
    sol = spla.spsolve(A.tocsr(), rhs)
    phi_c = np.full((nc, nc, nc), np.nan)
    for k in cids:
        phi_c[k] = sol[cids[k]]
    phi_f = sol[fb:].reshape(nfx, nfy, nfz)
    return phi_c, phi_f


def manufactured_error_3d(nc, r=2):
    """max error of the 3-D composite solve vs phi=prod sin(2pi x_i); (err, hc)."""
    hc = 1.0 / nc
    lo = (nc // 3, nc // 3, nc // 3); hi = (2 * nc // 3, 2 * nc // 3, 2 * nc // 3)
    (ci0, cj0, ck0), (ci1, cj1, ck1) = lo, hi
    s = lambda x: np.sin(2 * np.pi * x)
    xc = (np.arange(nc) + 0.5) * hc
    ex_c = s(xc)[:, None, None] * s(xc)[None, :, None] * s(xc)[None, None, :]
    fx = ci0 * hc + (np.arange(r * (ci1 - ci0)) + 0.5) * (hc / r)
    fy = cj0 * hc + (np.arange(r * (cj1 - cj0)) + 0.5) * (hc / r)
    fz = ck0 * hc + (np.arange(r * (ck1 - ck0)) + 0.5) * (hc / r)
    ex_f = s(fx)[:, None, None] * s(fy)[None, :, None] * s(fz)[None, None, :]
    pc, pf = solve_3d(-12 * np.pi ** 2 * ex_c, -12 * np.pi ** 2 * ex_f,
                      nc, r, lo, hi, anchor_value=ex_c[0, 0, 0])
    return float(max(np.nanmax(np.abs(pc - ex_c)), np.abs(pf - ex_f).max())), hc


def _interface_flux_2d(nc, r, ci0, ci1, cj0, cj1):
    """The single-valued coarse-fine interface flux operator used by :func:`solve_2d`
    (quadratic normal ghost + linear tangential coarse interp, oriented d(phi)/d+axis).
    Returns ``flux(edge, b) -> {key: weight}`` with ``key`` an ``("f",a,b)`` fine or
    ``("c",I,J)`` coarse cell -- the exact operator the composite Laplacian is built
    from, so a divergence/gradient built on it satisfies ``L = div . grad``."""
    hf = (1.0 / nc) / r
    nfx, nfy = r * (ci1 - ci0), r * (cj1 - cj0)
    a0, a1, a2 = _ghost_weights(r)

    def tang(colfix_is_x, fixed, b):
        base = cj0 if colfix_is_x else ci0
        Jc = base + b // r; t = ((b % r) + 0.5) / r - 0.5
        Jn = Jc + (1 if t > 0 else -1); w0, w1 = 1 - abs(t), abs(t)
        return ({(fixed, Jc): w0, (fixed, Jn % nc): w1} if colfix_is_x
                else {(Jc, fixed): w0, (Jn % nc, fixed): w1})

    def ghost(edge, b):
        if edge == "R": d = {("f", nfx-2, b): a0, ("f", nfx-1, b): a1}; cc = tang(True, ci1, b)
        elif edge == "L": d = {("f", 1, b): a0, ("f", 0, b): a1}; cc = tang(True, ci0-1, b)
        elif edge == "T": d = {("f", b, nfy-2): a0, ("f", b, nfy-1): a1}; cc = tang(False, cj1, b)
        else: d = {("f", b, 1): a0, ("f", b, 0): a1}; cc = tang(False, cj0-1, b)
        for (I, J), w in cc.items(): d[("c", I, J)] = d.get(("c", I, J), 0.0) + a2 * w
        return d

    def flux(edge, b):
        fb = {"R": (nfx-1, b), "L": (0, b), "T": (b, nfy-1), "B": (b, 0)}[edge]
        key = ("f", fb[0], fb[1])
        if edge in ("R", "T"):
            d = {k: v/hf for k, v in ghost(edge, b).items()}; d[key] = d.get(key, 0.0) - 1.0/hf
        else:
            d = {k: -v/hf for k, v in ghost(edge, b).items()}; d[key] = d.get(key, 0.0) + 1.0/hf
        return d

    return flux


def project_divergence_2d(nc, r=2, seed=0):
    """Two-level composite MAC **projection** demonstrator (the composite Poisson in
    its anelastic role).  A random face-flux velocity ``u*`` (fine interior faces,
    single-valued interface faces, coarse faces) is made discretely divergence-free:

        ``f = div(u*)`` -> ``solve_2d`` gives ``p`` (L p = f) -> ``u = u* - grad(p)``

    Because the divergence, the gradient and the composite Laplacian ``L`` all use the
    **same** single-valued interface flux, ``L = div . grad`` and ``div(u) = f - L p``
    vanishes to the solve tolerance -- *including at the coarse-fine interface*.  The
    divergence/gradient here are built independently of :func:`solve_2d`, so a nonzero
    result would expose any inconsistency (self-validating).

    Returns ``(max|div u| coarse, max|div u| fine, max|div u| on fine interface cells)``.
    """
    ci0, ci1, cj0, cj1 = nc // 3, 2 * nc // 3, nc // 3, 2 * nc // 3
    hc = 1.0 / nc; hf = hc / r
    nfx, nfy = r * (ci1 - ci0), r * (cj1 - cj0)
    cov = lambda I, J: ci0 <= I < ci1 and cj0 <= J < cj1
    flux = _interface_flux_2d(nc, r, ci0, ci1, cj0, cj1)

    def evalflux(edge, b, pc, pf):
        return sum(w * (pf[k[1], k[2]] if k[0] == "f" else pc[k[1], k[2]])
                   for k, w in flux(edge, b).items())

    rng = np.random.default_rng(seed)
    ufx = rng.standard_normal((nfx - 1, nfy))        # fine interior x-faces
    ufy = rng.standard_normal((nfx, nfy - 1))        # fine interior y-faces
    ui = {e: rng.standard_normal(nfx) for e in ("R", "L", "T", "B")}   # interface faces
    cfx = rng.standard_normal((nc, nc)); cfy = rng.standard_normal((nc, nc))  # coarse faces

    def cface_x(ui, cfx, I, J):                      # +x coarse face; interface -> mean fine
        if cov((I+1) % nc, J):
            base = (J-cj0)*r; return np.mean([ui["L"][b] for b in range(base, base+r)])
        if cov(I, J):
            base = (J-cj0)*r; return np.mean([ui["R"][b] for b in range(base, base+r)])
        return cfx[I, J]

    def cface_y(ui, cfy, I, J):
        if cov(I, (J+1) % nc):
            base = (I-ci0)*r; return np.mean([ui["B"][a] for a in range(base, base+r)])
        if cov(I, J):
            base = (I-ci0)*r; return np.mean([ui["T"][a] for a in range(base, base+r)])
        return cfy[I, J]

    def divergence(ufx, ufy, ui, cfx, cfy):
        dc = np.zeros((nc, nc)); df = np.zeros((nfx, nfy))
        for a in range(nfx):
            for b in range(nfy):
                xR = ui["R"][b] if a == nfx-1 else ufx[a, b]
                xL = ui["L"][b] if a == 0 else ufx[a-1, b]
                yT = ui["T"][a] if b == nfy-1 else ufy[a, b]
                yB = ui["B"][a] if b == 0 else ufy[a, b-1]
                df[a, b] = (xR - xL)/hf + (yT - yB)/hf
        for I in range(nc):
            for J in range(nc):
                if cov(I, J): continue
                dc[I, J] = ((cface_x(ui, cfx, I, J) - cface_x(ui, cfx, (I-1) % nc, J))/hc
                            + (cface_y(ui, cfy, I, J) - cface_y(ui, cfy, I, (J-1) % nc))/hc)
        return dc, df

    fc, ff = divergence(ufx, ufy, ui, cfx, cfy)
    pc, pf = solve_2d(fc, ff, nc, r, ci0, ci1, cj0, cj1, anchor_value=0.0)
    pc = np.nan_to_num(pc)

    gfx = np.zeros_like(ufx); gfy = np.zeros_like(ufy)
    gfx[:] = (pf[1:, :] - pf[:-1, :]) / hf
    gfy[:] = (pf[:, 1:] - pf[:, :-1]) / hf
    gi = {e: np.array([evalflux(e, b, pc, pf) for b in range(nfx)]) for e in ("R", "L", "T", "B")}
    gcfx = np.zeros((nc, nc)); gcfy = np.zeros((nc, nc))
    for I in range(nc):
        for J in range(nc):
            if not (cov(I, J) or cov((I+1) % nc, J)): gcfx[I, J] = (pc[(I+1) % nc, J]-pc[I, J])/hc
            if not (cov(I, J) or cov(I, (J+1) % nc)): gcfy[I, J] = (pc[I, (J+1) % nc]-pc[I, J])/hc

    ufx2, ufy2 = ufx - gfx, ufy - gfy
    ui2 = {e: ui[e] - gi[e] for e in ui}
    cfx2, cfy2 = cfx - gcfx, cfy - gcfy
    dc2, df2 = divergence(ufx2, ufy2, ui2, cfx2, cfy2)
    dc2[0, 0] = 0.0                                  # exclude the anchor cell
    border = np.zeros((nfx, nfy), bool)
    border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    return (float(np.abs(dc2).max()), float(np.abs(df2).max()), float(np.abs(df2[border]).max()))


def _interface_flux_3d(nc, r, lo, hi):
    """The single-valued interface flux operator used by :func:`solve_3d` (quadratic
    normal ghost + bilinear tangential coarse interp, oriented d(phi)/d+axis).
    Returns ``flux(face, u, v) -> {key: weight}`` (face in Xp/Xm/Yp/Ym/Zp/Zm)."""
    (ci0, cj0, ck0), (ci1, cj1, ck1) = lo, hi
    hf = (1.0 / nc) / r
    nfx, nfy, nfz = r * (ci1 - ci0), r * (cj1 - cj0), r * (ck1 - ck0)
    a0, a1, a2 = _ghost_weights(r)

    def tangw(t):
        Jc = t // r; s = ((t % r) + 0.5) / r - 0.5
        Jn = Jc + (1 if s > 0 else -1)
        return [(Jc, 1 - abs(s)), (Jn, abs(s))]

    def coarse_at(axis, fixedI, u, v, b0, b1):
        out = {}
        for (Ju, wU) in [(b0 + jj, w) for jj, w in tangw(u)]:
            for (Jv, wV) in [(b1 + kk, w) for kk, w in tangw(v)]:
                key = ((fixedI, Ju % nc, Jv % nc) if axis == 0 else
                       (Ju % nc, fixedI, Jv % nc) if axis == 1 else
                       (Ju % nc, Jv % nc, fixedI))
                out[key] = out.get(key, 0.0) + wU * wV
        return out

    def ghost(face, u, v):
        if face == "Xp": d = {("f", nfx-2, u, v): a0, ("f", nfx-1, u, v): a1}; cc = coarse_at(0, ci1, u, v, cj0, ck0)
        elif face == "Xm": d = {("f", 1, u, v): a0, ("f", 0, u, v): a1}; cc = coarse_at(0, ci0-1, u, v, cj0, ck0)
        elif face == "Yp": d = {("f", u, nfy-2, v): a0, ("f", u, nfy-1, v): a1}; cc = coarse_at(1, cj1, u, v, ci0, ck0)
        elif face == "Ym": d = {("f", u, 1, v): a0, ("f", u, 0, v): a1}; cc = coarse_at(1, cj0-1, u, v, ci0, ck0)
        elif face == "Zp": d = {("f", u, v, nfz-2): a0, ("f", u, v, nfz-1): a1}; cc = coarse_at(2, ck1, u, v, ci0, cj0)
        else: d = {("f", u, v, 1): a0, ("f", u, v, 0): a1}; cc = coarse_at(2, ck0-1, u, v, ci0, cj0)
        for key, w in cc.items(): d[("c",) + key] = d.get(("c",) + key, 0.0) + a2 * w
        return d

    def flux(face, u, v):
        bc = {"Xp": (nfx-1, u, v), "Xm": (0, u, v), "Yp": (u, nfy-1, v),
              "Ym": (u, 0, v), "Zp": (u, v, nfz-1), "Zm": (u, v, 0)}[face]
        key = ("f",) + bc
        if face in ("Xp", "Yp", "Zp"):
            d = {k: w/hf for k, w in ghost(face, u, v).items()}; d[key] = d.get(key, 0.0) - 1.0/hf
        else:
            d = {k: -w/hf for k, w in ghost(face, u, v).items()}; d[key] = d.get(key, 0.0) + 1.0/hf
        return d

    return flux


def project_divergence_3d(nc, r=2, seed=0):
    """Two-level composite MAC projection in **3-D** (the mechanical analogue of
    :func:`project_divergence_2d`, built on :func:`solve_3d`).  A random face-flux
    velocity is made discretely divergence-free by ``div(u*)`` → ``solve_3d`` →
    ``u = u* - grad(p)``; ``div(u)`` vanishes to the solve tolerance *including at the
    coarse-fine interface*.  Returns ``(max|div u| coarse, fine, fine-interface)``."""
    lo = (nc // 3, nc // 3, nc // 3); hi = (2 * nc // 3, 2 * nc // 3, 2 * nc // 3)
    (ci0, cj0, ck0), (ci1, cj1, ck1) = lo, hi
    hc = 1.0 / nc; hf = hc / r
    nfx, nfy, nfz = r * (ci1 - ci0), r * (cj1 - cj0), r * (ck1 - ck0)
    cov = lambda I, J, K: ci0 <= I < ci1 and cj0 <= J < cj1 and ck0 <= K < ck1
    flux = _interface_flux_3d(nc, r, lo, hi)

    def evalflux(face, u, v, pc, pf):
        return sum(w * (pf[k[1], k[2], k[3]] if k[0] == "f" else pc[k[1], k[2], k[3]])
                   for k, w in flux(face, u, v).items())

    rng = np.random.default_rng(seed)
    ufx = rng.standard_normal((nfx - 1, nfy, nfz))
    ufy = rng.standard_normal((nfx, nfy - 1, nfz))
    ufz = rng.standard_normal((nfx, nfy, nfz - 1))
    ui = {"Xp": rng.standard_normal((nfy, nfz)), "Xm": rng.standard_normal((nfy, nfz)),
          "Yp": rng.standard_normal((nfx, nfz)), "Ym": rng.standard_normal((nfx, nfz)),
          "Zp": rng.standard_normal((nfx, nfy)), "Zm": rng.standard_normal((nfx, nfy))}
    cfx = rng.standard_normal((nc, nc, nc)); cfy = rng.standard_normal((nc, nc, nc))
    cfz = rng.standard_normal((nc, nc, nc))

    def mean_iface(ui, side, p0, p1):  # mean of the r*r fine interface faces at a coarse face
        b0, b1 = p0 * r, p1 * r
        return ui[side][b0:b0 + r, b1:b1 + r].mean()

    def cface_x(ui, cfx, I, J, K):
        if cov((I + 1) % nc, J, K): return mean_iface(ui, "Xm", J - cj0, K - ck0)
        if cov(I, J, K): return mean_iface(ui, "Xp", J - cj0, K - ck0)
        return cfx[I, J, K]

    def cface_y(ui, cfy, I, J, K):
        if cov(I, (J + 1) % nc, K): return mean_iface(ui, "Ym", I - ci0, K - ck0)
        if cov(I, J, K): return mean_iface(ui, "Yp", I - ci0, K - ck0)
        return cfy[I, J, K]

    def cface_z(ui, cfz, I, J, K):
        if cov(I, J, (K + 1) % nc): return mean_iface(ui, "Zm", I - ci0, J - cj0)
        if cov(I, J, K): return mean_iface(ui, "Zp", I - ci0, J - cj0)
        return cfz[I, J, K]

    def divergence(ufx, ufy, ufz, ui, cfx, cfy, cfz):
        dc = np.zeros((nc, nc, nc)); df = np.zeros((nfx, nfy, nfz))
        for a in range(nfx):
            for b in range(nfy):
                for c in range(nfz):
                    xR = ui["Xp"][b, c] if a == nfx-1 else ufx[a, b, c]
                    xL = ui["Xm"][b, c] if a == 0 else ufx[a-1, b, c]
                    yT = ui["Yp"][a, c] if b == nfy-1 else ufy[a, b, c]
                    yB = ui["Ym"][a, c] if b == 0 else ufy[a, b-1, c]
                    zF = ui["Zp"][a, b] if c == nfz-1 else ufz[a, b, c]
                    zN = ui["Zm"][a, b] if c == 0 else ufz[a, b, c-1]
                    df[a, b, c] = ((xR-xL) + (yT-yB) + (zF-zN)) / hf
        for I in range(nc):
            for J in range(nc):
                for K in range(nc):
                    if cov(I, J, K): continue
                    dc[I, J, K] = ((cface_x(ui, cfx, I, J, K) - cface_x(ui, cfx, (I-1) % nc, J, K))
                                   + (cface_y(ui, cfy, I, J, K) - cface_y(ui, cfy, I, (J-1) % nc, K))
                                   + (cface_z(ui, cfz, I, J, K) - cface_z(ui, cfz, I, J, (K-1) % nc))) / hc
        return dc, df

    fc, ff = divergence(ufx, ufy, ufz, ui, cfx, cfy, cfz)
    pc, pf = solve_3d(fc, ff, nc, r, lo, hi, anchor_value=0.0)
    pc = np.nan_to_num(pc)

    gfx = (pf[1:, :, :] - pf[:-1, :, :]) / hf
    gfy = (pf[:, 1:, :] - pf[:, :-1, :]) / hf
    gfz = (pf[:, :, 1:] - pf[:, :, :-1]) / hf
    gi = {}
    for side, (n0, n1) in (("Xp", (nfy, nfz)), ("Xm", (nfy, nfz)), ("Yp", (nfx, nfz)),
                           ("Ym", (nfx, nfz)), ("Zp", (nfx, nfy)), ("Zm", (nfx, nfy))):
        gi[side] = np.array([[evalflux(side, u, v, pc, pf) for v in range(n1)]
                             for u in range(n0)])
    gcfx = np.zeros((nc, nc, nc)); gcfy = np.zeros((nc, nc, nc)); gcfz = np.zeros((nc, nc, nc))
    for I in range(nc):
        for J in range(nc):
            for K in range(nc):
                if not (cov(I, J, K) or cov((I+1) % nc, J, K)): gcfx[I, J, K] = (pc[(I+1) % nc, J, K]-pc[I, J, K])/hc
                if not (cov(I, J, K) or cov(I, (J+1) % nc, K)): gcfy[I, J, K] = (pc[I, (J+1) % nc, K]-pc[I, J, K])/hc
                if not (cov(I, J, K) or cov(I, J, (K+1) % nc)): gcfz[I, J, K] = (pc[I, J, (K+1) % nc]-pc[I, J, K])/hc

    u2 = (ufx - gfx, ufy - gfy, ufz - gfz)
    ui2 = {e: ui[e] - gi[e] for e in ui}
    dc2, df2 = divergence(*u2, ui2, cfx - gcfx, cfy - gcfy, cfz - gcfz)
    dc2[0, 0, 0] = 0.0                               # exclude the anchor cell
    border = np.zeros((nfx, nfy, nfz), bool)
    border[0, :, :] = border[-1, :, :] = True
    border[:, 0, :] = border[:, -1, :] = True
    border[:, :, 0] = border[:, :, -1] = True
    return (float(np.abs(dc2).max()), float(np.abs(df2).max()), float(np.abs(df2[border]).max()))


__all__ = ["solve_1d", "manufactured_error", "solve_2d", "manufactured_error_2d",
           "solve_3d", "manufactured_error_3d", "project_divergence_2d",
           "project_divergence_3d"]


if __name__ == "__main__":
    for dim, fn, grids in (("1-D", manufactured_error, (48, 96, 192, 384)),
                           ("2-D", manufactured_error_2d, (24, 48, 96)),
                           ("3-D", manufactured_error_3d, (10, 20))):
        print("composite 2-level %s Poisson (2nd-order interface stencil):" % dim)
        prev = None
        for nc in grids:
            err, hc = fn(nc)
            ratio = "" if prev is None else "  ratio=%.2f" % (prev / err)
            print("  nc=%3d  hc=%.4f  max|err|=%.3e%s" % (nc, hc, err, ratio))
            prev = err

    for dim, fn, grids in (("2-D", project_divergence_2d, (12, 24)),
                           ("3-D", project_divergence_3d, (8, 12))):
        print("two-level composite %s MAC projection (div u -> 0 across the interface):" % dim)
        for nc in grids:
            dc, df, di = fn(nc, 2)
            print("  nc=%3d  max|div u|: coarse=%.2e  fine=%.2e  fine-interface=%.2e"
                  % (nc, dc, df, di))
