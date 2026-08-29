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
:func:`project_divergence_2d` / :func:`project_divergence_3d` then wire it into a
**two-level MAC projection**: a face-flux velocity is made discretely
divergence-free to machine precision *including at the coarse-fine interface*,
because the divergence, gradient and Laplacian share the single-valued interface
flux (``L = div . grad``).

**Anelastic interpretation.** The storm's projection is ``u = u* - grad(p)/rho0``
with ``div(rho0 u) = 0``.  In mass-flux variables ``m = rho0 u`` this is exactly
``m = m* - grad(p)``, ``div(m) = 0`` -- i.e. the projection here, with the face
field read as the anelastic mass flux ``rho0 u`` (the density weight cancels in the
divergence constraint and only re-enters when recovering ``u = m/rho0``).  So these
functions *are* the anelastic two-level projection across a refinement interface;
no new algorithm is needed to make the storm anelastic.

**Final assembly (the storm's real operator).** :func:`solve_composite_hz` assembles the
unified operator for the storm's nest geometry -- the horizontal composite interface at
every z-level (x,y refined by ``r``) plus a variable-dz finite-volume vertical coupling
(walls top/bottom, matched z, not refined) -- verified 2nd-order.
:func:`composite_project_massflux_hz` is the full 3-D projection on the storm's staggered
C-grid mass fluxes (parent + nest, ``u:(nc+1,nc,nz)`` etc.): it makes ``div(m)=0`` across
the interface to ~1e-13, combining all three plumbing pieces (wall BC, the face bridge,
the stretched-z metric).  The remaining work is only the call site: form ``rho0 u*`` on
both levels, call this, and recover ``u = m/rho0`` -- see ``docs/amr_design.md``.
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


def solve_2d(f_coarse, f_fine, nc, r, ci0, ci1, cj0, cj1, anchor_value=0.0, periodic=True):
    """Composite 2-D Poisson: coarse grid (nc x nc) with a fine patch (refinement
    ``r``) over coarse cells ``[ci0,ci1) x [cj0,cj1)``.

    Normal direction: the same 2nd-order ghost as :func:`solve_1d`. Tangential:
    linear interpolation of the coarse solution to the fine cell position. The
    interface flux is oriented ``d(phi)/d(+axis)`` and is **single-valued**
    (conservative). Assembled sparsely, solved directly.  Returns
    ``(phi_coarse, phi_fine)`` (covered coarse cells are NaN).

    ``periodic`` (default True) wraps the outer coarse boundary; ``False`` applies a
    **solid-wall Neumann** BC there (``dphi/dn = 0`` -- the storm's wall condition),
    i.e. a boundary coarse cell simply drops its outward face.  The patch is interior
    so the interface stencil is unchanged either way.
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
            ni, nj = I + dI, J + dJ
            if not periodic and (ni < 0 or ni >= nc or nj < 0 or nj >= nc):
                continue                        # solid-wall Neumann: drop the face
            In, Jn = ni % nc, nj % nc
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


def manufactured_error_2d_wall(nc, r=2):
    """max error of the **solid-wall** (Neumann) composite solve vs
    phi=cos(pi x)cos(pi y) (which has dphi/dn=0 on the domain walls); (err, hc)."""
    hc = 1.0 / nc
    ci0, ci1, cj0, cj1 = nc // 3, 2 * nc // 3, nc // 3, 2 * nc // 3
    nfx, nfy = r * (ci1 - ci0), r * (cj1 - cj0)
    xc = (np.arange(nc) + 0.5) * hc
    Xc, Yc = np.meshgrid(xc, xc, indexing="ij")
    ex_c = np.cos(np.pi * Xc) * np.cos(np.pi * Yc)
    xf = ci0 * hc + (np.arange(nfx) + 0.5) * (hc / r)
    yf = cj0 * hc + (np.arange(nfy) + 0.5) * (hc / r)
    Xf, Yf = np.meshgrid(xf, yf, indexing="ij")
    ex_f = np.cos(np.pi * Xf) * np.cos(np.pi * Yf)
    pc, pf = solve_2d(-2 * np.pi ** 2 * ex_c, -2 * np.pi ** 2 * ex_f,
                      nc, r, ci0, ci1, cj0, cj1, anchor_value=ex_c[0, 0], periodic=False)
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


def _interface_flux_2d(nc, r, ci0, ci1, cj0, cj1, hx=None):
    """The single-valued coarse-fine interface flux operator used by :func:`solve_2d`
    (quadratic normal ghost + linear tangential coarse interp, oriented d(phi)/d+axis).
    Returns ``flux(edge, b) -> {key: weight}`` with ``key`` an ``("f",a,b)`` fine or
    ``("c",I,J)`` coarse cell -- the exact operator the composite Laplacian is built
    from, so a divergence/gradient built on it satisfies ``L = div . grad``.
    ``hx`` is the physical coarse horizontal spacing (default the unit-domain ``1/nc``)."""
    hf = (hx if hx is not None else 1.0 / nc) / r
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


def project_divergence_2d(nc, r=2, seed=0, periodic=True):
    """Two-level composite MAC **projection** demonstrator (the composite Poisson in
    its anelastic role).  A random face-flux velocity ``u*`` (fine interior faces,
    single-valued interface faces, coarse faces) is made discretely divergence-free:

        ``f = div(u*)`` -> ``solve_2d`` gives ``p`` (L p = f) -> ``u = u* - grad(p)``

    Because the divergence, the gradient and the composite Laplacian ``L`` all use the
    **same** single-valued interface flux, ``L = div . grad`` and ``div(u) = f - L p``
    vanishes to the solve tolerance -- *including at the coarse-fine interface*.  The
    divergence/gradient here are built independently of :func:`solve_2d`, so a nonzero
    result would expose any inconsistency (self-validating).

    ``periodic`` (default True) wraps the outer boundary; ``False`` applies the
    **solid-wall** BC there (zero normal velocity: the outward boundary face carries no
    flux and gets no pressure correction -- the storm/nest wall condition).

    Returns ``(max|div u| coarse, max|div u| fine, max|div u| on fine interface cells)``.
    """
    ci0, ci1, cj0, cj1 = nc // 3, 2 * nc // 3, nc // 3, 2 * nc // 3
    hc = 1.0 / nc; hf = hc / r
    nfx, nfy = r * (ci1 - ci0), r * (cj1 - cj0)
    cov = lambda I, J: ci0 <= I < ci1 and cj0 <= J < cj1
    flux = _interface_flux_2d(nc, r, ci0, ci1, cj0, cj1)
    wxp = lambda I: (not periodic) and I == nc - 1     # +x face of cell I is a wall
    wxm = lambda I: (not periodic) and I == 0
    wyp = lambda J: (not periodic) and J == nc - 1
    wym = lambda J: (not periodic) and J == 0

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
                fxp = 0.0 if wxp(I) else cface_x(ui, cfx, I, J)
                fxm = 0.0 if wxm(I) else cface_x(ui, cfx, (I-1) % nc, J)
                fyp = 0.0 if wyp(J) else cface_y(ui, cfy, I, J)
                fym = 0.0 if wym(J) else cface_y(ui, cfy, I, (J-1) % nc)
                dc[I, J] = ((fxp - fxm) + (fyp - fym)) / hc
        return dc, df

    fc, ff = divergence(ufx, ufy, ui, cfx, cfy)
    pc, pf = solve_2d(fc, ff, nc, r, ci0, ci1, cj0, cj1, anchor_value=0.0, periodic=periodic)
    pc = np.nan_to_num(pc)

    gfx = np.zeros_like(ufx); gfy = np.zeros_like(ufy)
    gfx[:] = (pf[1:, :] - pf[:-1, :]) / hf
    gfy[:] = (pf[:, 1:] - pf[:, :-1]) / hf
    gi = {e: np.array([evalflux(e, b, pc, pf) for b in range(nfx)]) for e in ("R", "L", "T", "B")}
    gcfx = np.zeros((nc, nc)); gcfy = np.zeros((nc, nc))
    for I in range(nc):
        for J in range(nc):
            if not (cov(I, J) or cov((I+1) % nc, J)) and not wxp(I): gcfx[I, J] = (pc[(I+1) % nc, J]-pc[I, J])/hc
            if not (cov(I, J) or cov(I, (J+1) % nc)) and not wyp(J): gcfy[I, J] = (pc[I, (J+1) % nc]-pc[I, J])/hc

    ufx2, ufy2 = ufx - gfx, ufy - gfy
    ui2 = {e: ui[e] - gi[e] for e in ui}
    cfx2, cfy2 = cfx - gcfx, cfy - gcfy
    dc2, df2 = divergence(ufx2, ufy2, ui2, cfx2, cfy2)
    dc2[0, 0] = 0.0                                  # exclude the anchor cell
    border = np.zeros((nfx, nfy), bool)
    border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    return (float(np.abs(dc2).max()), float(np.abs(df2).max()), float(np.abs(df2[border]).max()))


def composite_project_massflux_2d(mu_c, mv_c, mu_f, mv_f, nc, r,
                                   ci0, ci1, cj0, cj1, periodic=False):
    """**Item (b) of the step-2 integration bridge**: project the storm's staggered
    C-grid **mass fluxes** ``m = rho0 u`` on a two-level (coarse parent + fine nest)
    grid so that ``div(m) = 0`` across the coarse-fine interface, operating directly
    on the storm's face-array convention.

    Arrays (modified in place):
      ``mu_c`` (nc+1, nc), ``mv_c`` (nc, nc+1)   -- coarse parent x-/y-face mass fluxes;
      ``mu_f`` (nfx+1, nfy), ``mv_f`` (nfx, nfy+1) -- fine nest face mass fluxes,
      with ``nfx=r*(ci1-ci0)``, ``nfy=r*(cj1-cj0)`` over coarse cells ``[ci0,ci1)x[cj0,cj1)``.

    The bridge maps these to the composite decomposition, runs the composite projection
    (``div`` -> :func:`solve_2d` -> correct with the single-valued interface flux),
    writes the corrected fluxes back, and **refluxes** the parent's interface faces to
    the single-valued fine-face mean.  In mass-flux variables this *is* the anelastic
    projection (the caller recovers ``u = m/rho0``); the density weight is the caller's.

    ``periodic`` wraps the outer boundary; ``False`` (default, the storm/nest case)
    applies the solid-wall BC.  Returns ``(max|div m| coarse, fine, fine-interface)``
    recomputed straight from the written-back arrays (self-validating).
    """
    hc = 1.0 / nc; hf = hc / r
    nfx, nfy = r * (ci1 - ci0), r * (cj1 - cj0)
    cov = lambda I, J: ci0 <= I < ci1 and cj0 <= J < cj1
    flux = _interface_flux_2d(nc, r, ci0, ci1, cj0, cj1)
    wxp = lambda I: (not periodic) and I == nc - 1
    wxm = lambda I: (not periodic) and I == 0
    wyp = lambda J: (not periodic) and J == nc - 1
    wym = lambda J: (not periodic) and J == 0
    if not periodic:
        mu_c[0, :] = mu_c[nc, :] = 0.0; mv_c[:, 0] = mv_c[:, nc] = 0.0

    ufx = mu_f[1:nfx, :].copy(); ufy = mv_f[:, 1:nfy].copy()
    ui = {"L": mu_f[0, :].copy(), "R": mu_f[nfx, :].copy(),
          "B": mv_f[:, 0].copy(), "T": mv_f[:, nfy].copy()}
    cfx = np.array([[mu_c[I + 1, J] for J in range(nc)] for I in range(nc)])
    cfy = np.array([[mv_c[I, J + 1] for J in range(nc)] for I in range(nc)])

    def evalflux(edge, b, pc, pf):
        return sum(w * (pf[k[1], k[2]] if k[0] == "f" else pc[k[1], k[2]])
                   for k, w in flux(edge, b).items())

    def cface_x(ui, cfx, I, J):
        if cov((I + 1) % nc, J): base = (J - cj0) * r; return np.mean(ui["L"][base:base + r])
        if cov(I, J): base = (J - cj0) * r; return np.mean(ui["R"][base:base + r])
        return cfx[I, J]

    def cface_y(ui, cfy, I, J):
        if cov(I, (J + 1) % nc): base = (I - ci0) * r; return np.mean(ui["B"][base:base + r])
        if cov(I, J): base = (I - ci0) * r; return np.mean(ui["T"][base:base + r])
        return cfy[I, J]

    def divergence(ufx, ufy, ui, cfx, cfy):
        dc = np.full((nc, nc), np.nan); df = np.zeros((nfx, nfy))
        for a in range(nfx):
            for b in range(nfy):
                xR = ui["R"][b] if a == nfx - 1 else ufx[a, b]
                xL = ui["L"][b] if a == 0 else ufx[a - 1, b]
                yT = ui["T"][a] if b == nfy - 1 else ufy[a, b]
                yB = ui["B"][a] if b == 0 else ufy[a, b - 1]
                df[a, b] = (xR - xL) / hf + (yT - yB) / hf
        for I in range(nc):
            for J in range(nc):
                if cov(I, J): continue
                fxp = 0.0 if wxp(I) else cface_x(ui, cfx, I, J)
                fxm = 0.0 if wxm(I) else cface_x(ui, cfx, (I - 1) % nc, J)
                fyp = 0.0 if wyp(J) else cface_y(ui, cfy, I, J)
                fym = 0.0 if wym(J) else cface_y(ui, cfy, I, (J - 1) % nc)
                dc[I, J] = ((fxp - fxm) + (fyp - fym)) / hc
        return dc, df

    fc, ff = divergence(ufx, ufy, ui, cfx, cfy)
    fc = np.nan_to_num(fc)
    pc, pf = solve_2d(fc, ff, nc, r, ci0, ci1, cj0, cj1, anchor_value=0.0, periodic=periodic)
    pc = np.nan_to_num(pc)

    gfx = (pf[1:, :] - pf[:-1, :]) / hf
    gfy = (pf[:, 1:] - pf[:, :-1]) / hf
    gi = {e: np.array([evalflux(e, b, pc, pf) for b in range(nfx)]) for e in ("R", "L", "T", "B")}
    ufx -= gfx; ufy -= gfy
    for e in ui: ui[e] -= gi[e]
    for I in range(nc):
        for J in range(nc):
            if not (cov(I, J) or cov((I + 1) % nc, J)) and not wxp(I):
                cfx[I, J] -= (pc[(I + 1) % nc, J] - pc[I, J]) / hc
            if not (cov(I, J) or cov(I, (J + 1) % nc)) and not wyp(J):
                cfy[I, J] -= (pc[I, (J + 1) % nc] - pc[I, J]) / hc

    # write back into the storm arrays
    mu_f[1:nfx, :] = ufx; mv_f[:, 1:nfy] = ufy
    mu_f[0, :] = ui["L"]; mu_f[nfx, :] = ui["R"]; mv_f[:, 0] = ui["B"]; mv_f[:, nfy] = ui["T"]
    for I in range(nc):
        for J in range(nc):
            if not cov(I, J) and not cov((I + 1) % nc, J): mu_c[I + 1, J] = cfx[I, J]
            if not cov(I, J) and not cov(I, (J + 1) % nc): mv_c[I, J + 1] = cfy[I, J]
    for J in range(cj0, cj1):                     # reflux: parent interface faces = fine mean
        base = (J - cj0) * r
        mu_c[ci0, J] = np.mean(ui["L"][base:base + r])
        mu_c[ci1, J] = np.mean(ui["R"][base:base + r])
    for I in range(ci0, ci1):
        base = (I - ci0) * r
        mv_c[I, cj0] = np.mean(ui["B"][base:base + r])
        mv_c[I, cj1] = np.mean(ui["T"][base:base + r])
    if periodic:
        mu_c[0, :] = mu_c[nc, :]; mv_c[:, 0] = mv_c[:, nc]

    # self-validating: div(m) recomputed straight from the written-back arrays
    df_chk = (mu_f[1:, :] - mu_f[:-1, :]) / hf + (mv_f[:, 1:] - mv_f[:, :-1]) / hf
    dc_chk = np.full((nc, nc), np.nan)
    for I in range(nc):
        for J in range(nc):
            if cov(I, J): continue
            fxp = 0.0 if wxp(I) else mu_c[I + 1, J]
            fxm = 0.0 if wxm(I) else mu_c[I, J]
            fyp = 0.0 if wyp(J) else mv_c[I, J + 1]
            fym = 0.0 if wym(J) else mv_c[I, J]
            dc_chk[I, J] = ((fxp - fxm) + (fyp - fym)) / hc
    dc_chk[0, 0] = 0.0
    border = np.zeros((nfx, nfy), bool)
    border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    return (float(np.nanmax(np.abs(dc_chk))), float(np.abs(df_chk).max()),
            float(np.abs(df_chk[border]).max()))


def manufactured_error_metric_z(nc, nz, r=2, s=1.05, Lz=2.0):
    """**Item (c) of the step-2 bridge**: the **stretched vertical metric** composed with
    the horizontal composite interface.  Solves an (x-z) composite matching the storm's
    nest (x horizontally refined by ``r`` -> the composite interface, uniform ``hx``; z
    stretched by ratio ``s`` with a finite-volume variable-dz operator, walls top/bottom,
    NOT refined -- matched levels) on the manufactured ``phi = cos(2 pi x) cos(pi z/Lz)``.
    Returns max|error|.

    Composes correctly: for uniform z (``s=1``) the error is clean 2nd order (~4x per
    2x refinement); moderate stretching is supraconvergent (~1.8-2 order), the standard
    cell-centred finite-volume behaviour on a non-uniform grid (the storm's own solver
    has it too).  The vertical metric is orthogonal to the interface stencil.
    """
    hc = 1.0 / nc; hf = hc / r
    ci0, ci1 = nc // 3, 2 * nc // 3
    nfx = r * (ci1 - ci0)
    dz = s ** np.arange(nz); dz *= Lz / dz.sum()
    zf = np.concatenate([[0.0], np.cumsum(dz)])
    zc = 0.5 * (zf[:-1] + zf[1:]); dzc = zf[1:] - zf[:-1]; dzf = np.diff(zc)
    incov = lambda I: ci0 <= I < ci1
    a0, a1, a2 = _ghost_weights(r)
    Fright = {("f", nfx - 2): a0 / hf, ("f", nfx - 1): (a1 - 1.0) / hf, ("c", ci1): a2 / hf}
    Fleft = {("f", 1): -a0 / hf, ("f", 0): (1.0 - a1) / hf, ("c", ci0 - 1): -a2 / hf}
    cxs = [I for I in range(nc) if not incov(I)]
    cid = {(I, k): idx for idx, (I, k) in enumerate((I, k) for I in cxs for k in range(nz))}
    fb = len(cid)
    fid = lambda a, k: fb + a * nz + k
    N = fb + nfx * nz
    gcol = lambda kind, i, k: (fid(i, k) if kind == "f" else cid[(i, k)])
    rows, cols, data = [], [], []
    add = lambda rr, cc, vv: (rows.append(rr), cols.append(cc), data.append(vv))

    def vertical(row, colof, k):                       # variable-dz FV, walls at 0, nz-1
        if k < nz - 1:
            w = 1.0 / (dzf[k] * dzc[k]); add(row, colof(k + 1), w); add(row, row, -w)
        if k > 0:
            w = 1.0 / (dzf[k - 1] * dzc[k]); add(row, colof(k - 1), w); add(row, row, -w)

    for I in cxs:
        for k in range(nz):
            row = cid[(I, k)]
            if I == ci0 - 1:
                for (kk, idx), wv in Fleft.items(): add(row, gcol(kk, idx, k), wv / hc)
            else:
                add(row, cid[((I + 1) % nc, k)], 1.0 / hc ** 2); add(row, row, -1.0 / hc ** 2)
            if I == ci1:
                for (kk, idx), wv in Fright.items(): add(row, gcol(kk, idx, k), -wv / hc)
            else:
                add(row, cid[((I - 1) % nc, k)], 1.0 / hc ** 2); add(row, row, -1.0 / hc ** 2)
            vertical(row, lambda kk, I=I: cid[(I, kk)], k)
    for a in range(nfx):
        for k in range(nz):
            row = fid(a, k)
            if a == nfx - 1:
                for (kk, idx), wv in Fright.items(): add(row, gcol(kk, idx, k), wv / hf)
            else:
                add(row, fid(a + 1, k), 1.0 / hf ** 2); add(row, row, -1.0 / hf ** 2)
            if a == 0:
                for (kk, idx), wv in Fleft.items(): add(row, gcol(kk, idx, k), -wv / hf)
            else:
                add(row, fid(a - 1, k), 1.0 / hf ** 2); add(row, row, -1.0 / hf ** 2)
            vertical(row, lambda kk, a=a: fid(a, kk), k)

    kz = np.pi / Lz
    xf = ci0 * hc + (np.arange(nfx) + 0.5) * hf
    exc = lambda I, k: np.cos(2 * np.pi * (I + 0.5) * hc) * np.cos(kz * zc[k])
    exf = lambda a, k: np.cos(2 * np.pi * xf[a]) * np.cos(kz * zc[k])
    lam = -(4 * np.pi ** 2 + kz ** 2)
    rhs = np.zeros(N)
    for (I, k), idx in cid.items(): rhs[idx] = lam * exc(I, k)
    for a in range(nfx):
        for k in range(nz): rhs[fid(a, k)] = lam * exf(a, k)
    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N)).tolil()
    anc = cid[(cxs[0], 0)]
    A[anc, :] = 0; A[anc, anc] = 1.0; rhs[anc] = exc(cxs[0], 0)
    sol = spla.spsolve(A.tocsr(), rhs)
    ec = max(abs(sol[cid[(I, k)]] - exc(I, k)) for I in cxs for k in range(nz))
    ef = max(abs(sol[fid(a, k)] - exf(a, k)) for a in range(nfx) for k in range(nz))
    return float(max(ec, ef))


def solve_composite_hz(f_c, f_f, nc, nz, r, ci0, ci1, cj0, cj1, dzc, dzf,
                       anchor_value=0.0, periodic_h=True, hx=None):
    """**Final-assembly unified operator** for the storm's nest geometry: the horizontal
    composite interface (x,y refined by ``r`` at each z-level, uniform ``hx=1/nc``) plus a
    variable-dz **finite-volume vertical** coupling (walls top/bottom, NOT refined --
    matched z-levels), on mass-flux/pressure variables.

    ``f_c`` (nc,nc,nz), ``f_f`` (nfx,nfy,nz) are the RHS on each level (``nfx=r*(ci1-ci0)``,
    ``nfy=r*(cj1-cj0)``).  ``dzc`` (nz,) cell heights, ``dzf`` (nz-1,) centre spacings.
    ``periodic_h`` wraps the horizontal boundary (the parent); ``False`` gives solid walls
    there.  Returns ``(p_c, p_f)`` with covered coarse cells NaN.  Assembled sparsely,
    solved directly.  ``hx`` is the physical coarse horizontal spacing (default ``1/nc``).
    """
    hc = hx if hx is not None else 1.0 / nc
    hf = hc / r
    nfx, nfy = r * (ci1 - ci0), r * (cj1 - cj0)
    cov = lambda I, J: ci0 <= I < ci1 and cj0 <= J < cj1
    flux = _interface_flux_2d(nc, r, ci0, ci1, cj0, cj1, hx=hc)
    cxy = [(I, J) for I in range(nc) for J in range(nc) if not cov(I, J)]
    cid = {}
    ix = 0
    for (I, J) in cxy:
        for k in range(nz):
            cid[(I, J, k)] = ix; ix += 1
    fb = ix
    fid = lambda a, b, k: fb + (a * nfy + b) * nz + k
    N = fb + nfx * nfy * nz
    gcol = lambda key, k: (fid(key[1], key[2], k) if key[0] == "f" else cid[(key[1], key[2], k)])
    rows, cols, data, rhs = [], [], [], np.zeros(N)
    add = lambda rr, cc, vv: (rows.append(rr), cols.append(cc), data.append(vv))

    def vertical(row, colof, k):                       # variable-dz FV, walls at 0, nz-1
        if k < nz - 1:
            w = 1.0 / (dzf[k] * dzc[k]); add(row, colof(k + 1), w); add(row, row, -w)
        if k > 0:
            w = 1.0 / (dzf[k - 1] * dzc[k]); add(row, colof(k - 1), w); add(row, row, -w)

    for a in range(nfx):                               # fine cells
        for b in range(nfy):
            for k in range(nz):
                row = fid(a, b, k)
                if a == nfx - 1:
                    for key, w in flux("R", b).items(): add(row, gcol(key, k), w / hf)
                else:
                    add(row, fid(a + 1, b, k), 1 / hf ** 2); add(row, row, -1 / hf ** 2)
                if a == 0:
                    for key, w in flux("L", b).items(): add(row, gcol(key, k), -w / hf)
                else:
                    add(row, fid(a - 1, b, k), 1 / hf ** 2); add(row, row, -1 / hf ** 2)
                if b == nfy - 1:
                    for key, w in flux("T", a).items(): add(row, gcol(key, k), w / hf)
                else:
                    add(row, fid(a, b + 1, k), 1 / hf ** 2); add(row, row, -1 / hf ** 2)
                if b == 0:
                    for key, w in flux("B", a).items(): add(row, gcol(key, k), -w / hf)
                else:
                    add(row, fid(a, b - 1, k), 1 / hf ** 2); add(row, row, -1 / hf ** 2)
                vertical(row, lambda kk, a=a, b=b: fid(a, b, kk), k)
                rhs[row] = f_f[a, b, k]

    for (I, J) in cxy:                                 # coarse cells
        for k in range(nz):
            row = cid[(I, J, k)]
            rhs[row] = f_c[I, J, k]
            for (dI, dJ, edge, along) in ((1, 0, "L", J), (-1, 0, "R", J),
                                          (0, 1, "B", I), (0, -1, "T", I)):
                ni, nj = I + dI, J + dJ
                if not periodic_h and (ni < 0 or ni >= nc or nj < 0 or nj >= nc):
                    continue                           # solid-wall Neumann
                In, Jn = ni % nc, nj % nc
                if cov(In, Jn):
                    base = (along - (cj0 if edge in ("L", "R") else ci0)) * r
                    s2 = -(1.0 if edge in ("R", "T") else -1.0)
                    for bb in range(base, base + r):
                        for key, w in flux(edge, bb).items(): add(row, gcol(key, k), s2 * (w / r) / hc)
                else:
                    add(row, cid[(In, Jn, k)], 1 / hc ** 2); add(row, row, -1 / hc ** 2)
            vertical(row, lambda kk, I=I, J=J: cid[(I, J, kk)], k)

    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N)).tolil()
    anc = cid[(cxy[0][0], cxy[0][1], 0)]
    A[anc, :] = 0; A[anc, anc] = 1.0; rhs[anc] = anchor_value
    sol = spla.spsolve(A.tocsr(), rhs)
    p_c = np.full((nc, nc, nz), np.nan)
    for (I, J) in cxy:
        for k in range(nz):
            p_c[I, J, k] = sol[cid[(I, J, k)]]
    p_f = sol[fb:].reshape(nfx, nfy, nz)
    return p_c, p_f


def manufactured_error_hz(nc, nz, r=2, s=1.05, Lz=2.0, periodic_h=True):
    """max error of :func:`solve_composite_hz` vs phi=cos(2pi x)cos(2pi y)cos(pi z/Lz)
    (periodic in x,y, dphi/dz=0 on the z-walls) -- the final unified operator."""
    hc = 1.0 / nc; ci0, ci1, cj0, cj1 = nc // 3, 2 * nc // 3, nc // 3, 2 * nc // 3
    nfx, nfy = r * (ci1 - ci0), r * (cj1 - cj0)
    dz = s ** np.arange(nz); dz *= Lz / dz.sum()
    zf = np.concatenate([[0.0], np.cumsum(dz)]); zc = 0.5 * (zf[:-1] + zf[1:])
    dzc = zf[1:] - zf[:-1]; dzf = np.diff(zc); kz = np.pi / Lz; hf = hc / r
    xc = (np.arange(nc) + 0.5) * hc
    xf = ci0 * hc + (np.arange(nfx) + 0.5) * hf
    yf = cj0 * hc + (np.arange(nfy) + 0.5) * hf
    Cx = np.cos(2 * np.pi * xc); Zc = np.cos(kz * zc)
    ex_c = Cx[:, None, None] * Cx[None, :, None] * Zc[None, None, :]
    ex_f = (np.cos(2 * np.pi * xf)[:, None, None] * np.cos(2 * np.pi * yf)[None, :, None]
            * Zc[None, None, :])
    lam = -(8 * np.pi ** 2 + kz ** 2)
    pc, pf = solve_composite_hz(lam * ex_c, lam * ex_f, nc, nz, r, ci0, ci1, cj0, cj1,
                                dzc, dzf, anchor_value=ex_c[0, 0, 0], periodic_h=periodic_h)
    return float(max(np.nanmax(np.abs(pc - ex_c)), np.abs(pf - ex_f).max()))


def composite_project_massflux_hz(mu_c, mv_c, mw_c, mu_f, mv_f, mw_f, nc, nz, r,
                                  ci0, ci1, cj0, cj1, dzc, dzf, periodic_h=False, hx=None):
    """**Final-assembly projection**: project the storm's full 3-D staggered mass fluxes
    ``m = rho0 u`` (coarse parent + fine nest, horizontal composite interface + variable-dz
    vertical, walls top/bottom) so that ``div(m) = 0`` across the coarse-fine interface,
    combining all three plumbing pieces (wall BC, face bridge, stretched z metric) with the
    unified operator :func:`solve_composite_hz`.

    Arrays (modified in place), storm C-grid convention:
      parent  ``mu_c`` (nc+1,nc,nz), ``mv_c`` (nc,nc+1,nz), ``mw_c`` (nc,nc,nz+1);
      nest    ``mu_f`` (nfx+1,nfy,nz), ``mv_f`` (nfx,nfy+1,nz), ``mw_f`` (nfx,nfy,nz+1),
      ``nfx=r*(ci1-ci0)``, ``nfy=r*(cj1-cj0)`` over coarse cells ``[ci0,ci1)x[cj0,cj1)``.
    ``dzc`` (nz,), ``dzf`` (nz-1,) the vertical metric; ``periodic_h`` wraps the horizontal
    boundary (parent) else solid walls (nest).  The vertical uses solid walls (w=0 at
    ``k=0,nz``).  Returns ``(max|div m| coarse, fine, fine-interface)`` recomputed straight
    from the written-back arrays (self-validating).  In mass-flux variables this is exactly
    the anelastic constraint ``div(rho0 u)=0``; the caller recovers ``u = m/rho0``.
    """
    hc = hx if hx is not None else 1.0 / nc
    hf = hc / r
    nfx, nfy = r * (ci1 - ci0), r * (cj1 - cj0)
    cov = lambda I, J: ci0 <= I < ci1 and cj0 <= J < cj1
    flux = _interface_flux_2d(nc, r, ci0, ci1, cj0, cj1, hx=hc)
    wxp = lambda I: (not periodic_h) and I == nc - 1
    wxm = lambda I: (not periodic_h) and I == 0
    wyp = lambda J: (not periodic_h) and J == nc - 1
    wym = lambda J: (not periodic_h) and J == 0
    anchor = next((I, J) for I in range(nc) for J in range(nc) if not cov(I, J))
    if not periodic_h:
        mu_c[0] = mu_c[nc] = 0.0; mv_c[:, 0] = mv_c[:, nc] = 0.0
    mw_c[:, :, 0] = mw_c[:, :, nz] = 0.0; mw_f[:, :, 0] = mw_f[:, :, nz] = 0.0

    def evalflux(edge, b, k, pc, pf):
        return sum(w * (pf[key[1], key[2], k] if key[0] == "f" else pc[key[1], key[2], k])
                   for key, w in flux(edge, b).items())

    def cfx_at(mu_c, ui, I, J, k):
        if cov((I + 1) % nc, J): base = (J - cj0) * r; return np.mean(ui["L"][base:base + r, k])
        if cov(I, J): base = (J - cj0) * r; return np.mean(ui["R"][base:base + r, k])
        return mu_c[I + 1, J, k]

    def cfy_at(mv_c, ui, I, J, k):
        if cov(I, (J + 1) % nc): base = (I - ci0) * r; return np.mean(ui["B"][base:base + r, k])
        if cov(I, J): base = (I - ci0) * r; return np.mean(ui["T"][base:base + r, k])
        return mv_c[I, J + 1, k]

    def divergence(mu_c, mv_c, mw_c, mu_f, mv_f, mw_f, ui):
        dc = np.full((nc, nc, nz), np.nan); df = np.zeros((nfx, nfy, nz))
        for a in range(nfx):
            for b in range(nfy):
                for k in range(nz):
                    xR = ui["R"][b, k] if a == nfx - 1 else mu_f[a + 1, b, k]
                    xL = ui["L"][b, k] if a == 0 else mu_f[a, b, k]
                    yT = ui["T"][a, k] if b == nfy - 1 else mv_f[a, b + 1, k]
                    yB = ui["B"][a, k] if b == 0 else mv_f[a, b, k]
                    df[a, b, k] = ((xR - xL) + (yT - yB)) / hf + (mw_f[a, b, k + 1] - mw_f[a, b, k]) / dzc[k]
        for I in range(nc):
            for J in range(nc):
                if cov(I, J): continue
                for k in range(nz):
                    fxp = 0.0 if wxp(I) else cfx_at(mu_c, ui, I, J, k)
                    fxm = 0.0 if wxm(I) else cfx_at(mu_c, ui, (I - 1) % nc, J, k)
                    fyp = 0.0 if wyp(J) else cfy_at(mv_c, ui, I, J, k)
                    fym = 0.0 if wym(J) else cfy_at(mv_c, ui, I, (J - 1) % nc, k)
                    dc[I, J, k] = ((fxp - fxm) + (fyp - fym)) / hc + (mw_c[I, J, k + 1] - mw_c[I, J, k]) / dzc[k]
        return dc, df

    ui = {"L": mu_f[0].copy(), "R": mu_f[nfx].copy(), "B": mv_f[:, 0].copy(), "T": mv_f[:, nfy].copy()}
    fc, ff = divergence(mu_c, mv_c, mw_c, mu_f, mv_f, mw_f, ui)
    pc, pf = solve_composite_hz(np.nan_to_num(fc), ff, nc, nz, r, ci0, ci1, cj0, cj1,
                                dzc, dzf, anchor_value=0.0, periodic_h=periodic_h, hx=hc)
    pc = np.nan_to_num(pc)

    mu_f[1:nfx] -= (pf[1:] - pf[:-1]) / hf
    mv_f[:, 1:nfy] -= (pf[:, 1:] - pf[:, :-1]) / hf
    mw_f[:, :, 1:nz] -= (pf[:, :, 1:] - pf[:, :, :-1]) / dzf[None, None, :]
    for k in range(nz):
        for b in range(nfy):
            ui["L"][b, k] -= evalflux("L", b, k, pc, pf); ui["R"][b, k] -= evalflux("R", b, k, pc, pf)
        for a in range(nfx):
            ui["B"][a, k] -= evalflux("B", a, k, pc, pf); ui["T"][a, k] -= evalflux("T", a, k, pc, pf)
    mu_f[0] = ui["L"]; mu_f[nfx] = ui["R"]; mv_f[:, 0] = ui["B"]; mv_f[:, nfy] = ui["T"]
    mw_c[:, :, 1:nz] -= (pc[:, :, 1:] - pc[:, :, :-1]) / dzf[None, None, :]
    for I in range(nc):
        for J in range(nc):
            for k in range(nz):
                if not cov(I, J) and not cov((I + 1) % nc, J) and not wxp(I):
                    mu_c[I + 1, J, k] -= (pc[(I + 1) % nc, J, k] - pc[I, J, k]) / hc
                if not cov(I, J) and not cov(I, (J + 1) % nc) and not wyp(J):
                    mv_c[I, J + 1, k] -= (pc[I, (J + 1) % nc, k] - pc[I, J, k]) / hc
    for J in range(cj0, cj1):                          # reflux parent interface faces
        base = (J - cj0) * r
        mu_c[ci0, J] = ui["L"][base:base + r].mean(axis=0); mu_c[ci1, J] = ui["R"][base:base + r].mean(axis=0)
    for I in range(ci0, ci1):
        base = (I - ci0) * r
        mv_c[I, cj0] = ui["B"][base:base + r].mean(axis=0); mv_c[I, cj1] = ui["T"][base:base + r].mean(axis=0)
    if periodic_h:
        mu_c[0] = mu_c[nc]; mv_c[:, 0] = mv_c[:, nc]

    df_chk = ((mu_f[1:] - mu_f[:-1]) / hf + (mv_f[:, 1:] - mv_f[:, :-1]) / hf
              + (mw_f[:, :, 1:] - mw_f[:, :, :-1]) / dzc[None, None, :])
    dc_chk = np.full((nc, nc, nz), np.nan)
    for I in range(nc):
        for J in range(nc):
            if cov(I, J): continue
            for k in range(nz):
                fxp = 0.0 if wxp(I) else mu_c[I + 1, J, k]; fxm = 0.0 if wxm(I) else mu_c[I, J, k]
                fyp = 0.0 if wyp(J) else mv_c[I, J + 1, k]; fym = 0.0 if wym(J) else mv_c[I, J, k]
                dc_chk[I, J, k] = ((fxp - fxm) + (fyp - fym)) / hc + (mw_c[I, J, k + 1] - mw_c[I, J, k]) / dzc[k]
    dc_chk[anchor[0], anchor[1], 0] = 0.0
    border = np.zeros((nfx, nfy, nz), bool)
    border[0] = border[-1] = True; border[:, 0] = border[:, -1] = True
    return (float(np.nanmax(np.abs(dc_chk))), float(np.abs(df_chk).max()),
            float(np.abs(df_chk[border]).max()))


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
           "manufactured_error_2d_wall", "solve_3d", "manufactured_error_3d",
           "project_divergence_2d", "project_divergence_3d",
           "composite_project_massflux_2d", "manufactured_error_metric_z",
           "solve_composite_hz", "manufactured_error_hz", "composite_project_massflux_hz"]


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

    print("FINAL ASSEMBLY -- unified operator (horizontal interface + stretched-z), manufactured:")
    for s in (1.0, 1.05):
        prev = None
        for nc in (12, 24):
            err = manufactured_error_hz(nc, nc, s=s)
            ratio = "" if prev is None else "  ratio=%.2f" % (prev / err)
            print("  s=%.2f nc=%3d nz=%3d  max|err|=%.3e%s" % (s, nc, nc, err, ratio)); prev = err
    print("FINAL ASSEMBLY -- full 3-D storm mass-flux projection (div m -> 0 across interface):")
    for periodic_h in (False, True):
        nc, nz, rr = 12, 8, 2
        ci0, ci1, cj0, cj1 = nc // 3, 2 * nc // 3, nc // 3, 2 * nc // 3
        nfx, nfy = rr * (ci1 - ci0), rr * (cj1 - cj0)
        dz = 1.05 ** np.arange(nz); dz *= 2.0 / dz.sum()
        zf = np.concatenate([[0.0], np.cumsum(dz)]); zc = 0.5 * (zf[:-1] + zf[1:])
        dzc = zf[1:] - zf[:-1]; dzf = np.diff(zc)
        rng = np.random.default_rng(0)
        mu_c = rng.standard_normal((nc + 1, nc, nz)); mv_c = rng.standard_normal((nc, nc + 1, nz))
        mw_c = rng.standard_normal((nc, nc, nz + 1))
        mu_f = rng.standard_normal((nfx + 1, nfy, nz)); mv_f = rng.standard_normal((nfx, nfy + 1, nz))
        mw_f = rng.standard_normal((nfx, nfy, nz + 1))
        dc, df, di = composite_project_massflux_hz(mu_c, mv_c, mw_c, mu_f, mv_f, mw_f,
                                                   nc, nz, rr, ci0, ci1, cj0, cj1, dzc, dzf,
                                                   periodic_h=periodic_h)
        print("  periodic_h=%-5s  max|div m|: coarse=%.2e  fine=%.2e  interface=%.2e"
              % (periodic_h, dc, df, di))
