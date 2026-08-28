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

This 1-D kernel nails the normal-direction interface stencil; the 2-D/3-D extension
applies the same normal stencil with a tangential coarse interpolation for the
ghost -- the remaining step before wiring it into the anelastic projection.
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


__all__ = ["solve_1d", "manufactured_error", "solve_2d", "manufactured_error_2d"]


if __name__ == "__main__":
    for dim, fn, grids in (("1-D", manufactured_error, (48, 96, 192, 384)),
                           ("2-D", manufactured_error_2d, (24, 48, 96))):
        print("composite 2-level %s Poisson (2nd-order interface stencil):" % dim)
        prev = None
        for nc in grids:
            err, hc = fn(nc)
            ratio = "" if prev is None else "  ratio=%.2f" % (prev / err)
            print("  nc=%3d  hc=%.4f  max|err|=%.3e%s" % (nc, hc, err, ratio))
            prev = err
