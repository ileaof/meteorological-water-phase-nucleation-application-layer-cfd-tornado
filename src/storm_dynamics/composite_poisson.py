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


__all__ = ["solve_1d", "manufactured_error"]


if __name__ == "__main__":
    print("composite 2-level 1-D Poisson (2nd-order interface stencil):")
    prev = None
    for nc in (48, 96, 192, 384):
        err, hc = manufactured_error(nc)
        ratio = "" if prev is None else "  ratio=%.2f" % (prev / err)
        print("  nc=%3d  hc=%.4f  max|err|=%.3e%s" % (nc, hc, err, ratio))
        prev = err
