"""Geometric multigrid Poisson solver (the "bring-our-own-solver" for AMR).

pyAMReX does not expose AMReX's `MLMG` composite multigrid to Python (see
`docs/amr_design.md`), so the AMR pressure/projection solve has to come from us.
This is the solver **kernel**: a 2-D, cell-centred, periodic geometric multigrid
V-cycle (red-black Gauss-Seidel smoother, full-weighting restriction, cell-centred
bilinear prolongation) for ``lap(phi) = f``.

Verified (see ``test_poisson_multigrid_*``) against a manufactured solution
``phi = sin(2pi x) sin(2pi y)``:

* **h-independent convergence** -- the same handful of V-cycles drives the residual
  to ~1e-11 at any resolution (the multigrid property);
* **2nd-order accuracy** -- ``max|phi - phi_exact|`` falls by 4x when the grid is
  refined 2x (error ∝ h^2).

The **composite** (2-level, coarse+fine) version -- this smoother/V-cycle on each
level plus the coarse-fine interface coupling (matched flux + the refluxing already
in :mod:`storm_dynamics.amr`) -- is the remaining AMR-projection step; this kernel
is its building block.
"""
from __future__ import annotations

import numpy as np


def laplacian(phi, h):
    """Periodic 5-point Laplacian of a cell-centred field."""
    return (np.roll(phi, 1, 0) + np.roll(phi, -1, 0) +
            np.roll(phi, 1, 1) + np.roll(phi, -1, 1) - 4 * phi) / h ** 2


def _smooth(phi, f, h, n=2):
    """Red-black Gauss-Seidel relaxation of lap(phi)=f (periodic)."""
    idx_parity = (np.indices(phi.shape).sum(0) % 2)
    for _ in range(n):
        for color in (0, 1):
            nb = (np.roll(phi, 1, 0) + np.roll(phi, -1, 0) +
                  np.roll(phi, 1, 1) + np.roll(phi, -1, 1))
            phi = np.where(idx_parity == color, (nb - h ** 2 * f) / 4.0, phi)
    return phi


def _restrict(r):
    """Full-weighting restriction (2x2 average) fine -> coarse."""
    n = r.shape[0]
    return 0.25 * (r[0:n:2, 0:n:2] + r[1:n:2, 0:n:2] + r[0:n:2, 1:n:2] + r[1:n:2, 1:n:2])


def _prolong(e):
    """Cell-centred bilinear prolongation coarse -> fine (periodic)."""
    n = e.shape[0]
    out = np.empty((2 * n, 2 * n))
    eL = np.roll(e, 1, 0); eR = np.roll(e, -1, 0)
    eD = np.roll(e, 1, 1); eU = np.roll(e, -1, 1)
    eLD = np.roll(eL, 1, 1); eLU = np.roll(eL, -1, 1)
    eRD = np.roll(eR, 1, 1); eRU = np.roll(eR, -1, 1)
    a, b, c = 9 / 16, 3 / 16, 1 / 16
    out[0::2, 0::2] = a * e + b * eL + b * eD + c * eLD
    out[1::2, 0::2] = a * e + b * eR + b * eD + c * eRD
    out[0::2, 1::2] = a * e + b * eL + b * eU + c * eLU
    out[1::2, 1::2] = a * e + b * eR + b * eU + c * eRU
    return out


def _vcycle(phi, f, h, npre=3, npost=3):
    phi = _smooth(phi, f, h, npre)
    if phi.shape[0] <= 4:                    # coarsest: smooth hard, pin the mean
        for _ in range(30):
            phi = _smooth(phi, f, h, 1); phi -= phi.mean()
        return phi
    res = f - laplacian(phi, h)
    corr = _vcycle(np.zeros_like(_restrict(res)), _restrict(res), 2 * h, npre, npost)
    phi = phi + _prolong(corr)
    phi = _smooth(phi, f, h, npost)
    return phi - phi.mean()                  # periodic null space: pin the mean


def solve(f, h, ncyc=30, tol=1e-9):
    """Solve ``lap(phi) = f`` on a periodic n x n grid (n a power of two).

    ``f`` must have zero mean (periodic compatibility).  Returns ``(phi, res_hist)``
    where ``res_hist`` is the relative residual after each V-cycle.
    """
    f = f - f.mean()
    phi = np.zeros_like(f)
    r0 = np.linalg.norm(f - laplacian(phi, h)) or 1.0
    hist = []
    for _ in range(ncyc):
        phi = _vcycle(phi, f, h)
        rel = np.linalg.norm(f - laplacian(phi, h)) / r0
        hist.append(rel)
        if rel < tol:
            break
    return phi, hist


__all__ = ["solve", "laplacian"]


if __name__ == "__main__":
    print("Geometric multigrid Poisson -- manufactured-solution check")
    for n in (64, 128, 256):
        h = 1.0 / n
        xs = (np.arange(n) + 0.5) * h
        X, Y = np.meshgrid(xs, xs, indexing="ij")
        exact = np.sin(2 * np.pi * X) * np.sin(2 * np.pi * Y)
        phi, hist = solve(-8 * np.pi ** 2 * exact, h)
        err = np.abs((phi - phi.mean()) - (exact - exact.mean())).max()
        print("  n=%3d  V-cycles=%2d  res=%.1e  max err=%.2e  (h^2=%.1e)"
              % (n, len(hist), hist[-1], err, h ** 2))
