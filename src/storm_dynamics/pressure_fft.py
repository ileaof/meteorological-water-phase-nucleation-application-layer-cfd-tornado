"""Low-memory anelastic pressure projection (ROADMAP §3f) — the fine-AMR-nest OOM fix.

The storm's grid is **uniform in x,y and stretched in z**, so the pressure Poisson
``lap(phi) = div(rho0 u*)`` separates: a transform in (x,y) diagonalises the horizontal
Laplacian, leaving, for each horizontal wavenumber, a **tridiagonal** system in z (the
variable-dz vertical operator shifted by the horizontal eigenvalue).  Solved by a batched
Thomas sweep — ``O(N)`` memory and ``O(N log N)`` time, **no stored LU factorisation**, so a
fine nest (e.g. 48³) that OOMs the direct ``splu`` solve fits easily.

Transform: **FFT** for a periodic-horizontal domain (the parent), **DCT-II** for
solid-wall / free-slip horizontal walls (the nest) — the cell-centred Neumann basis.
Vertical: solid walls (``w=0`` top/bottom), i.e. Neumann for ``phi``.

The projection is exact for the discrete operator: ``u = u* - grad(phi)/rho0`` makes
``div(rho0 u) = div(rho0 u*) - lap(phi) = 0`` to the transform's round-off.  This module is
self-contained (NumPy + scipy.fft); it does **not** modify ``meteorological_flow``'s
``PressureSolver``.  Wire it into a nest via :func:`project_anelastic_fft`.

**Assumptions / generality (see docs/ROADMAP.md §3f).**  Exact but only for the *separable*
case: **uniform horizontal spacing** (constant dx, dy), coefficients **homogeneous in x,y**
(anelastic ``rho0 = rho0(z)`` only — the z-tridiagonal is the same for every column; dz(z)
and rho0(z) may vary in z, not in x,y), and **separable homogeneous BCs** (periodic -> FFT,
uniform wall -> DCT).  The current storm and all its AMR nests are in this class (no loss).
It does **not** handle terrain-following coordinates, horizontal grid stretching, map
projections, x,y-varying reference states, or irregular/immersed boundaries — those need the
general direct / multigrid-preconditioned solvers; select this one only when ``grid`` is
separable, never as a silent default.
"""
from __future__ import annotations

import numpy as np

try:
    from scipy.fft import dctn, idctn
    _HAVE_DCT = True
except Exception:                                    # pragma: no cover
    _HAVE_DCT = False


def _z_tridiag(dzc, dzf):
    """Variable-dz vertical finite-volume operator (Neumann walls) as (lower, diag, upper).

    ``(L_z phi)_k = [ (phi_{k+1}-phi_k)/dzf_k - (phi_k-phi_{k-1})/dzf_{k-1} ] / dzc_k`` with
    no flux through the top/bottom walls.  ``dzc`` (nz,) cell heights, ``dzf`` (nz-1,) centre
    spacings."""
    nz = len(dzc)
    lower = np.zeros(nz); upper = np.zeros(nz); diag = np.zeros(nz)
    for k in range(nz):
        if k < nz - 1:
            u = 1.0 / (dzf[k] * dzc[k]); upper[k] = u; diag[k] -= u
        if k > 0:
            l = 1.0 / (dzf[k - 1] * dzc[k]); lower[k] = l; diag[k] -= l
    return lower, diag, upper


def _thomas_batch(lower, diag, upper, rhs):
    """Solve a batch of tridiagonals sharing off-diagonals: ``lower``/``upper`` (nz,),
    ``diag`` (B, nz), ``rhs`` (B, nz) -> x (B, nz).  Vectorised over the batch B."""
    nz = diag.shape[1]
    cp = np.empty_like(diag); dp = np.empty_like(rhs)
    cp[:, 0] = upper[0] / diag[:, 0]; dp[:, 0] = rhs[:, 0] / diag[:, 0]
    for k in range(1, nz):
        m = diag[:, k] - lower[k] * cp[:, k - 1]
        cp[:, k] = (upper[k] / m) if k < nz - 1 else 0.0
        dp[:, k] = (rhs[:, k] - lower[k] * dp[:, k - 1]) / m
    x = np.empty_like(rhs); x[:, -1] = dp[:, -1]
    for k in range(nz - 2, -1, -1):
        x[:, k] = dp[:, k] - cp[:, k] * x[:, k + 1]
    return x


def _solve_poisson(f, dx, dy, dzc, dzf, periodic_h):
    """Solve ``lap(phi) = f`` (uniform x,y + variable-dz z, Neumann z-walls) by transform +
    batched Thomas.  ``periodic_h`` -> FFT; else DCT-II (wall).  Returns phi (nx,ny,nz)."""
    nx, ny, nz = f.shape
    lower, diagz, upper = _z_tridiag(dzc, dzf)
    if periodic_h:
        fh = np.fft.fft2(f, axes=(0, 1))
        lamx = -2.0 * (1.0 - np.cos(2 * np.pi * np.arange(nx) / nx)) / dx ** 2
        lamy = -2.0 * (1.0 - np.cos(2 * np.pi * np.arange(ny) / ny)) / dy ** 2
    else:
        if not _HAVE_DCT:
            raise RuntimeError("scipy.fft.dctn required for wall (non-periodic) domains")
        fh = dctn(f, axes=(0, 1), type=2, norm="ortho")
        lamx = -2.0 * (1.0 - np.cos(np.pi * np.arange(nx) / nx)) / dx ** 2
        lamy = -2.0 * (1.0 - np.cos(np.pi * np.arange(ny) / ny)) / dy ** 2
    shift = (lamx[:, None] + lamy[None, :]).reshape(-1)          # (nx*ny,)
    rhs = fh.reshape(nx * ny, nz)
    diag = diagz[None, :] + shift[:, None]                       # (B, nz)
    # the shift==0 mode(s) (the (0,0) horizontal mean) give the pure vertical Neumann
    # operator L_z, singular ONLY in the additive constant -> solve pinned, don't zero.
    sing = np.abs(shift) < 1e-12 * (np.abs(diagz).max() + 1e-30)
    diag_safe = diag.copy(); diag_safe[sing, :] = 1.0           # dummy (overwritten below)
    ph = _thomas_batch(lower, diag_safe, upper, rhs)
    if sing.any():                                              # pinned vertical solve: phi[0]=0
        di = diagz.copy(); up = upper.copy(); lo = lower
        di[0] = 1.0; up[0] = 0.0
        A = np.diag(di) + np.diag(up[:-1], 1) + np.diag(lo[1:], -1)
        for idx in np.where(sing)[0]:
            rr = rhs[idx].copy(); rr[0] = 0.0
            ph[idx] = np.linalg.solve(A, rr)
    ph = ph.reshape(nx, ny, nz)
    if periodic_h:
        return np.real(np.fft.ifft2(ph, axes=(0, 1)))
    return idctn(ph, axes=(0, 1), type=2, norm="ortho")


def anelastic_divergence(u, v, w, rho0_c, rho0_wface, dx, dy, dzc):
    """``div(rho0 u)`` at cell centres from the staggered mass flux (rho0=rho0(z))."""
    rc = np.asarray(rho0_c)[None, None, :]; rw = np.asarray(rho0_wface)[None, None, :]
    return (rc * (u[1:, :, :] - u[:-1, :, :]) / dx
            + rc * (v[:, 1:, :] - v[:, :-1, :]) / dy
            + (rw[:, :, 1:] * w[:, :, 1:] - rw[:, :, :-1] * w[:, :, :-1]) / np.asarray(dzc)[None, None, :])


def project_anelastic_fft(u, v, w, rho0_c, rho0_wface, dx, dy, dzc, dzf, periodic_h=True):
    """Make the staggered velocity anelastically divergence-free, ``div(rho0 u)=0``, with the
    **low-memory** FFT/DCT + tridiagonal solve.  ``u`` (nx+1,ny,nz), ``v`` (nx,ny+1,nz),
    ``w`` (nx,ny,nz+1) (host NumPy, modified in place); ``rho0_c`` (nz,), ``rho0_wface``
    (nz+1,); ``dzc`` (nz,), ``dzf`` (nz-1,).  ``periodic_h`` wraps x,y (parent) else solid
    walls (nest).  Returns ``max|div(rho0 u)|`` after projection (should be ~round-off)."""
    rc = np.asarray(rho0_c); rw = np.asarray(rho0_wface)
    w[:, :, 0] = 0.0; w[:, :, -1] = 0.0                         # solid top/bottom
    if periodic_h:                                             # sync the redundant wrap face
        u[-1, :, :] = u[0, :, :]; v[:, -1, :] = v[:, 0, :]
    else:                                                      # solid lateral walls
        u[0, :, :] = u[-1, :, :] = 0.0; v[:, 0, :] = v[:, -1, :] = 0.0
    f = anelastic_divergence(u, v, w, rc, rw, dx, dy, dzc)
    phi = _solve_poisson(f, dx, dy, dzc, dzf, periodic_h)
    # correct: u -= grad(phi)/rho0   (rho0 at the face)
    gx = (phi[1:, :, :] - phi[:-1, :, :]) / dx
    gy = (phi[:, 1:, :] - phi[:, :-1, :]) / dy
    gz = (phi[:, :, 1:] - phi[:, :, :-1]) / np.asarray(dzf)[None, None, :]
    if periodic_h:
        u[1:-1, :, :] -= gx / rc[None, None, :]
        u[0, :, :] -= (phi[0, :, :] - phi[-1, :, :]) / dx / rc[None, :]        # periodic wrap
        u[-1, :, :] = u[0, :, :]
        v[:, 1:-1, :] -= gy / rc[None, None, :]
        v[:, 0, :] -= (phi[:, 0, :] - phi[:, -1, :]) / dy / rc[None, :]
        v[:, -1, :] = v[:, 0, :]
    else:
        u[1:-1, :, :] -= gx / rc[None, None, :]                 # wall faces stay 0
        v[:, 1:-1, :] -= gy / rc[None, None, :]
    # interior z-faces (k = 1..nz-1): w -= (dphi/dz)/rho0_wface  (walls at k=0,nz stay 0)
    w[:, :, 1:-1] -= gz / rw[None, None, 1:-1]
    return float(np.abs(anelastic_divergence(u, v, w, rc, rw, dx, dy, dzc)).max())
