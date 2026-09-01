"""General low-memory anelastic pressure projection (ROADMAP §3f increment 3) — the
**non-separable fallback** for the FFT+tridiag solver.

Assembles the sparse anelastic Poisson matrix and solves it with a **Jacobi-preconditioned
conjugate gradient** (scipy ``cg``).  Two properties the other solvers lack:

* **General** — no separability assumption.  CG works on *any* sparse SPD operator, so when
  the model gains terrain-following coordinates, horizontal grid stretching, or x,y-varying
  coefficients (the forecast-model generalisation), only the matrix assembly changes; the
  solver is unchanged.  (The assembly here is the current uniform-horizontal / variable-dz /
  ``rho0(z)`` operator — the same one FFT+tridiag solves, so the two cross-check.)
* **Low memory** — CG stores a handful of vectors, no factorisation at all: far below the
  full direct ``splu`` that OOMs a fine nest, and below ILU too.

Two subtleties make CG actually converge here (both were bugs during development):

* **Sign.**  The discrete Laplacian ``L = div grad`` is *negative* definite; CG requires an
  SPD matrix, so we assemble and solve the **negated** operator ``-L`` (``lap(phi)=f`` becomes
  ``(-L) phi = -f``; the recovered ``phi`` is identical).  On the wrong-sign operator CG
  stalls — which is exactly why the plain Jacobi-CG in ``meteorological_flow.PressureSolver``
  "diverges" on the stretched grid; the operator here is symmetric *and* positive-definite.
* **Preconditioner must be SPD.**  ``scipy.sparse.linalg.spilu`` (ILU) is *not* symmetric, so
  using it as a CG preconditioner breaks the Krylov recursion (CG then stalls at ~1e-3).  We
  use a **Jacobi (diagonal)** preconditioner instead — SPD, free, and it trims iterations.
  For very large stiff problems an algebraic-multigrid preconditioner (``pyamg``, optional)
  would cut iterations further; Jacobi-CG already reaches ~1e-11 in a fraction of a second on
  uniform and stretched (s up to 1.10) grids, periodic or wall.

The null space (periodic / pure-Neumann) is removed by a **symmetric pin** of cell 0 (zero
its row and column, unit diagonal), keeping the matrix SPD.  Self-contained (scipy.sparse);
does not touch ``PressureSolver``.  Use it wherever a grid is *not* separable (FFT+tridiag is
the fast path when it is).
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .pressure_fft import anelastic_divergence, _z_tridiag


def build_anelastic_matrix(nx, ny, nz, dx, dy, dzc, dzf, periodic_h):
    """Sparse **negated** anelastic Poisson operator ``-L`` (constant-coeff horizontal +
    variable-dz vertical, walls in z; periodic or wall in x,y), assembled **SPD**: each row is
    volume-weighted by the cell height ``dzc_k`` (the variable-dz FV Laplacian is self-adjoint
    only in the volume inner product) and the whole operator is negated so the diagonal is
    positive and the matrix is positive-definite — CG needs SPD.  Symmetrically pinned at cell
    0 (row+col zeroed, unit diagonal) to remove the constant null space.  Solve it against the
    **negated, volume-weighted RHS** ``-dzc_k * f``.  This is the seam a terrain /
    variable-coefficient assembly would replace."""
    lower, diagz, upper = _z_tridiag(dzc, dzf)
    idx = lambda i, j, k: (i * ny + j) * nz + k
    N = nx * ny * nz
    hx2, hy2 = 1.0 / dx ** 2, 1.0 / dy ** 2
    rows, cols, data = [], [], []
    add = lambda r, c, v: (rows.append(r), cols.append(c), data.append(v))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                r = idx(i, j, k); vol = dzc[k]; diag = 0.0    # row scaled by the cell volume
                for di in (-1, 1):
                    ii = (i + di) % nx if periodic_h else i + di
                    if not periodic_h and (ii < 0 or ii >= nx):
                        continue                              # wall Neumann: drop the face
                    add(r, idx(ii, j, k), -hx2 * vol); diag -= hx2   # negated: off-diag < 0
                for dj in (-1, 1):
                    jj = (j + dj) % ny if periodic_h else j + dj
                    if not periodic_h and (jj < 0 or jj >= ny):
                        continue
                    add(r, idx(i, jj, k), -hy2 * vol); diag -= hy2
                if k < nz - 1:
                    add(r, idx(i, j, k + 1), -upper[k] * vol)  # = -1/dzf[k]  (symmetric)
                if k > 0:
                    add(r, idx(i, j, k - 1), -lower[k] * vol)  # = -1/dzf[k-1]
                add(r, r, -(diag + diagz[k]) * vol)            # positive diagonal (SPD)
    A = sp.csr_matrix((data, (rows, cols)), shape=(N, N)).tolil()
    A[0, :] = 0; A[:, 0] = 0; A[0, 0] = 1.0                   # symmetric pin (phi[0]=0), +1 diag
    return A.tocsr()


def project_anelastic_iterative(u, v, w, rho0_c, rho0_wface, dx, dy, dzc, dzf,
                                periodic_h=True, tol=1e-10, maxiter=20000):
    """Make the staggered velocity anelastically divergence-free with a Jacobi-preconditioned
    CG solve of the SPD operator ``-L`` (general, low-memory).  Same signature/semantics as
    :func:`pressure_fft.project_anelastic_fft`; modifies ``u,v,w`` in place, returns
    ``max|div(rho0 u)|`` (should reach the CG tolerance)."""
    rc, rw = np.asarray(rho0_c), np.asarray(rho0_wface)
    nx, ny, nz = u.shape[0] - 1, v.shape[1] - 1, len(dzc)
    w[:, :, 0] = 0.0; w[:, :, -1] = 0.0
    if periodic_h:
        u[-1] = u[0]; v[:, -1] = v[:, 0]
    else:
        u[0] = u[-1] = 0.0; v[:, 0] = v[:, -1] = 0.0
    # negated, volume-weighted RHS (-dzc_k * f) to match the SPD, volume-weighted operator -L
    f = (-anelastic_divergence(u, v, w, rc, rw, dx, dy, dzc)
         * np.asarray(dzc)[None, None, :]).reshape(-1).copy()
    f[0] = 0.0                                                # consistent with the pinned row
    A = build_anelastic_matrix(nx, ny, nz, dx, dy, dzc, dzf, periodic_h)
    d = A.diagonal()                                          # SPD Jacobi (diagonal) precond
    M = spla.LinearOperator(A.shape, lambda x: x / d)
    phi, _ = spla.cg(A, f, M=M, rtol=tol, atol=0.0, maxiter=maxiter)
    phi = phi.reshape(nx, ny, nz)
    gz = (phi[:, :, 1:] - phi[:, :, :-1]) / np.asarray(dzf)[None, None, :]
    if periodic_h:
        u[1:-1] -= (phi[1:] - phi[:-1]) / dx / rc[None, None, :]
        u[0] -= (phi[0] - phi[-1]) / dx / rc[None, :]; u[-1] = u[0]
        v[:, 1:-1] -= (phi[:, 1:] - phi[:, :-1]) / dy / rc[None, None, :]
        v[:, 0] -= (phi[:, 0] - phi[:, -1]) / dy / rc[None, :]; v[:, -1] = v[:, 0]
    else:
        u[1:-1] -= (phi[1:] - phi[:-1]) / dx / rc[None, None, :]
        v[:, 1:-1] -= (phi[:, 1:] - phi[:, :-1]) / dy / rc[None, None, :]
    w[:, :, 1:-1] -= gz / rw[None, None, 1:-1]
    return float(np.abs(anelastic_divergence(u, v, w, rc, rw, dx, dy, dzc)).max())
