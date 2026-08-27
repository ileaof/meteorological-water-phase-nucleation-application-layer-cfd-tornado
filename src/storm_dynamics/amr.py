"""AMR Milestone 1 -- flux-conservative **refluxing** on a static 2-level hierarchy.

This is a small, self-contained **reference implementation** of the Berger-Colella
flux-register correction (`docs/amr_design.md`, §2.2), for flux-form advection on a
coarse periodic grid with one refined patch.  It is *not* the storm solver -- it is
the algorithm that makes a coarse-fine interface **conserve**, verified here in
pure NumPy so the logic is nailed down before it is ported onto a framework's
`FluxRegister` (AMReX/Chombo).

The demonstration (see :func:`demo` and ``test_amr_refluxing_conserves``):

* advect a passive tracer with a uniform velocity on the coarse grid and, inside a
  patch, on a grid ``refine`` times finer, sub-cycled ``refine`` times;
* **without** refluxing the total mass drifts (coarse and fine fluxes across the
  interface disagree -> an interface leak);
* **with** refluxing the flux mismatch on the coarse-fine faces is applied back to
  the coarse cells just outside the patch, restoring conservation to **machine
  precision**.

Restriction (average-down of the covered cells) is the companion operator, already
in :func:`storm_dynamics.nesting.conservative_restrict`.  What remains for full AMR
(a multilevel Poisson solve and adaptive regridding) is in ``docs/amr_design.md``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def _upwind(q, u, v, dt, h, periodic=True, ghost=None):
    """1st-order upwind flux-form advection (u, v > 0).  Returns (q_new, Fx, Fy)
    with face fluxes Fx (nx+1, ny), Fy (nx, ny+1).  Non-periodic patches take the
    incoming (left column, bottom row) values from ``ghost``."""
    nx, ny = q.shape
    Fx = np.zeros((nx + 1, ny)); Fy = np.zeros((nx, ny + 1))
    Fx[1:, :] = u * q                                  # face i+1/2 <- cell i (u>0)
    Fx[0, :] = u * (q[-1, :] if periodic else ghost[0])
    Fy[:, 1:] = v * q
    Fy[:, 0] = v * (q[:, -1] if periodic else ghost[1])
    qn = q - dt / h * (Fx[1:, :] - Fx[:-1, :]) - dt / h * (Fy[:, 1:] - Fy[:, :-1])
    return qn, Fx, Fy


@dataclass
class TwoLevelReflux:
    """Static 2-level advection with optional Berger-Colella refluxing.

    A coarse ``N x N`` periodic grid (spacing ``dx``) with one patch over coarse
    cells ``[ci0, ci1) x [cj0, cj1)`` refined by ``refine`` and sub-cycled
    ``refine`` times per coarse step.  Uniform velocity ``(U, V) > 0``.
    """
    N: int = 24
    dx: float = 1.0
    U: float = 1.0
    V: float = 1.0
    refine: int = 3
    ci0: int = 8
    ci1: int = 16
    cj0: int = 8
    cj1: int = 16
    cfl: float = 0.4
    qc: np.ndarray = field(default=None, repr=False)
    qf: np.ndarray = field(default=None, repr=False)

    def __post_init__(self):
        r = self.refine
        i = np.arange(self.N)[:, None]; j = np.arange(self.N)[None, :]
        self.qc = 1.0 + 0.3 * np.sin(2 * np.pi * i / self.N) * np.cos(2 * np.pi * j / self.N)
        self.qf = np.repeat(np.repeat(self.qc[self.ci0:self.ci1, self.cj0:self.cj1],
                                      r, axis=0), r, axis=1).copy()
        self.dxf = self.dx / r
        self.dtc = self.cfl * self.dx / max(self.U, self.V)

    # ---- total mass = valid coarse (outside patch) + fine patch ----
    def total_mass(self) -> float:
        mask = np.ones_like(self.qc, bool)
        mask[self.ci0:self.ci1, self.cj0:self.cj1] = False
        return float(self.qc[mask].sum() * self.dx ** 2 + self.qf.sum() * self.dxf ** 2)

    def step(self, reflux: bool = True) -> None:
        r, dx, dxf = self.refine, self.dx, self.dxf
        ci0, ci1, cj0, cj1 = self.ci0, self.ci1, self.cj0, self.cj1
        # coarse step over the whole grid; record coarse interface fluxes
        qc_new, Fxc, Fyc = _upwind(self.qc, self.U, self.V, self.dtc, dx, periodic=True)
        reg_L = Fxc[ci0, cj0:cj1] * self.dtc * dx
        reg_R = Fxc[ci1, cj0:cj1] * self.dtc * dx
        reg_B = Fyc[ci0:ci1, cj0] * self.dtc * dx
        reg_T = Fyc[ci0:ci1, cj1] * self.dtc * dx
        # fine sub-steps; accumulate fine interface fluxes at coarse-cell granularity
        dtf = self.dtc / r
        fL = np.zeros(ci1 - ci0); fR = np.zeros(ci1 - ci0)
        fB = np.zeros(ci1 - ci0); fT = np.zeros(ci1 - ci0)
        for _ in range(r):
            gl = np.repeat(self.qc[ci0 - 1, cj0:cj1], r)
            gb = np.repeat(self.qc[ci0:ci1, cj0 - 1], r)
            self.qf, Fxf, Fyf = _upwind(self.qf, self.U, self.V, dtf, dxf,
                                        periodic=False, ghost=(gl, gb))
            fL += Fxf[0, :].reshape(ci1 - ci0, r).sum(1) * dtf * dxf
            fR += Fxf[-1, :].reshape(ci1 - ci0, r).sum(1) * dtf * dxf
            fB += Fyf[:, 0].reshape(ci1 - ci0, r).sum(1) * dtf * dxf
            fT += Fyf[:, -1].reshape(ci1 - ci0, r).sum(1) * dtf * dxf
        self.qc = qc_new
        # average-down (conservative restriction) onto the covered coarse cells
        self.qc[ci0:ci1, cj0:cj1] = self.qf.reshape(ci1 - ci0, r, cj1 - cj0, r).mean(axis=(1, 3))
        # reflux: the flux mismatch on each interface face corrects the coarse cell
        # just OUTSIDE the patch, so the mass leaving the coarse side equals the mass
        # entering the fine side -> conservation.
        if reflux:
            Vc = dx ** 2
            self.qc[ci0 - 1, cj0:cj1] += (reg_L - fL) / Vc
            self.qc[ci1, cj0:cj1] -= (reg_R - fR) / Vc
            self.qc[ci0:ci1, cj0 - 1] += (reg_B - fB) / Vc
            self.qc[ci0:ci1, cj1] -= (reg_T - fT) / Vc

    def run(self, nsteps: int = 40, reflux: bool = True) -> float:
        """Run ``nsteps`` and return the relative total-mass drift |ΔM|/M0."""
        m0 = self.total_mass()
        for _ in range(nsteps):
            self.step(reflux=reflux)
        return abs(self.total_mass() - m0) / abs(m0)


def demo(nsteps: int = 40) -> dict:
    """Relative mass drift with and without refluxing (fresh state each)."""
    return {"no_reflux": TwoLevelReflux().run(nsteps, reflux=False),
            "reflux": TwoLevelReflux().run(nsteps, reflux=True)}


__all__ = ["TwoLevelReflux", "demo"]


if __name__ == "__main__":
    d = demo()
    print("AMR Milestone 1 -- refluxing on a static 2-level hierarchy")
    print("  total-mass relative drift over 40 steps:")
    print("    WITHOUT refluxing: %.3e  (coarse-fine interface leak)" % d["no_reflux"])
    print("    WITH    refluxing: %.3e  (conserved to machine precision)" % d["reflux"])
