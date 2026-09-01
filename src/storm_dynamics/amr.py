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


    # ---- free-stream preservation: a uniform field must stay uniform ----
    def free_stream_error(self, nsteps: int = 30, reflux: bool = True) -> float:
        """Set a uniform field and return the max deviation from uniform after
        ``nsteps`` -- an AMR correctness test (a sign error in refluxing or a bad
        coarse-fine transfer shows up as spurious structure at the interface)."""
        self.qc[:] = 2.0
        self.qf[:] = 2.0
        for _ in range(nsteps):
            self.step(reflux=reflux)
        return float(max(np.abs(self.qc - 2.0).max(), np.abs(self.qf - 2.0).max()))


def _burgers_upwind(u, dt, h, periodic=True, ghost=None):
    """Flux-form 2-D inviscid-Burgers step, the **nonlinear** momentum flux ``f(u)=u²/2`` in
    x and y, upwind for ``u > 0``.  Returns ``(u_new, Fx, Fy)`` with the face-flux arrays
    ``Fx`` (nx+1, ny), ``Fy`` (nx, ny+1).  Non-periodic patches take the incoming (left column,
    bottom row) **u-values** from ``ghost`` and square them.  Total momentum ``sum(u)`` is
    conserved (flux-form); unlike a passive scalar under a *uniform* velocity, the interface
    flux ``u²/2`` genuinely differs between the coarse and fine sides of a resolved gradient,
    so refluxing has real work to do — the momentum-conservation test."""
    nx, ny = u.shape
    fc = 0.5 * u * u                                   # cell-centred flux f(u)=u²/2
    Fx = np.zeros((nx + 1, ny)); Fy = np.zeros((nx, ny + 1))
    Fx[1:, :] = fc                                     # face i+1/2 <- upwind cell i (u>0)
    Fx[0, :] = 0.5 * (u[-1, :] ** 2 if periodic else ghost[0] ** 2)
    Fy[:, 1:] = fc
    Fy[:, 0] = 0.5 * (u[:, -1] ** 2 if periodic else ghost[1] ** 2)
    un = u - dt / h * (Fx[1:, :] - Fx[:-1, :]) - dt / h * (Fy[:, 1:] - Fy[:, :-1])
    return un, Fx, Fy


@dataclass
class TwoLevelBurgersReflux:
    """Static 2-level **nonlinear-momentum** advection (inviscid Burgers, flux ``u²/2``) with
    optional Berger–Colella refluxing — the momentum-conservation companion to
    :class:`TwoLevelReflux` (ROADMAP §2b).  Same geometry: a coarse ``N×N`` periodic grid with
    one patch over ``[ci0,ci1)×[cj0,cj1)`` refined by ``refine`` and sub-cycled ``refine``
    times.  Because the flux is nonlinear the coarse interface flux ``u_c²/2`` differs from the
    average of the fine interface fluxes ``⟨u_f²/2⟩`` across a gradient, so **without** reflux
    the total momentum ``∑u`` drifts; **with** reflux the flux mismatch on each coarse–fine
    face corrects the coarse cell just outside the patch and momentum is conserved to machine
    precision.  This is the algorithm the storm's staggered momentum reflux ports onto a
    ``FluxRegister``; here it is nailed down in pure NumPy first."""
    N: int = 24
    dx: float = 1.0
    refine: int = 3
    ci0: int = 8
    ci1: int = 16
    cj0: int = 8
    cj1: int = 16
    cfl: float = 0.4
    uc: np.ndarray = field(default=None, repr=False)
    uf: np.ndarray = field(default=None, repr=False)

    def __post_init__(self):
        r = self.refine
        i = np.arange(self.N)[:, None]; j = np.arange(self.N)[None, :]
        # strictly positive momentum (upwind for u>0), with a resolved gradient
        self.uc = 1.0 + 0.3 * np.sin(2 * np.pi * i / self.N) * np.cos(2 * np.pi * j / self.N)
        self.uf = np.repeat(np.repeat(self.uc[self.ci0:self.ci1, self.cj0:self.cj1],
                                      r, axis=0), r, axis=1).copy()
        self.dxf = self.dx / r
        self.dtc = self.cfl * self.dx / float(self.uc.max())

    def total_momentum(self) -> float:
        mask = np.ones_like(self.uc, bool)
        mask[self.ci0:self.ci1, self.cj0:self.cj1] = False
        return float(self.uc[mask].sum() * self.dx ** 2 + self.uf.sum() * self.dxf ** 2)

    def step(self, reflux: bool = True) -> None:
        r, dx, dxf = self.refine, self.dx, self.dxf
        ci0, ci1, cj0, cj1 = self.ci0, self.ci1, self.cj0, self.cj1
        uc_new, Fxc, Fyc = _burgers_upwind(self.uc, self.dtc, dx, periodic=True)
        reg_L = Fxc[ci0, cj0:cj1] * self.dtc * dx        # coarse interface momentum fluxes
        reg_R = Fxc[ci1, cj0:cj1] * self.dtc * dx
        reg_B = Fyc[ci0:ci1, cj0] * self.dtc * dx
        reg_T = Fyc[ci0:ci1, cj1] * self.dtc * dx
        dtf = self.dtc / r
        fL = np.zeros(ci1 - ci0); fR = np.zeros(ci1 - ci0)
        fB = np.zeros(ci1 - ci0); fT = np.zeros(ci1 - ci0)
        for _ in range(r):
            gl = np.repeat(self.uc[ci0 - 1, cj0:cj1], r)  # ghost u-values (squared in the flux)
            gb = np.repeat(self.uc[ci0:ci1, cj0 - 1], r)
            self.uf, Fxf, Fyf = _burgers_upwind(self.uf, dtf, dxf, periodic=False, ghost=(gl, gb))
            fL += Fxf[0, :].reshape(ci1 - ci0, r).sum(1) * dtf * dxf
            fR += Fxf[-1, :].reshape(ci1 - ci0, r).sum(1) * dtf * dxf
            fB += Fyf[:, 0].reshape(ci1 - ci0, r).sum(1) * dtf * dxf
            fT += Fyf[:, -1].reshape(ci1 - ci0, r).sum(1) * dtf * dxf
        self.uc = uc_new
        self.uc[ci0:ci1, cj0:cj1] = self.uf.reshape(ci1 - ci0, r, cj1 - cj0, r).mean(axis=(1, 3))
        if reflux:
            Vc = dx ** 2
            self.uc[ci0 - 1, cj0:cj1] += (reg_L - fL) / Vc
            self.uc[ci1, cj0:cj1] -= (reg_R - fR) / Vc
            self.uc[ci0:ci1, cj0 - 1] += (reg_B - fB) / Vc
            self.uc[ci0:ci1, cj1] -= (reg_T - fT) / Vc

    def run(self, nsteps: int = 40, reflux: bool = True) -> float:
        """Run ``nsteps`` and return the relative total-momentum drift ``|Δ∑u|/∑u₀``."""
        m0 = self.total_momentum()
        for _ in range(nsteps):
            self.step(reflux=reflux)
        return abs(self.total_momentum() - m0) / abs(m0)

    def free_stream_error(self, nsteps: int = 30, reflux: bool = True) -> float:
        """A uniform momentum field (uniform flux -> zero divergence) must stay uniform; a
        reflux sign error shows up as spurious interface structure."""
        self.uc[:] = 2.0; self.uf[:] = 2.0
        for _ in range(nsteps):
            self.step(reflux=reflux)
        return float(max(np.abs(self.uc - 2.0).max(), np.abs(self.uf - 2.0).max()))


def conservative_prolong(coarse_block: np.ndarray, refine: int) -> np.ndarray:
    """Coarse -> fine conservative interpolation (the regridding operator, the
    companion to :func:`storm_dynamics.nesting.conservative_restrict`).

    Piecewise-constant injection: each coarse cell fills its ``refine x refine``
    fine cells.  This is exactly conservative (the fine integral equals the coarse
    integral) and is the inverse of average-down on a constant block:
    ``restrict(prolong(x)) == x``.
    """
    r = refine
    return np.repeat(np.repeat(coarse_block, r, axis=0), r, axis=1)


def demo(nsteps: int = 40) -> dict:
    """Relative mass drift with and without refluxing, and the free-stream error."""
    return {"no_reflux": TwoLevelReflux().run(nsteps, reflux=False),
            "reflux": TwoLevelReflux().run(nsteps, reflux=True),
            "free_stream": TwoLevelReflux().free_stream_error(nsteps, reflux=True),
            "momentum_no_reflux": TwoLevelBurgersReflux().run(nsteps, reflux=False),
            "momentum_reflux": TwoLevelBurgersReflux().run(nsteps, reflux=True),
            "momentum_free_stream": TwoLevelBurgersReflux().free_stream_error(nsteps, reflux=True)}


__all__ = ["TwoLevelReflux", "TwoLevelBurgersReflux", "conservative_prolong", "demo"]


if __name__ == "__main__":
    d = demo()
    print("AMR Milestone 1 -- refluxing on a static 2-level hierarchy")
    print("  total-mass relative drift over 40 steps:")
    print("    WITHOUT refluxing: %.3e  (coarse-fine interface leak)" % d["no_reflux"])
    print("    WITH    refluxing: %.3e  (conserved to machine precision)" % d["reflux"])
