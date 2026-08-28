"""Port scaffold -- storm physics on an AMReX (pyAMReX) data model.

This is the first concrete step of the M3 "full AMR" port (`docs/amr_design.md`):
the field lives in an **AMReX `MultiFab`** (the framework's block-structured,
ghost-celled, distributable data model) and **AMReX does the halo exchange**
(`fill_boundary`), while **our physics** (a NumPy flux-form stencil, the same kind
the storm core uses) computes the per-cell update on the MultiFab-backed array.
It proves the binding "AMReX infrastructure + our RHS" that the full port rests
on, and runs on the validated pyAMReX build (`scripts/build_pyamrex_wsl.sh`).

Requires ``pyamrex`` (WSL: ``conda activate amr312``); on a plain CPU/Windows env
it is import-safe but raises a clear message when used.  The two-level, refluxed
version reuses the verified operators in :mod:`storm_dynamics.amr`
(``TwoLevelReflux`` refluxing, conservative restriction/prolongation) as the
coarse-fine coupling once ``MLMG``/``FluxRegister`` bindings are enabled (see the
design doc); this scaffold nails the single-level AMReX<->physics binding first.
"""
from __future__ import annotations

import numpy as np

try:                                        # pyAMReX is WSL-only here
    import amrex.space3d as amr
    _HAVE_AMREX = True
except Exception:                           # pragma: no cover - env dependent
    amr = None
    _HAVE_AMREX = False


def have_amrex() -> bool:
    return _HAVE_AMREX


def _require():
    if not _HAVE_AMREX:
        raise RuntimeError(
            "pyAMReX not available. Build it in WSL with "
            "scripts/build_pyamrex_wsl.sh, then run with the amr312 conda python "
            "(see docs/amr_design.md).")


class AmrexAdvect3D:
    """Single-level, periodic, uniform-velocity flux-form advection with the field
    stored in an AMReX ``MultiFab``.

    AMReX owns the grid (``Geometry``/``BoxArray``) and the ghost exchange; the
    update is our upwind stencil on the ``to_numpy`` view.  Periodic + flux form =>
    the total mass is conserved to machine precision (the correctness check).
    """

    def __init__(self, n=32, L=1.0, U=1.0, V=0.7, W=0.4, cfl=0.4):
        _require()
        self.n, self.L = n, L
        self.U, self.V, self.W = U, V, W
        self.dx = L / n
        self.dt = cfl * self.dx / max(U, V, W)
        lo = amr.IntVect(0, 0, 0)
        hi = amr.IntVect(n - 1, n - 1, n - 1)
        self.domain = amr.Box(lo, hi)
        self.ba = amr.BoxArray(self.domain)          # single box
        self.dm = amr.DistributionMapping(self.ba)
        self.geom = amr.Geometry(self.domain, amr.RealBox(0, 0, 0, L, L, L),
                                 0, [1, 1, 1])        # periodic in x,y,z
        self.mf = amr.MultiFab(self.ba, self.dm, 1, 1)  # 1 comp, 1 ghost
        self.mf.set_val(0.0)
        # initialise a smooth tracer on the valid interior
        a = self._view()
        xs = (np.arange(n) + 0.5) * self.dx
        X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
        a[1:-1, 1:-1, 1:-1] = 1.0 + 0.3 * np.sin(2 * np.pi * X / L) * np.cos(2 * np.pi * Y / L) \
            * np.sin(2 * np.pi * Z / L)

    def _view(self):
        """Live (nx+2, ny+2, nz+2) view of the MultiFab (single box, ghost=1)."""
        return self.mf.to_numpy()[0][..., 0]

    def mass(self) -> float:
        return float(self.mf.sum(0)) * self.dx ** 3

    def step(self) -> None:
        self.mf.fill_boundary(self.geom.periodicity())   # AMReX halo exchange
        a = self._view()
        av = a[1:-1, 1:-1, 1:-1]
        lx = a[0:-2, 1:-1, 1:-1]; ly = a[1:-1, 0:-2, 1:-1]; lz = a[1:-1, 1:-1, 0:-2]
        dt, dx = self.dt, self.dx
        new = av - dt / dx * (self.U * (av - lx) + self.V * (av - ly) + self.W * (av - lz))
        a[1:-1, 1:-1, 1:-1] = new                        # write back into the MultiFab

    def run(self, nsteps=50) -> float:
        m0 = self.mass()
        for _ in range(nsteps):
            self.step()
        return abs(self.mass() - m0) / abs(m0)


def demo(n=32, nsteps=50) -> dict:
    """Initialise/step/finalise AMReX; return the relative mass drift."""
    _require()
    amr.initialize([])
    try:
        drift = AmrexAdvect3D(n=n).run(nsteps)
    finally:
        amr.finalize()
    return {"amrex_version": amr.__version__ if _HAVE_AMREX else None,
            "n": n, "nsteps": nsteps, "mass_rel_drift": drift}


__all__ = ["AmrexAdvect3D", "demo", "have_amrex"]


if __name__ == "__main__":
    d = demo()
    print("AMReX port scaffold -- storm physics on a MultiFab")
    print("  pyAMReX:", d["amrex_version"])
    print("  %d^3 periodic advection, %d steps" % (d["n"], d["nsteps"]))
    print("  total-mass relative drift: %.3e  (periodic + flux form -> ~0)"
          % d["mass_rel_drift"])
