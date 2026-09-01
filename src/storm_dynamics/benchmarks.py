"""Standard idealised benchmarks for quantitative verification (ROADMAP §3e).

Beyond the model's own conservation checks, these are the community reference problems that
turn "indicative" into "quantitative".  First one: the **Straka et al. (1993) density
current** -- a cold bubble collapsing in a neutral, dry, non-rotating atmosphere with a
*fixed* viscosity, so the solution is grid-converged and every code should agree on the front
position and the Kelvin-Helmholtz rotor structure.

Straka, J. M., R. B. Wilhelmson, L. J. Wicker, J. R. Anderson, K. K. Droegemeier (1993),
"Numerical solutions of a non-linear density current: A benchmark solution and comparisons",
Int. J. Numer. Methods Fluids 17, 1-22.  Setup: domain 51.2 km x 6.4 km, neutral theta=300 K,
cold bubble Delta_theta = -15 K * [cos(pi L)+1]/2 (L the elliptic radius about x_c, z=3 km,
half-axes 4 km x 2 km), constant nu = kappa = 75 m^2/s, integrate to 900 s.  Reference front
(leading edge of the surface outflow) reaches ~15.5 km from the centre at the 100 m resolution.
"""
from __future__ import annotations

import numpy as np

from meteorological_flow import thermodynamics as th
from meteorological_flow.base_state import BaseState
from meteorological_flow.grid import Grid

from .config import build_storm_config
from .core import StormSimulation


def neutral_dry_base(grid: Grid, theta_K: float = 300.0, p_sfc_Pa: float = 1.0e5) -> BaseState:
    """A neutral (``dtheta/dz = 0``), dry (``qv = 0``), resting base state -- the Straka
    environment.  Pressure is integrated hydrostatically for the constant-theta column."""
    zc = np.asarray(grid.backend.to_cpu(grid.zc), float)
    nz = zc.size
    theta0 = np.full(nz, float(theta_K))
    qv0 = np.zeros(nz)
    p0 = np.empty(nz); T0 = np.empty(nz)
    p_prev, z_prev = float(p_sfc_Pa), 0.0
    for k in range(nz):
        dz = zc[k] - z_prev; p_new = p_prev
        for _ in range(3):
            Tk = float(th.T_from_theta(theta0[k], p_new, th.P0_REF))
            p_new = p_prev * float(np.exp(-th.g0 * dz / (th.R_d * Tk)))
        p0[k] = p_new; T0[k] = float(th.T_from_theta(theta0[k], p0[k], th.P0_REF))
        p_prev, z_prev = p_new, zc[k]
    rho0 = p0 / (th.R_d * T0)
    return BaseState(zc=zc, theta0=theta0, qv0=qv0, p0=p0, T0=T0, rho0=rho0,
                     u0=np.zeros(nz), v0=np.zeros(nz))


def _straka_bubble(grid, xc, zc_b=3000.0, xr=4000.0, zr=2000.0, dtheta=-15.0):
    """The Straka elliptic cold-bubble potential-temperature perturbation (nx,ny,nz)."""
    xp = grid.xp
    X = grid.xc.reshape(-1, 1, 1); Z = grid.zc.reshape(1, 1, -1)
    L = xp.sqrt(((X - xc) / xr) ** 2 + ((Z - zc_b) / zr) ** 2)
    dth = xp.where(L <= 1.0, dtheta * 0.5 * (xp.cos(np.pi * L) + 1.0), 0.0)
    return xp.broadcast_to(dth, grid.center_shape).copy()


def straka_simulation(nz: int = 32, Lx: float = 51200.0, Lz: float = 6400.0, ny: int = 4,
                      duration: float = 900.0, nu: float = 75.0, device: str = "cpu",
                      dt_max: float = 2.0) -> StormSimulation:
    """Build the Straka density-current :class:`StormSimulation` (``dx = dz``, ``nx = 8 nz``;
    ``nz=32`` -> 200 m coarse test, ``nz=64`` -> 100 m reference).  Neutral/dry base, no
    Coriolis, no drag, constant viscosity ``nu`` (``les_model='none'`` + ``nu_background=nu``,
    ``Pr_t=1`` so ``kappa=nu``), a resting cold bubble.  Effectively 2-D (thin, y-periodic)."""
    nx = 8 * nz                                            # dx == dz
    dx = Lx / nx
    Ly = dx * ny                                          # dy == dx (thin y)
    scfg = build_storm_config(preset="storm", nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz,
                              duration=duration, dt_max=dt_max, coriolis=False, drag=False,
                              les_model="none", z_stretch=1.0, device=device)
    scfg.dyn.les.model = "none"; scfg.dyn.les.nu_background = float(nu); scfg.dyn.les.Pr_t = 1.0
    scfg.sim.physics.bubble_dtheta = 0.0                  # no warm bubble; we inject the cold one
    grid = Grid(nx=nx, ny=ny, nz=nz, Lx=Lx, Ly=Ly, Lz=Lz, z_stretch=1.0, periodic=True)
    base = neutral_dry_base(grid)
    sim = StormSimulation(scfg, base=base)
    sim.state.theta = sim.theta0_field + _straka_bubble(sim.grid, xc=0.5 * Lx)
    sim.state.diagnose(sim.cfg)
    return sim


def straka_front_position(sim: StormSimulation, theta_thresh: float = -1.0) -> float:
    """Distance [m] from the domain centre to the leading edge of the surface cold outflow --
    the furthest surface cell whose potential-temperature perturbation is below
    ``theta_thresh`` [K].  The Straka reference front is ~15.5 km at 900 s."""
    to = sim.grid.backend.to_cpu
    dth_sfc = np.asarray(to(sim.state.theta - sim.theta0_field))[:, 0, 0]   # surface, y=0
    x = np.asarray(to(sim.grid.xc)); xc = 0.5 * float(sim.grid.Lx)
    cold = np.where(dth_sfc < theta_thresh)[0]
    return float(np.max(np.abs(x[cold] - xc))) if cold.size else 0.0
