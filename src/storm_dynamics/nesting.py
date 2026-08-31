"""M3 (phase 1) -- static one-way nested-grid refinement.

Full adaptive mesh refinement (AMR) is a separate project; this module delivers
the classical idealised-tornado approach instead: mature the storm on the coarse
**parent** domain, then integrate a finer **nest** over the low-level vortex
region, with the nest's lateral boundaries relaxed toward the parent solution
(Davies-style nudging).  This resolves the near-surface vortex at O(100 m) in the
region of interest without refining the whole domain.

Scope / honesty (see docs/storm_dynamics_guide.md):

* **one-way** -- the parent drives the nest; the nest does not feed back;
* **static** -- a fixed nest region and refinement (not adaptive);
* **frozen-parent boundary** -- the nest border is nudged toward the parent state
  captured at the nest start time, so the nest is valid for a **short window**
  (minutes) after which the parent boundary would have evolved away;
* still **idealised** and, at O(100 m), only *approaching* a resolved vortex --
  not a converged tornado, and never a forecast.

The nest reuses the whole :class:`~storm_dynamics.core.StormSimulation` stepping
machinery (momentum advection, LES, surface drag, projection, microphysics,
rotation diagnostics); only the initial state (interpolated from the parent), the
base state (interpolated in z), and the boundary treatment (walls + relaxation
nudging, instead of periodic) differ.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from meteorological_flow.base_state import BaseState
from meteorological_flow.grid import Grid
from meteorological_flow.state import FlowState


@dataclass
class NestSpec:
    """A static nest: a finer sub-domain of the parent over the vortex region.

    ``x0, y0`` are the parent-relative offsets [m] of the nest's south-west
    corner; ``Lx, Ly`` its horizontal extent [m].  The nest spans the **full**
    parent depth (so the deep updraft is kept) at ``refine`` times finer
    horizontal spacing, and ``nz`` vertical levels with optional near-surface
    clustering ``z_stretch``.  ``relax_width`` cells of boundary nudging at
    ``relax_rate`` [1/s] tie the border to the parent.
    """
    x0: float
    y0: float
    Lx: float
    Ly: float
    refine: int = 3                 # horizontal refinement factor
    nz: int | None = None           # nest vertical cells (default: parent nz)
    z_stretch: float = 1.03         # cluster levels near the surface
    relax_width: int = 4            # boundary relaxation band (cells)
    relax_rate: float = 0.02        # max nudging rate at the outermost cell [1/s]

    @classmethod
    def centered(cls, parent_grid: Grid, frac: float = 0.4, **kw) -> "NestSpec":
        """A nest covering the central ``frac`` of the parent (x, y)."""
        Lx = frac * parent_grid.Lx
        Ly = frac * parent_grid.Ly
        x0 = 0.5 * (parent_grid.Lx - Lx)
        y0 = 0.5 * (parent_grid.Ly - Ly)
        return cls(x0=x0, y0=y0, Lx=Lx, Ly=Ly, **kw)

    @classmethod
    def around(cls, parent_grid: Grid, xc: float, yc: float, half: float, **kw) -> "NestSpec":
        """A square nest of half-width ``half`` [m] centred on (xc, yc), clipped
        to the parent domain."""
        x0 = max(0.0, min(xc - half, parent_grid.Lx - 2 * half))
        y0 = max(0.0, min(yc - half, parent_grid.Ly - 2 * half))
        return cls(x0=x0, y0=y0, Lx=2 * half, Ly=2 * half, **kw)

    @classmethod
    def aligned(cls, parent_grid: Grid, i0: int, j0: int, ncx: int, ncy: int,
                refine: int = 3, **kw) -> "NestSpec":
        """A nest **cell-aligned** to the parent: its SW corner sits on parent cell
        ``(i0, j0)`` and it covers ``ncx x ncy`` parent cells, refined ``refine``.
        With the parent's vertical grid (matched z, set here) the fine cells nest
        *exactly* into the coarse ones, so :func:`conservative_restrict` (average
        down) preserves the overlap integral to machine precision."""
        kw.setdefault("nz", parent_grid.nz)
        kw.setdefault("z_stretch", parent_grid.z_stretch)
        return cls(x0=i0 * parent_grid.dx, y0=j0 * parent_grid.dy,
                   Lx=ncx * parent_grid.dx, Ly=ncy * parent_grid.dy, refine=refine, **kw)


def build_nest_grid(spec: NestSpec, parent_grid: Grid, backend=None) -> Grid:
    """Build the (non-periodic) finer nest :class:`Grid` for ``spec``."""
    dxp = parent_grid.dx
    nx = max(4, int(round(spec.Lx / dxp * spec.refine)))
    ny = max(4, int(round(spec.Ly / parent_grid.dy * spec.refine)))
    nz = spec.nz or parent_grid.nz
    return Grid(nx=nx, ny=ny, nz=nz, Lx=spec.Lx, Ly=spec.Ly, Lz=parent_grid.Lz,
                z_stretch=spec.z_stretch, periodic=False,
                backend=backend or parent_grid.backend)


# ---------------------------------------------------------------------------
# interpolation parent -> nest  (host/NumPy; one-time at nest init)
# ---------------------------------------------------------------------------
def _axis(parent_grid, kind):
    """Parent coordinate axes (host) for a field of the given staggered kind."""
    to = parent_grid.backend.to_cpu
    xc = np.asarray(to(parent_grid.xc)); yc = np.asarray(to(parent_grid.yc))
    zc = np.asarray(to(parent_grid.zc)); xf = np.asarray(to(parent_grid.xf))
    yf = np.asarray(to(parent_grid.yf)); zf = np.asarray(to(parent_grid.zf))
    return {"c": (xc, yc, zc), "u": (xf, yc, zc),
            "v": (xc, yf, zc), "w": (xc, yc, zf)}[kind]


def _target(spec, nest_grid, kind):
    """Nest sample points in PARENT physical coordinates for a field kind."""
    to = nest_grid.backend.to_cpu
    xc = spec.x0 + np.asarray(to(nest_grid.xc)); yc = spec.y0 + np.asarray(to(nest_grid.yc))
    zc = np.asarray(to(nest_grid.zc)); xf = spec.x0 + np.asarray(to(nest_grid.xf))
    yf = spec.y0 + np.asarray(to(nest_grid.yf)); zf = np.asarray(to(nest_grid.zf))
    return {"c": (xc, yc, zc), "u": (xf, yc, zc),
            "v": (xc, yf, zc), "w": (xc, yc, zf)}[kind]


def _interp(parent_grid, field_host, kind, spec, nest_grid):
    ax = _axis(parent_grid, kind)
    tgt = _target(spec, nest_grid, kind)
    rgi = RegularGridInterpolator(ax, field_host, method="linear",
                                  bounds_error=False, fill_value=None)
    # clamp targets to the parent range (avoid wild linear extrapolation past edges)
    gx = np.clip(tgt[0], ax[0][0], ax[0][-1])
    gy = np.clip(tgt[1], ax[1][0], ax[1][-1])
    gz = np.clip(tgt[2], ax[2][0], ax[2][-1])
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=-1)
    return rgi(pts).reshape(X.shape)


def interpolate_state_to_nest(parent_state: FlowState, parent_grid: Grid,
                              nest_grid: Grid, spec: NestSpec) -> FlowState:
    """Trilinearly interpolate the parent prognostic state onto the nest grid."""
    xp = nest_grid.xp
    to = parent_grid.backend.to_cpu
    st = FlowState.zeros(nest_grid)
    st.u = xp.asarray(_interp(parent_grid, np.asarray(to(parent_state.u)), "u", spec, nest_grid))
    st.v = xp.asarray(_interp(parent_grid, np.asarray(to(parent_state.v)), "v", spec, nest_grid))
    st.w = xp.asarray(_interp(parent_grid, np.asarray(to(parent_state.w)), "w", spec, nest_grid))
    for nm in ("theta", "qv", "ql", "qi", "qr", "qs", "qg", "qh"):
        a = getattr(parent_state, nm, None)
        if a is not None:
            setattr(st, nm, xp.asarray(np.maximum(
                _interp(parent_grid, np.asarray(to(a)), "c", spec, nest_grid),
                0.0 if nm != "theta" else -1e30)))
    st.w[:, :, 0] = 0.0                          # ground
    return st


def interpolate_base_to_nest(parent_base: BaseState, nest_grid: Grid) -> BaseState:
    """Interpolate the parent base-state (z-only) profiles onto the nest levels."""
    zc_p = np.asarray(parent_base.zc)
    zc_n = np.asarray(nest_grid.backend.to_cpu(nest_grid.zc))
    f = lambda a: np.interp(zc_n, zc_p, np.asarray(a))
    return BaseState(zc=zc_n, theta0=f(parent_base.theta0), qv0=f(parent_base.qv0),
                     p0=f(parent_base.p0), T0=f(parent_base.T0), rho0=f(parent_base.rho0),
                     u0=f(parent_base.u0), v0=f(parent_base.v0))


def relaxation_weight(nest_grid: Grid, spec: NestSpec):
    """A (nx,ny,1) nudging-weight field: ``relax_rate`` at the outermost cell,
    ramping quadratically to 0 at ``relax_width`` cells in; 0 in the interior."""
    xp = nest_grid.xp
    nx, ny = nest_grid.nx, nest_grid.ny
    ix = xp.arange(nx); iy = xp.arange(ny)
    dx_edge = xp.minimum(ix, nx - 1 - ix)          # (nx,) distance to nearest x-edge
    dy_edge = xp.minimum(iy, ny - 1 - iy)          # (ny,)
    dist = xp.minimum(dx_edge[:, None], dy_edge[None, :])   # (nx,ny)
    w = xp.clip(1.0 - dist / max(spec.relax_width, 1), 0.0, 1.0) ** 2
    return (spec.relax_rate * w)[:, :, None]


# ---------------------------------------------------------------------------
# the nested simulation
# ---------------------------------------------------------------------------
import copy as _copy
import dataclasses as _dc
import time as _time

from meteorological_flow import advection as _adv
from meteorological_flow import buoyancy as _buo
from meteorological_flow import diagnostics as _diag
from meteorological_flow import thermodynamics as _th
from meteorological_flow.pressure_solver import PressureSolver

from . import momentum as _mom
from . import rotation as _rot
from . import surface_drag as _sfc
from . import turbulence as _les


def _zero_gradient_velocity(st):
    """Outflow (zero normal gradient) lateral velocity + solid ground/top-w=0."""
    st.u[0, :, :] = st.u[1, :, :]; st.u[-1, :, :] = st.u[-2, :, :]
    st.v[:, 0, :] = st.v[:, 1, :]; st.v[:, -1, :] = st.v[:, -2, :]
    st.w[:, :, 0] = 0.0; st.w[:, :, -1] = 0.0


def _zero_gradient_scalars(st):
    for nm in ("theta", "qv", "ql", "qi", "qr", "qs", "qg", "qh"):
        a = getattr(st, nm, None)
        if a is None:
            continue
        a[0, :, :] = a[1, :, :]; a[-1, :, :] = a[-2, :, :]
        a[:, 0, :] = a[:, 1, :]; a[:, -1, :] = a[:, -2, :]


class NestedStormSimulation:
    """One-way static nest: a finer StormSimulation over the vortex region, with
    its border relaxed toward the (frozen) interpolated parent state.

    Build from a matured parent :class:`~storm_dynamics.core.StormSimulation` and
    a :class:`NestSpec`; then :meth:`run` for a short window.
    """

    def __init__(self, parent, spec: NestSpec, les_boost: float = 1.25,
                 cfl: float = 0.25):
        from .core import _pressure_method
        self.spec = spec
        self.scfg = parent.scfg
        # deep-copy dyn so nest tweaks never mutate the parent; the finer,
        # potentially under-resolved nest gets a little more SGS dissipation
        # (les_boost x C_s) and a tighter CFL for stability as its vortex sharpens.
        self.dyn = _copy.deepcopy(parent.dyn)
        self.dyn.les.C_s = parent.dyn.les.C_s * float(les_boost)
        self.backend = parent.backend
        # nest config: reuse the parent's, but non-periodic walls + nest geometry
        cfg = _copy.deepcopy(parent.cfg)
        self.grid = build_nest_grid(spec, parent.grid, backend=self.backend)
        g = self.grid
        xp = g.xp
        cfg.grid.nx, cfg.grid.ny, cfg.grid.nz = g.nx, g.ny, g.nz
        cfg.domain.Lx, cfg.domain.Ly, cfg.domain.Lz = g.Lx, g.Ly, g.Lz
        cfg.boundaries.x_west = cfg.boundaries.x_east = "wall"
        cfg.boundaries.y = "free_slip"
        # finer dx -> tighter CFL + smaller dt cap for the fine vertical spacing
        cfg.time.cfl = min(cfg.time.cfl, float(cfl))
        cfg.time.dt_max = min(cfg.time.dt_max, 2.0)
        self.cfg = cfg
        # interpolated base state + initial condition
        self.base = interpolate_base_to_nest(parent.base, g)
        self.theta0_field = self.base.field(self.base.theta0, g.center_shape, xp=xp)
        self.qv0_field = self.base.field(self.base.qv0, g.center_shape, xp=xp)
        self.f = self.dyn.f_value()
        self._u0_face = xp.broadcast_to(xp.asarray(self.base.u0, float)[None, None, :],
                                        g.u_shape).copy()
        self._v0_face = xp.broadcast_to(xp.asarray(self.base.v0, float)[None, None, :],
                                        g.v_shape).copy()
        self.state = interpolate_state_to_nest(parent.state, parent.grid, g, spec)
        self.state.p0_field = self.base.field(self.base.p0, g.center_shape, xp=xp)
        self.state.diagnose(cfg)
        # frozen relaxation target (the initial interpolated parent state) + weight
        self._tgt = {nm: getattr(self.state, nm).copy()
                     for nm in ("u", "v", "w", "theta", "qv")}
        self._relax = relaxation_weight(g, spec)
        # references / density (anelastic reference from the interpolated base)
        self.T_ref = parent.T_ref
        self.qv_ref = parent.qv_ref
        self.rho0 = parent.rho0
        self.dynamics = parent.dynamics
        rho0_prof = np.asarray(self.base.rho0, float)
        rho0_wface_h = np.interp(g.backend.to_cpu(g.zf), g.backend.to_cpu(g.zc), rho0_prof)
        self.rho0_c = xp.asarray(rho0_prof); self.rho0_wface = xp.asarray(rho0_wface_h)
        if self.dynamics == "anelastic":
            self._transport_rho_c = self.rho0_c; self._transport_rho_wf = self.rho0_wface
            self.rho_ref = self.rho0_c
        else:
            self._transport_rho_c = xp.ones(g.nz); self._transport_rho_wf = xp.ones(g.nz + 1)
            self.rho_ref = self.rho0
        self.pressure = PressureSolver(g, method=_pressure_method(g))
        from meteorological_flow.microphysics_coupling import MicrophysicsCoupler
        self.coupler = MicrophysicsCoupler()
        self.couple_nucleation = False
        self.tracker = _rot.TornadogenesisTracker()
        self.history = []
        self.step = 0; self.t = 0.0
        self._last_res = 0.0; self._last_iters = 0
        self._Km = None
        self._composite = False          # composite_projection mode replaces the sponge
        self._t0 = _time.perf_counter()

    # ---- boundary / relaxation ----
    def set_target(self, tgt: dict) -> None:
        """Set the relaxation target (parent state interpolated onto the nest).
        For phase 1 this is the frozen initial state; for phase 2
        (:func:`run_concurrent_nest`) it is updated each nest step."""
        self._tgt = tgt

    def _apply_nest_bcs(self):
        _zero_gradient_velocity(self.state)
        _zero_gradient_scalars(self.state)

    def _relax_to_parent(self, dt):
        """Nudge the border band toward the frozen parent target (sponge).  In
        composite-projection mode the velocity sponge is DROPPED -- the composite
        solve at the coarse-fine interface replaces it (docs/ROADMAP.md §1); only
        the scalar border still relaxes (open-wall BC needs a scalar tie)."""
        st = self.state
        wu = self._relax                     # (nx,ny,1) centre-cell nudging weight
        xp = self.grid.xp
        if not getattr(self, "_composite", False):
            # velocity relaxation weights averaged onto the staggered faces
            wru = xp.zeros(self.grid.u_shape); wru[1:-1] = 0.5 * (wu[:-1] + wu[1:])
            wru[0] = wu[0]; wru[-1] = wu[-1]
            wrv = xp.zeros(self.grid.v_shape); wrv[:, 1:-1] = 0.5 * (wu[:, :-1] + wu[:, 1:])
            wrv[:, 0] = wu[:, 0]; wrv[:, -1] = wu[:, -1]
            st.u += wru * dt * (self._tgt["u"] - st.u)
            st.v += wrv * dt * (self._tgt["v"] - st.v)
            st.w += wu * dt * (self._tgt["w"] - st.w)
        st.theta += wu * dt * (self._tgt["theta"] - st.theta)
        st.qv = xp.maximum(st.qv + wu * dt * (self._tgt["qv"] - st.qv), 0.0)

    # ---- adaptive dt (advective + LES-diffusive) ----
    def _dt(self):
        from .core import StormSimulation
        return StormSimulation._dt(self)

    # ---- one step ----
    # Split into phases (mirrors core.StormSimulation) so the composite
    # parent+nest projection can replace this level's own projection.
    def _step(self, dt):
        Km = self._predictor(dt)
        self._project(dt)
        self._transport(dt, Km)

    def _predictor(self, dt):
        cfg = self.cfg; g = self.grid; xp = g.xp; st = self.state
        self._apply_nest_bcs(); st.diagnose(cfg)
        Km = _les.strain_and_viscosity(st, g, self.dyn.les, theta0=self.theta0_field)
        self._Km = Km
        st.u -= self._u0_face; st.v -= self._v0_face
        _les.apply_les_momentum(st, g, Km, dt)
        st.u += self._u0_face; st.v += self._v0_face
        if self.dyn.momentum_advection:
            _mom.add_momentum_advection(st, g, dt, order=self.dyn.momentum_order, periodic=False)
        Bf = _buo.buoyancy_w_tendency(st, g, cfg, self.T_ref, self.qv_ref,
                                      theta0=self.theta0_field, qv0=self.qv0_field)
        st.w += dt * Bf
        if self.dyn.coriolis:
            from .coriolis import add_coriolis
            add_coriolis(st, g, dt, self.f, self._u0_face, self._v0_face)
        _sfc.apply_surface_drag(st, g, dt, self.dyn.drag)
        vg = self.dyn.v_guard
        xp.clip(st.u, -vg, vg, out=st.u); xp.clip(st.v, -vg, vg, out=st.v); xp.clip(st.w, -vg, vg, out=st.w)
        self._apply_nest_bcs()
        return Km

    def _project(self, dt):
        "This level's own anelastic projection (skipped in composite mode)."
        st = self.state
        if self.dynamics == "anelastic":
            res, it = self.pressure.project_anelastic(st, dt, self.rho0_c, self.rho0_wface)
        else:
            res, it = self.pressure.project(st, dt, self.rho0)
        self._last_res, self._last_iters = res, it
        self._apply_nest_bcs()

    def _transport(self, dt, Km):
        cfg = self.cfg; g = self.grid; xp = g.xp; st = self.state
        if not getattr(self, "_composite", False):
            # In composite mode the boundary faces were just written by the
            # composite solve -- re-applying the zero-gradient BC here would
            # clobber the interface coupling (normal modes: harmless).
            self._apply_nest_bcs()
        trc, trwf = self._transport_rho_c, self._transport_rho_wf
        adv = lambda f: _adv.advect_center_massflux(f, st.u, st.v, st.w, g, dt, trc, trwf, order=2)
        adv = lambda f: _adv.advect_center_massflux(f, st.u, st.v, st.w, g, dt, trc, trwf, order=2)
        st.theta = adv(st.theta); st.qv = xp.maximum(adv(st.qv), 0.0)
        st.ensure_hydrometeors()
        for nm in ("ql", "qi", "qr", "qs", "qg", "qh"):
            setattr(st, nm, xp.maximum(adv(getattr(st, nm)), 0.0))
        st.theta = _les.les_scalar_diffusion(st.theta, Km, g, self.dyn.les, dt, base=self.theta0_field)
        st.qv = xp.maximum(_les.les_scalar_diffusion(st.qv, Km, g, self.dyn.les, dt, base=self.qv0_field), 0.0)
        self._relax_to_parent(dt)                    # sponge border -> frozen parent
        self._apply_nest_bcs(); st.diagnose(cfg)
        self.coupler.apply(st, g, dt, nf=None)
        self.coupler.sediment(st, g, dt, rho_ref=self._transport_rho_c)
        self._apply_nest_bcs(); st.diagnose(cfg)
        Tc = xp.clip(st.T, 180.0, 335.0)
        if not bool(xp.array_equal(Tc, st.T)):
            st.theta = _th.theta_from_T(Tc, st.P_total, _th.P0_REF, xp=xp); st.diagnose(cfg)
        st.t = self.t + dt

    def run(self, progress=None, record_interval=None, capture_frames=False) -> dict:
        from .core import StormSimulation
        self._capture_frames = bool(capture_frames); self.frames = []
        initial = _diag.initial_budgets(self.state, self.rho_ref)
        duration = self.cfg.time.duration
        interval = record_interval or max(1, self.cfg.output.interval_steps)
        self.tracker.update(self.t, self.state, self.grid)
        StormSimulation._record(self, initial)
        while self.t < duration - 1e-9:
            dt = self._dt()
            if self.t + dt > duration:
                dt = duration - self.t
            self._step(dt); self.step += 1; self.t = float(self.state.t)
            if self.step % interval == 0 or self.t >= duration - 1e-9:
                self.tracker.update(self.t, self.state, self.grid)
                StormSimulation._record(self, initial)
            if progress and (self.step % max(1, min(interval, 10)) == 0):
                progress(self.t, duration, self.step)
        rep = StormSimulation._finalise(self, initial)
        rep["nest"] = {"dx_m": self.grid.dx, "dz0_m": float(self.grid.dz_c[0]),
                       "nx": self.grid.nx, "ny": self.grid.ny, "nz": self.grid.nz,
                       "region": {"x0": self.spec.x0, "y0": self.spec.y0,
                                  "Lx": self.spec.Lx, "Ly": self.spec.Ly},
                       "refine": self.spec.refine}
        return rep


# reuse the parent class's 2-D frame capture (generic over self.state / self.grid)
from .core import StormSimulation as _SS
NestedStormSimulation._capture_frame = _SS._capture_frame


# ---------------------------------------------------------------------------
# M3 phase 2 -- concurrent one-way nesting (time-evolving parent boundaries)
# ---------------------------------------------------------------------------
def _interp_targets(parent_state: FlowState, parent_grid: Grid,
                    nest_grid: Grid, spec: NestSpec) -> dict:
    """Interpolate just the relaxation fields (u, v, w, theta, qv) onto the nest."""
    xp = nest_grid.xp
    to = parent_grid.backend.to_cpu
    out = {}
    for nm, kind in (("u", "u"), ("v", "v"), ("w", "w"),
                     ("theta", "c"), ("qv", "c")):
        out[nm] = xp.asarray(_interp(parent_grid, np.asarray(to(getattr(parent_state, nm))),
                                     kind, spec, nest_grid))
    return out


def _blend(a: dict, b: dict, alpha: float) -> dict:
    return {k: (1.0 - alpha) * a[k] + alpha * b[k] for k in a}


def interior_near_surface_zeta(nest, margin: int | None = None,
                               z_near: float = 500.0) -> float:
    """Peak |near-surface ζ| in the nest INTERIOR, excluding the boundary
    relaxation band (the sponge, where nudging can generate edge vorticity that
    is not a physical vortex).  This is the honest low-level-rotation measure for
    a nest."""
    g = nest.grid
    m = margin if margin is not None else (nest.spec.relax_width + 2)
    _, _, zeta = _rot.vorticity_3d(nest.state, g)
    zl = g.backend.to_cpu(zeta)[:, :, int(np.argmin(np.abs(np.asarray(
        g.backend.to_cpu(g.zc)) - z_near)))]
    if 2 * m >= min(g.nx, g.ny):
        return float(np.max(np.abs(zl)))
    return float(np.max(np.abs(zl[m:-m, m:-m])))


def _interp_grid(src_grid, field_host, kind, tx, ty, tz):
    """Interpolate a source-grid field to arbitrary target points (host)."""
    to = src_grid.backend.to_cpu
    xc = np.asarray(to(src_grid.xc)); yc = np.asarray(to(src_grid.yc)); zc = np.asarray(to(src_grid.zc))
    xf = np.asarray(to(src_grid.xf)); yf = np.asarray(to(src_grid.yf)); zf = np.asarray(to(src_grid.zf))
    ax = {"c": (xc, yc, zc), "u": (xf, yc, zc), "v": (xc, yf, zc), "w": (xc, yc, zf)}[kind]
    rgi = RegularGridInterpolator(ax, field_host, method="linear",
                                  bounds_error=False, fill_value=None)
    gx = np.clip(tx, ax[0][0], ax[0][-1]); gy = np.clip(ty, ax[1][0], ax[1][-1])
    gz = np.clip(tz, ax[2][0], ax[2][-1])
    X, Y, Z = np.meshgrid(gx, gy, gz, indexing="ij")
    return rgi(np.stack([X.ravel(), Y.ravel(), Z.ravel()], -1)).reshape(X.shape)


def _overlap_taper(coord, lo, hi, margin):
    """1-D weight: 1 inside [lo+margin, hi-margin], ramping to 0 at [lo, hi], 0 outside."""
    c = np.asarray(coord, dtype=float)
    d = np.minimum(c - lo, hi - c)                 # distance into the overlap
    return np.clip(d / max(margin, 1e-9), 0.0, 1.0) * ((c >= lo) & (c <= hi))


def restrict_nest_to_parent(nest, parent, spec: NestSpec, cx: float, cy: float,
                            t: float, rate: float = 0.5, margin_frac: float = 0.25) -> None:
    """M3 phase 3a -- **approximate two-way** feedback: blend the nest's finer
    solution back onto the parent cells it overlaps (converted to the ground
    frame), with a taper that fades to 0 at the overlap edge (no discontinuity).

    This is *injection* feedback (sample fine at coarse points), NOT rigorous
    flux-conservative refluxing -- it demonstrates the nest→parent coupling; strict
    interface conservation (Berger-Colella refluxing + multilevel Poisson) is the
    full-AMR project.
    """
    pg = parent.grid; xp = pg.xp; to = pg.backend.to_cpu
    x0 = spec.x0 + cx * t; y0 = spec.y0 + cy * t
    Lx, Ly = spec.Lx, spec.Ly
    mgx = margin_frac * Lx; mgy = margin_frac * Ly
    pxc = np.asarray(to(pg.xc)); pyc = np.asarray(to(pg.yc)); pzc = np.asarray(to(pg.zc))
    pxf = np.asarray(to(pg.xf)); pyf = np.asarray(to(pg.yf)); pzf = np.asarray(to(pg.zf))

    def weight(px, py):
        wx = _overlap_taper(px, x0, x0 + Lx, mgx)
        wy = _overlap_taper(py, y0, y0 + Ly, mgy)
        return rate * (wx[:, None] * wy[None, :])[:, :, None]

    def blend(name, kind, px, py, pz, frame=0.0):
        fld = getattr(nest.state, name)
        vals = _interp_grid(nest.grid, np.asarray(nest.grid.backend.to_cpu(fld)), kind,
                            px - x0, py - y0, pz) + frame     # -> ground frame
        w = xp.asarray(weight(px, py))
        cur = getattr(parent.state, name)
        setattr(parent.state, name, cur + w * (xp.asarray(vals) - cur))

    blend("u", "u", pxf, pyc, pzc, frame=cx)
    blend("v", "v", pxc, pyf, pzc, frame=cy)
    blend("w", "w", pxc, pyc, pzf)
    for nm in ("theta", "qv", "ql", "qi", "qr", "qs", "qg", "qh"):
        if getattr(parent.state, nm, None) is not None and getattr(nest.state, nm, None) is not None:
            blend(nm, "c", pxc, pyc, pzc)
    for nm in ("qv", "ql", "qi", "qr", "qs", "qg", "qh"):
        a = getattr(parent.state, nm, None)
        if a is not None:
            setattr(parent.state, nm, xp.maximum(a, 0.0))


def conservative_restrict(nest, parent, spec: NestSpec) -> dict:
    """**Conservative** average-down of the nest scalars onto the parent overlap
    (the first rigorous piece of the AMR conservation machinery, distinct from the
    phase-3a *injection* feedback).

    For a **cell-aligned, matched-z** nest (:meth:`NestSpec.aligned`), each coarse
    cell contains exactly ``refine x refine`` fine cells in the same z-column, so
    replacing the coarse value by the fine block-mean preserves the scalar integral
    over the overlap **exactly** (to machine precision).  Returns the per-field
    overlap-integral change (should be ~0).

    NOTE: this is conservative *restriction* (average down of the covered region).
    Full flux-conservative **refluxing** at the coarse-fine interface, and a
    multilevel Poisson solve, are the remaining AMR pieces (see docs/amr_design.md).
    """
    pg = parent.grid; ng = nest.grid; xp = pg.xp
    r = spec.refine
    i0 = int(round(spec.x0 / pg.dx)); j0 = int(round(spec.y0 / pg.dy))
    ncx, ncy = ng.nx // r, ng.ny // r
    if ng.nx != ncx * r or ng.ny != ncy * r or ng.nz != pg.nz:
        raise ValueError("conservative_restrict needs a cell-aligned, matched-z nest "
                         "(build the spec with NestSpec.aligned)")
    dzc = pg.dz_c[None, None, :]
    cell_h = pg.dx * pg.dy
    out = {}
    for nm in ("theta", "qv", "ql", "qi", "qr", "qs", "qg", "qh"):
        cf = getattr(nest.state, nm, None); cp = getattr(parent.state, nm, None)
        if cf is None or cp is None:
            continue
        coarse = cf.reshape(ncx, r, ncy, r, ng.nz).mean(axis=(1, 3))   # exact average down
        before = float(xp.sum(cp[i0:i0 + ncx, j0:j0 + ncy, :] * dzc) * cell_h)
        cp[i0:i0 + ncx, j0:j0 + ncy, :] = coarse
        after = float(xp.sum(cp[i0:i0 + ncx, j0:j0 + ncy, :] * dzc) * cell_h)
        fine_int = float(xp.sum(cf * dzc) * (cell_h / r ** 2))
        out[nm] = {"coarse_after_minus_fine": after - fine_int, "coarse_change": after - before}
    return out


def composite_project_two_level(parent, nest, spec: NestSpec, periodic_h=None) -> dict:
    """**The composite two-level anelastic pressure projection** -- one solve over the
    parent (coarse) + nest (fine) mass fluxes so that ``div(rho0 u) = 0`` holds
    *consistently across the coarse-fine interface*, replacing the two independent
    per-level :meth:`PressureSolver.project_anelastic` calls.  This is the AMR-projection
    call site (:mod:`storm_dynamics.composite_poisson`, ``docs/amr_design.md`` Milestone 2).

    Operates **in place** on ``parent.state`` and ``nest.state`` (call it after both
    levels' momentum predictors, in place of the separate projections).  Requires a
    **cell-aligned, matched-z** nest (:meth:`NestSpec.aligned`); the parent may be
    **rectangular / anisotropic** (``nx != ny``, ``dx != dy``).  ``periodic_h`` overrides
    the horizontal BC (default: the parent grid's).  Returns ``max|div(rho0 u)|`` over the
    coarse, fine, and fine-interface cells (recomputed independently -- ~machine precision).

    In mass-flux variables ``m = rho0 u`` the anelastic correction ``u = u* - grad(p)/rho0``
    is exactly ``m = m* - grad(p)``, ``div(m)=0``; the density weight cancels in the
    constraint.  The nest's lateral wall/relaxation BC is *replaced* by the true interface
    coupling here -- that is the point of the composite solve.
    """
    from .composite_poisson import composite_project_massflux_hz
    pg = parent.grid; xp = pg.xp; to = pg.backend.to_cpu
    ncx, ncy, nz, r = pg.nx, pg.ny, pg.nz, spec.refine
    ci0 = int(round(spec.x0 / pg.dx)); ci1 = int(round((spec.x0 + spec.Lx) / pg.dx))
    cj0 = int(round(spec.y0 / pg.dy)); cj1 = int(round((spec.y0 + spec.Ly) / pg.dy))
    if nest.grid.nz != nz or nest.grid.nx != r * (ci1 - ci0) or nest.grid.ny != r * (cj1 - cj0):
        raise ValueError("needs a cell-aligned, matched-z nest (build the spec with "
                         "NestSpec.aligned)")
    zc = np.asarray(to(pg.zc)); zf = np.asarray(to(pg.zf))
    dzc = zf[1:] - zf[:-1]; dzf = np.diff(zc)
    if parent.dynamics == "anelastic":
        rc = np.asarray(to(parent.rho0_c)); rw = np.asarray(to(parent.rho0_wface))
    else:                                              # Boussinesq: div(u)=0
        rc = np.ones(nz); rw = np.ones(nz + 1)
    rc3, rw3 = rc[None, None, :], rw[None, None, :]

    # extract host face velocities -> anelastic mass fluxes m = rho0 u
    mu_c = np.asarray(to(parent.state.u)) * rc3; mv_c = np.asarray(to(parent.state.v)) * rc3
    mw_c = np.asarray(to(parent.state.w)) * rw3
    mu_f = np.asarray(to(nest.state.u)) * rc3; mv_f = np.asarray(to(nest.state.v)) * rc3
    mw_f = np.asarray(to(nest.state.w)) * rw3
    per = bool(getattr(pg, "periodic", False)) if periodic_h is None else periodic_h

    res = composite_project_massflux_hz(mu_c, mv_c, mw_c, mu_f, mv_f, mw_f, ncx, nz, r,
                                        ci0, ci1, cj0, cj1, dzc, dzf, periodic_h=per,
                                        hx=pg.dx, hy=pg.dy, ncy=ncy)

    # recover u = m / rho0 and write back into each FlowState
    parent.state.u = xp.asarray(mu_c / rc3); parent.state.v = xp.asarray(mv_c / rc3)
    parent.state.w = xp.asarray(mw_c / rw3)
    nest.state.u = xp.asarray(mu_f / rc3); nest.state.v = xp.asarray(mv_f / rc3)
    nest.state.w = xp.asarray(mw_f / rw3)
    return {"div_coarse": res[0], "div_fine": res[1], "div_interface": res[2]}


def _shift_to_relative_frame(nest, cx: float, cy: float) -> None:
    """Galilean shift the nest into the storm-relative frame (subtract C).

    Vorticity and w are unchanged by a constant velocity offset, so the rotation
    diagnostics are identical; the storm (moving at ~C in the ground frame) becomes
    quasi-stationary in the nest, so a *fixed* nest keeps it centred."""
    nest.state.u = nest.state.u - cx
    nest.state.v = nest.state.v - cy
    nest._u0_face = nest._u0_face - cx
    nest._v0_face = nest._v0_face - cy
    for nm in ("u", "v"):
        if nm in nest._tgt:
            nest._tgt[nm] = nest._tgt[nm] - (cx if nm == "u" else cy)


def _snap_aligned_spec(parent, spec: NestSpec) -> NestSpec:
    """Snap an arbitrary nest footprint onto parent cell boundaries (and to the
    parent's vertical grid) so the nest is **cell-aligned, matched-z** -- the
    prerequisite of the composite projection (:meth:`NestSpec.aligned`)."""
    pg = parent.grid
    i0 = max(0, min(int(round(spec.x0 / pg.dx)), pg.nx))
    j0 = max(0, min(int(round(spec.y0 / pg.dy)), pg.ny))
    ncx = max(2, min(int(round(spec.Lx / pg.dx)), pg.nx - i0))
    ncy = max(2, int(round(spec.Ly / pg.dy)))
    ncy = max(2, min(ncy, pg.ny - j0))
    return _dc.replace(spec, x0=i0 * pg.dx, y0=j0 * pg.dy,
                       Lx=ncx * pg.dx, Ly=ncy * pg.dy,
                       nz=pg.nz, z_stretch=getattr(pg, "z_stretch", 1.0))


def run_concurrent_nest(parent, spec: NestSpec, window: float,
                        record_interval=None, capture_frames=False,
                        follow=False, storm_motion=None, two_way=False,
                        two_way_rate=0.5, progress=None,
                        les_boost=1.25, cfl=0.25,
                        composite_projection=False):
    """M3 phase 2: integrate the nest with **time-evolving** parent boundaries.

    The parent keeps stepping alongside the nest; each parent step the parent
    state is re-interpolated onto the nest and the nest sub-cycles (finer dt) with
    its border relaxed toward the parent target **interpolated linearly in time**.
    Fresh inflow keeps entering, so the storm is sustained beyond the
    frozen-boundary window -- still one-way and at fixed refinement.

    ``follow`` (M3 phase 2b, **storm-following nest**): run the nest in the
    storm-relative frame (Galilean shift by the storm motion ``C``) and slide the
    *sampled* parent region at ``C`` so the storm stays centred in a fixed nest --
    a feature that would otherwise advect out is now tracked.  ``storm_motion``
    ``(cx, cy)`` overrides the Bunkers estimate.

    ``composite_projection`` (**M4**): replace the two independent per-level
    anelastic projections with ONE composite solve over both levels' staggered
    mass fluxes (:func:`composite_project_two_level`), so ``div(rho0 u) = 0``
    holds **consistently across the coarse-fine interface** every nest sub-step.
    Requires a **cell-aligned, matched-z** nest (``NestSpec.aligned``, or an
    ``around`` footprint that is snapped here); with ``follow=True`` the nest is
    reconciled back to the ground frame for the joint solve and shifted back to
    the storm-relative frame afterwards, and the nest's velocity sponge is
    dropped (the interface coupling replaces it).  The solve is host-side
    (NumPy) -- GPU runs pay a host round-trip each sub-step.  Composable with
    ``two_way`` (the phase-3a injection then feeds the projected fine solution).

    Returns ``(nest, report)``; the report carries the worst-case interface
    divergence under ``rep["nest"]["composite_div_interface"]``.
    """
    from .core import StormSimulation as SS
    nest = NestedStormSimulation(parent, spec, les_boost=les_boost, cfl=cfl)
    if composite_projection:
        if nest.grid.nz != parent.grid.nz or abs(nest.grid.z_stretch - parent.grid.z_stretch) > 1e-12:
            raise ValueError(
                "composite_projection needs a MATCHED-Z nest (nz == parent nz, same "
                "z_stretch): rebuild the spec with NestSpec.aligned(..., nz=parent.grid.nz, "
                "z_stretch=parent.grid.z_stretch)")
        spec = _snap_aligned_spec(parent, spec)
        r = spec.refine
        if (nest.grid.nx, nest.grid.ny) != (r * round(spec.Lx / parent.grid.dx),
                                            r * round(spec.Ly / parent.grid.dy)):
            raise ValueError("nest grid %dx%d does not match the snapped footprint --"
                             " rebuild with NestSpec.aligned" % (nest.grid.nx, nest.grid.ny))
        nest._composite = True
    nest._capture_frames = bool(capture_frames); nest.frames = []
    g = nest.grid
    cx = cy = 0.0
    if follow:
        if storm_motion is None:
            from .soundings import bunkers_storm_motion
            cx, cy = bunkers_storm_motion(parent.base)
        else:
            cx, cy = storm_motion

    def _targets_at(pstate, t):
        # sample the parent at the (possibly moving) nest region, in the nest frame
        sp = _dc.replace(spec, x0=spec.x0 + cx * t, y0=spec.y0 + cy * t) if follow else spec
        tg = _interp_targets(pstate, parent.grid, g, sp)
        if follow:
            tg["u"] = tg["u"] - cx; tg["v"] = tg["v"] - cy
        return tg

    initial = _diag.initial_budgets(nest.state, nest.rho_ref)
    interval = record_interval or max(1, nest.cfg.output.interval_steps)
    tgt_prev = _targets_at(parent.state, 0.0)
    nest.set_target(tgt_prev)
    if follow:
        _shift_to_relative_frame(nest, cx, cy)   # nest + initial target -> relative frame
        nest.state.diagnose(nest.cfg)
    nest.tracker.update(nest.t, nest.state, g); SS._record(nest, initial)
    parent._Km = getattr(parent, "_Km", None)
    t = 0.0
    div_if_max = 0.0
    while t < window - 1e-9:
        dtp = parent._dt()
        if t + dtp > window:
            dtp = window - t
        if composite_projection:
            # parent predictor only; the composite solve below replaces its
            # per-level projection (it projects BOTH levels jointly)
            Km_p = parent._predictor(dtp)
            tgt_next = _targets_at(parent.state, t + dtp)
            t_sub = 0.0
            while t_sub < dtp - 1e-9:
                dtn = min(nest._dt(), dtp - t_sub)
                alpha = (t_sub + dtn) / dtp
                nest.set_target(_blend(tgt_prev, tgt_next, alpha))
                Km_n = nest._predictor(dtn)
                # frame reconciliation: the parent is in the GROUND frame, the
                # storm-following nest in the storm-relative frame -- reconcile
                # to the ground frame before forming the joint mass fluxes.
                if follow:
                    nest.state.u = nest.state.u + cx
                    nest.state.v = nest.state.v + cy
                div = composite_project_two_level(parent, nest, spec)
                div_if = float(div["div_interface"])
                div_if_max = max(div_if_max, div_if)
                if follow:
                    nest.state.u = nest.state.u - cx
                    nest.state.v = nest.state.v - cy
                nest._last_res = div_if  # interface residual is the honest solver metric
                nest._transport(dtn, Km_n)
                nest.step += 1; nest.t = float(nest.state.t); t_sub += dtn
                if nest.step % interval == 0:
                    nest.tracker.update(nest.t, nest.state, g); SS._record(nest, initial)
            parent._transport(dtp, Km_p)
        else:
            parent._step(dtp)                   # advance the parent
            tgt_next = _targets_at(parent.state, t + dtp)
            t_sub = 0.0
            while t_sub < dtp - 1e-9:
                dtn = min(nest._dt(), dtp - t_sub)
                alpha = (t_sub + dtn) / dtp
                nest.set_target(_blend(tgt_prev, tgt_next, alpha))
                nest._step(dtn)
                nest.step += 1; nest.t = float(nest.state.t); t_sub += dtn
                if nest.step % interval == 0:
                    nest.tracker.update(nest.t, nest.state, g); SS._record(nest, initial)
        # phase 3a: feed the nest's finer solution back onto the parent overlap
        if two_way:
            restrict_nest_to_parent(nest, parent, spec, cx, cy, t + dtp, rate=two_way_rate)
        tgt_prev = tgt_next
        t += dtp
        if progress:
            progress(t, window, nest.step)
    nest.tracker.update(nest.t, nest.state, g); SS._record(nest, initial)
    rep = SS._finalise(nest, initial)
    base_mode = ("concurrent + storm-following (phase 2b)" if follow
                 else "concurrent (phase 2: time-evolving parent boundary)")
    if two_way:
        base_mode += " + approximate two-way feedback (phase 3a)"
    if composite_projection:
        base_mode += " + composite two-level projection"
    rep["nest"] = {"dx_m": g.dx, "dz0_m": float(g.dz_c[0]),
                   "nx": g.nx, "ny": g.ny, "nz": g.nz, "refine": spec.refine,
                   "storm_motion": [cx, cy], "two_way": bool(two_way), "mode": base_mode}
    if composite_projection:
        rep["nest"]["composite_div_interface"] = div_if_max
    return nest, rep


__all__ = [
    "NestSpec", "build_nest_grid", "interpolate_state_to_nest",
    "interpolate_base_to_nest", "relaxation_weight", "NestedStormSimulation",
    "run_concurrent_nest", "restrict_nest_to_parent", "conservative_restrict",
    "interior_near_surface_zeta", "composite_project_two_level",
]
