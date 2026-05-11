from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .mesh import AxisymmetricMesh, create_axisymmetric_mesh
from .transport import placeholder_diffusivity_m2_s

Array = np.ndarray


@dataclass
class LowMachState:
    """Cell-centered primitive/conserved state on the axisymmetric r-z mesh.

    Required array shape is ``(nr, nz)`` for scalar and velocity fields.
    Species mass fractions are stored as ``Y[species_name]``.
    """

    rho_kg_m3: Array
    u_r_m_s: Array
    u_z_m_s: Array
    T_K: Array
    Y: dict[str, Array]
    pressure_perturbation_Pa: Array | None = None
    time_s: float = 0.0


@dataclass
class TransportProperties:
    """Mixture transport fields used by the discrete operators."""

    mu_kg_m_s: Array
    lambda_W_m_K: Array
    cp_J_kg_K: Array
    species_diffusivity_m2_s: dict[str, Array]


@dataclass
class ChemistrySourceTerms:
    """Chemical source terms in conservative units.

    ``omega_dot_kg_m3_s[k]`` is the mass production rate of species k.
    ``heat_release_W_m3`` is positive for exothermic heat release.
    """

    omega_dot_kg_m3_s: dict[str, Array] = field(default_factory=dict)
    heat_release_W_m3: Array | None = None


@dataclass
class LowMachRHS:
    """Right-hand side returned by ``AxisymmetricLowMachSolver.compute_rhs``."""

    drho_dt: Array
    du_r_dt: Array
    du_z_dt: Array
    dT_dt: Array
    dY_dt: dict[str, Array]
    divergence_1_s: Array


ChemistrySourceCallback = Callable[[LowMachState, AxisymmetricMesh, dict], ChemistrySourceTerms]


def gravity_source_terms_enabled(config: dict) -> bool:
    return bool(config.get("gravity", {}).get("enabled", False))


class AxisymmetricLowMachSolver:
    r"""Finite-volume-like low-Mach axisymmetric reacting-flow framework.

    Continuous equations represented by this discretization, for no swirl:

    Continuity:
        ``d rho/dt + div_axisym(rho u) = 0``

    Momentum, with thermodynamic pressure removed from acoustics and optional
    perturbation pressure ``pi``:
        ``d(rho u_r)/dt + div_axisym(rho u_r u) = -d pi/dr
           + div(tau)_r + S_r``
        ``d(rho u_z)/dt + div_axisym(rho u_z u) = -d pi/dz
           + div(tau)_z + S_z``

    Temperature form of sensible enthalpy equation:
        ``rho cp (dT/dt + u.grad T) = div(lambda grad T)
           + qdot + S_T``

    Species:
        ``d(rho Y_k)/dt + div_axisym(rho Y_k u)
           = div_axisym(rho D_k grad Y_k) + omega_dot_k``

    The class currently performs an explicit predictor update. It does not yet
    solve the low-Mach pressure projection or enforce the thermodynamic
    divergence constraint; those extension points are isolated in
    ``pressure_gradient`` and ``low_mach_divergence_constraint``.
    """

    def __init__(
        self,
        config: dict,
        mesh: AxisymmetricMesh | None = None,
        chemistry_source_callback: ChemistrySourceCallback | None = None,
    ):
        self.config = config
        self.mesh = mesh if mesh is not None else create_axisymmetric_mesh(config["mesh"])
        self.chemistry_source_callback = chemistry_source_callback

    def step(self, state: LowMachState | dict, dt_s: float) -> LowMachState:
        """Advance one explicit predictor step.

        This is intentionally conservative and minimal. A production low-Mach
        DNS solver should replace this method with Strang chemistry splitting,
        pressure projection, boundary-condition handling, and implicit or
        semi-implicit diffusion as needed.
        """
        state_obj = self.as_state(state)
        rhs = self.compute_rhs(state_obj)
        rho = np.maximum(state_obj.rho_kg_m3 + dt_s * rhs.drho_dt, 1.0e-12)
        Y = {name: values + dt_s * rhs.dY_dt[name] for name, values in state_obj.Y.items()}
        Y = self.normalize_species(Y)
        return LowMachState(
            rho_kg_m3=rho,
            u_r_m_s=state_obj.u_r_m_s + dt_s * rhs.du_r_dt,
            u_z_m_s=state_obj.u_z_m_s + dt_s * rhs.du_z_dt,
            T_K=np.maximum(state_obj.T_K + dt_s * rhs.dT_dt, 1.0),
            Y=Y,
            pressure_perturbation_Pa=state_obj.pressure_perturbation_Pa,
            time_s=state_obj.time_s + dt_s,
        )

    def compute_rhs(self, state: LowMachState | dict) -> LowMachRHS:
        """Assemble continuity, momentum, energy, and species RHS terms."""
        state_obj = self.as_state(state)
        props = self.transport_properties(state_obj)
        chemistry = self.chemistry_source_terms(state_obj)

        rho = state_obj.rho_kg_m3
        ur = state_obj.u_r_m_s
        uz = state_obj.u_z_m_s
        continuity_flux_r = rho * ur
        continuity_flux_z = rho * uz
        drho_dt = -self.axisymmetric_divergence(continuity_flux_r, continuity_flux_z)

        dpi_dr, dpi_dz = self.pressure_gradient(state_obj)
        visc_r, visc_z = self.viscous_divergence(state_obj, props)
        body_r, body_z = self.body_force_acceleration(state_obj)
        adv_mom_r = self.axisymmetric_divergence(rho * ur * ur, rho * ur * uz)
        adv_mom_z = self.axisymmetric_divergence(rho * uz * ur, rho * uz * uz)
        mom_rhs_r = -adv_mom_r - dpi_dr + visc_r
        mom_rhs_z = -adv_mom_z - dpi_dz + visc_z
        du_r_dt = (mom_rhs_r - ur * drho_dt) / np.maximum(rho, 1.0e-12) + body_r
        du_z_dt = (mom_rhs_z - uz * drho_dt) / np.maximum(rho, 1.0e-12) + body_z

        dT_dr, dT_dz = self.gradient(state_obj.T_K)
        conductive = self.axisymmetric_divergence(props.lambda_W_m_K * dT_dr, props.lambda_W_m_K * dT_dz)
        heat_release = (
            chemistry.heat_release_W_m3
            if chemistry.heat_release_W_m3 is not None
            else np.zeros_like(state_obj.T_K)
        )
        dT_dt = (
            -ur * dT_dr
            - uz * dT_dz
            + (conductive + heat_release) / np.maximum(rho * props.cp_J_kg_K, 1.0e-12)
        )

        dY_dt: dict[str, Array] = {}
        for name, Yk in state_obj.Y.items():
            dY_dr, dY_dz = self.gradient(Yk)
            Dk = props.species_diffusivity_m2_s[name]
            diffusive = self.axisymmetric_divergence(rho * Dk * dY_dr, rho * Dk * dY_dz)
            omega = chemistry.omega_dot_kg_m3_s.get(name, np.zeros_like(Yk))
            conservative_rhs = -self.axisymmetric_divergence(rho * Yk * ur, rho * Yk * uz) + diffusive + omega
            dY_dt[name] = (conservative_rhs - Yk * drho_dt) / np.maximum(rho, 1.0e-12)

        divergence = self.axisymmetric_divergence(ur, uz)
        return LowMachRHS(
            drho_dt=drho_dt,
            du_r_dt=du_r_dt,
            du_z_dt=du_z_dt,
            dT_dt=dT_dt,
            dY_dt=dY_dt,
            divergence_1_s=divergence,
        )

    def chemistry_source_terms(self, state: LowMachState) -> ChemistrySourceTerms:
        """Return chemical source terms; Cantera coupling is inserted here.

        Expected future Cantera workflow:
        1. Build or reuse one ``ct.Solution`` per worker/thread.
        2. For each cell, set ``gas.TPY = T, p0, Y``.
        3. Read ``gas.net_production_rates`` and molecular weights.
        4. Convert to kg/m3/s and compute heat release from partial molar
           enthalpies.

        The current default is non-reacting, so Stage 1 and tests remain
        runnable without Cantera.
        """
        if self.chemistry_source_callback is not None:
            return self.chemistry_source_callback(state, self.mesh, self.config)
        zeros = np.zeros_like(state.T_K)
        return ChemistrySourceTerms(
            omega_dot_kg_m3_s={name: zeros.copy() for name in state.Y},
            heat_release_W_m3=zeros,
        )

    def transport_properties(self, state: LowMachState) -> TransportProperties:
        """Build transport fields for unity or non-unity Lewis-number modes."""
        T = state.T_K
        mu_ref = float(self.config.get("transport", {}).get("mu_ref_kg_m_s", 2.0e-5))
        cp_ref = float(self.config.get("transport", {}).get("cp_ref_J_kg_K", 14000.0))
        pr = float(self.config.get("transport", {}).get("prandtl_number", 0.7))
        mu = mu_ref * np.maximum(T / 298.0, 0.1) ** 0.7
        cp = np.full_like(T, cp_ref)
        conductivity = mu * cp / pr
        alpha = conductivity / np.maximum(state.rho_kg_m3 * cp, 1.0e-12)

        unity_lewis = bool(self.config.get("transport", {}).get("unity_lewis", True))
        diffusivity: dict[str, Array] = {}
        for name in state.Y:
            if unity_lewis:
                diffusivity[name] = alpha.copy()
            else:
                diffusivity[name] = np.full_like(T, placeholder_diffusivity_m2_s(float(np.mean(T))))
        return TransportProperties(
            mu_kg_m_s=mu,
            lambda_W_m_K=conductivity,
            cp_J_kg_K=cp,
            species_diffusivity_m2_s=diffusivity,
        )

    def pressure_gradient(self, state: LowMachState) -> tuple[Array, Array]:
        """Gradient of perturbation pressure ``pi``.

        The thermodynamic pressure is spatially uniform in the low-Mach model.
        If no perturbation pressure is supplied, the explicit predictor omits
        pressure forces. A projection method can compute and store ``pi`` here.
        """
        if state.pressure_perturbation_Pa is None:
            zeros = np.zeros_like(state.T_K)
            return zeros, zeros
        return self.gradient(state.pressure_perturbation_Pa)

    def low_mach_divergence_constraint(self, state: LowMachState) -> Array:
        """Reserved thermodynamic divergence constraint for projection.

        For reacting low-Mach flow, ``div(u)`` is set by heat release, species
        diffusion, and equation-of-state consistency rather than by acoustics.
        The current placeholder returns zero, i.e. incompressible projection.
        """
        return np.zeros_like(state.T_K)

    def body_force_acceleration(self, state: LowMachState) -> tuple[Array, Array]:
        """Return body-force acceleration components in m/s2."""
        zeros = np.zeros_like(state.T_K)
        if not gravity_source_terms_enabled(self.config):
            return zeros, zeros
        g = float(self.config.get("gravity", {}).get("acceleration_m_s2", 9.80665))
        return zeros, -g * np.ones_like(state.T_K)

    def viscous_divergence(self, state: LowMachState, props: TransportProperties) -> tuple[Array, Array]:
        """Divergence of Newtonian stress tensor for axisymmetric no-swirl flow."""
        ur = state.u_r_m_s
        uz = state.u_z_m_s
        mu = props.mu_kg_m_s
        dur_dr, dur_dz = self.gradient(ur)
        duz_dr, duz_dz = self.gradient(uz)
        div_u = self.axisymmetric_divergence(ur, uz)
        r = np.maximum(self.mesh.R, max(self.mesh.dr, 1.0e-12))
        ur_over_r = np.where(self.mesh.R > 0.0, ur / r, dur_dr)

        tau_rr = 2.0 * mu * dur_dr - (2.0 / 3.0) * mu * div_u
        tau_zz = 2.0 * mu * duz_dz - (2.0 / 3.0) * mu * div_u
        tau_tt = 2.0 * mu * ur_over_r - (2.0 / 3.0) * mu * div_u
        tau_rz = mu * (dur_dz + duz_dr)

        div_tau_r = self.axisymmetric_divergence(tau_rr, tau_rz) - np.where(
            self.mesh.R > 0.0,
            tau_tt / r,
            0.0,
        )
        div_tau_z = self.axisymmetric_divergence(tau_rz, tau_zz)
        div_tau_r[0, :] = 0.0
        return div_tau_r, div_tau_z

    def axisymmetric_divergence(self, flux_r: Array, flux_z: Array) -> Array:
        """Compute ``1/r d(r F_r)/dr + d(F_z)/dz`` on cell centers."""
        r = self.mesh.R
        r_flux = r * flux_r
        d_r_flux_dr = self.derivative_r(r_flux)
        d_flux_z_dz = self.derivative_z(flux_z)
        div = np.empty_like(flux_r)
        positive_r = r > 0.0
        div[positive_r] = d_r_flux_dr[positive_r] / r[positive_r] + d_flux_z_dz[positive_r]
        div[~positive_r] = 2.0 * self.derivative_r(flux_r)[~positive_r] + d_flux_z_dz[~positive_r]
        return div

    def gradient(self, scalar: Array) -> tuple[Array, Array]:
        return self.derivative_r(scalar), self.derivative_z(scalar)

    def derivative_r(self, scalar: Array) -> Array:
        if scalar.shape[0] < 2:
            return np.zeros_like(scalar)
        return np.gradient(scalar, self.mesh.dr, axis=0, edge_order=1)

    def derivative_z(self, scalar: Array) -> Array:
        if scalar.shape[1] < 2:
            return np.zeros_like(scalar)
        return np.gradient(scalar, self.mesh.dz, axis=1, edge_order=1)

    @staticmethod
    def normalize_species(Y: dict[str, Array]) -> dict[str, Array]:
        clipped = {name: np.clip(values, 0.0, 1.0) for name, values in Y.items()}
        total = np.zeros_like(next(iter(clipped.values()))) if clipped else None
        if total is None:
            return clipped
        for values in clipped.values():
            total += values
        total = np.maximum(total, 1.0e-30)
        return {name: values / total for name, values in clipped.items()}

    @staticmethod
    def as_state(state: LowMachState | dict) -> LowMachState:
        if isinstance(state, LowMachState):
            return state
        return LowMachState(
            rho_kg_m3=state["rho_kg_m3"],
            u_r_m_s=state["u_r_m_s"],
            u_z_m_s=state["u_z_m_s"],
            T_K=state["T_K"],
            Y=state["Y"],
            pressure_perturbation_Pa=state.get("pressure_perturbation_Pa"),
            time_s=float(state.get("time_s", 0.0)),
        )
