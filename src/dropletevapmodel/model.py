import numpy as np
from typing import Dict
from scipy.integrate import solve_ivp

from ._utils import setup_logger
from .models.phy_model import DropletEvapModel
from .phy_utils import layer_temp, sat_pressure_g, vap_pressure_g, eval_omega, eval_diff_coeff, eval_vap_heat, eval_viscosity, eval_m_viscosity, eval_cp, eval_k, eval_phi

logger = setup_logger(
    'modelling_droplet_evaporation',
    level='INFO',
    config_path='configs/log.yml'
)


class EvapModel():
    def __init__(
            self,
            model:DropletEvapModel
        ) -> None:
        self.model = model
        self.sigma12, self.epsilon12, self.mass_g = self.eval_const()
        logger.info('Completed evaluation of constants')
        self.init_params = self.initialize_params()
        self.plotting = self.model.simulation_config.output.plotting
        self.csv = self.model.simulation_config.output.csv
        logger.info('Completed initialization of t_zero state')


    def eval_const(self):
        sigma_d = self.model.fluid_properties.sigma_d
        sigma_g = self.model.gas_properties.sigma_g
        kb = self.model.simulation_config.environment.kb
        r = self.model.gas_properties.r
        eps_gas = self.model.gas_properties.eps_fact
        eps_liq = self.model.fluid_properties.eps_fact
        mm_g = self.model.gas_properties.mm_g
        pressure = self.model.simulation_config.environment.pressure
        t_g_zero = self.model.simulation_config.t_zero.t_g_zero
        drop_conc = self.model.simulation_config.environment.drop_conc

        sigma12 = (sigma_d + sigma_g)/2
        epsilon_g = kb*eps_gas
        epsilon_d = kb*eps_liq
        epsilon12 = np.sqrt(epsilon_g*epsilon_d)
        mass_g = pressure*mm_g / (drop_conc * r*t_g_zero)
        return sigma12, epsilon12, mass_g

    def initialize_params(self):
        params = {}
        params['d_d'] = self.model.simulation_config.t_zero.d_zero
        params['t_d'] = self.model.simulation_config.t_zero.t_d_zero
        params['t_g'] = self.model.simulation_config.t_zero.t_g_zero
        params['ppm'] = self.model.simulation_config.t_zero.water_content_g_zero
        params['rho_d'] = self.model.fluid_properties.rho_d
        params['t_layer'] = layer_temp(params['t_d'], params['t_g'])
        params['psat_g'] = sat_pressure_g(params['t_g'])
        params['pvap_g'] = vap_pressure_g(
            params['ppm'], self.model.simulation_config.environment.pressure,
            self.model.gas_properties.mm_g, self.model.fluid_properties.mm_d)
        params['rh'] = params['pvap_g']/params['psat_g']
        return params


    def _compute_intermediates(self, d_d: float, t_d: float, t_g: float, ppm: float, rho_d: float) -> Dict:
        """Compute all intermediate quantities and ODE derivatives for the given state."""
        kb = self.model.simulation_config.environment.kb
        r = self.model.gas_properties.r
        mm_g = self.model.gas_properties.mm_g
        mm_d = self.model.fluid_properties.mm_d
        pressure = self.model.simulation_config.environment.pressure
        patm = self.model.simulation_config.environment.pressure_atm
        boiling_t = self.model.fluid_properties.boiling_t
        u_g = self.model.gas_properties.speed_g
        u_d = self.model.fluid_properties.speed_d
        cp_d = self.model.fluid_properties.cp_l
        viscosity_a_g = self.model.gas_properties.viscosity_a
        viscosity_b_g = self.model.gas_properties.viscosity_b
        viscosity_a_l = self.model.fluid_properties.viscosity_a
        viscosity_b_l = self.model.fluid_properties.viscosity_b
        cp_a_g = self.model.gas_properties.cp_a
        cp_b_g = self.model.gas_properties.cp_b
        cp_c_g = self.model.gas_properties.cp_c
        cp_d_g = self.model.gas_properties.cp_d
        cp_e_g = self.model.gas_properties.cp_e
        cp_a_l = self.model.fluid_properties.cp_a
        cp_b_l = self.model.fluid_properties.cp_b
        cp_c_l = self.model.fluid_properties.cp_c
        cp_d_l = self.model.fluid_properties.cp_d
        cp_e_l = self.model.fluid_properties.cp_e
        k_a_g = self.model.gas_properties.k_a
        k_b_g = self.model.gas_properties.k_b
        k_c_g = self.model.gas_properties.k_c
        k_a_l = self.model.fluid_properties.k_a
        k_b_l = self.model.fluid_properties.k_b
        k_c_l = self.model.fluid_properties.k_c

        s = {}
        s['t_layer'] = layer_temp(t_d, t_g)
        s['kbtepsilon12'] = kb * s['t_layer'] / self.epsilon12
        s['omega'] = eval_omega(s['kbtepsilon12'])
        s['diff_vm'] = eval_diff_coeff(s['t_layer'], s['omega'], self.sigma12, pressure, mm_g, mm_d)
        s['vap_heat'] = eval_vap_heat(t_d)
        s['chi_vs'] = pressure/patm * np.exp(s['vap_heat']*mm_d/r*(1/boiling_t - 1/t_d))
        s['chi_gs'] = 1 - s['chi_vs']
        pvap_g = vap_pressure_g(ppm, pressure, mm_g, mm_d)
        s['pvap_g'] = pvap_g
        s['chi_vg'] = pvap_g / pressure
        s['y_vg'] = s['chi_vg']*mm_d / (s['chi_vg']*mm_d + (1 - s['chi_vg'])*mm_g)
        s['y_vs'] = s['chi_vs']*mm_d / (s['chi_vs']*mm_d + s['chi_gs']*mm_g)
        s['Bm'] = (s['y_vs'] - s['y_vg']) / (1 - s['y_vs'])
        s['rho_g_g'] = pressure * mm_g / (r * s['t_layer'])
        s['rho_l_g'] = pressure * mm_d / (r * s['t_layer'])
        s['rho_m_g'] = s['rho_g_g']*(1 - s['chi_vs']) + s['chi_vs']*s['rho_l_g']
        s['mu_g_g'] = eval_viscosity(s['t_layer'], a=viscosity_a_g, b=viscosity_b_g)
        s['mu_l_g'] = eval_viscosity(s['t_layer'], a=viscosity_a_l, b=viscosity_b_l)
        s['mu_m_g'] = eval_m_viscosity(s['mu_g_g'], s['mu_l_g'], s['chi_vs'], mm_g, mm_d)
        s['Red'] = pressure/(r*t_g)*mm_g * d_d*np.abs(u_g - u_d)/s['mu_m_g']*1.0e-3
        s['Scm'] = s['mu_m_g'] / (s['rho_m_g']*s['diff_vm']) * 1.0e7
        s['Shm'] = 2 + 0.6*s['Red']**0.5 * s['Scm']**(1/3)
        s['FM'] = (1 + s['Bm'])**0.7 * np.log(1 + s['Bm']) / s['Bm']
        s['Shm_star'] = 2 + (s['Shm'] - 2) / s['FM']
        s['md'] = np.pi*d_d*1.0e-4 * s['diff_vm']*s['rho_m_g']*s['Shm_star']*np.log(1 + s['Bm'])
        s['d_dt'] = -2*s['md'] / (np.pi*rho_d*d_d**2)
        s['cp_g_g'] = eval_cp(s['t_layer'], cp_a_g, cp_b_g, cp_c_g, cp_d_g, cp_e_g, mm_g)
        s['cp_l_g'] = eval_cp(s['t_layer'], cp_a_l, cp_b_l, cp_c_l, cp_d_l, cp_e_l, mm_d)
        s['cp_m_g'] = s['cp_l_g']*s['y_vs'] + (1 - s['y_vs'])*s['cp_g_g']
        s['k_g_g'] = eval_k(s['t_layer'], k_a_g, k_b_g, k_c_g)
        s['k_l_g'] = eval_k(s['t_layer'], k_a_l, k_b_l, k_c_l)
        s['phi_g_g'] = eval_phi(s['mu_g_g'], s['mu_l_g'], mm_g, mm_d)
        s['phi_l_g'] = eval_phi(s['mu_l_g'], s['mu_g_g'], mm_d, mm_g)
        s['k_m_g'] = (
            s['chi_vs']*s['k_l_g'] / ((1 - s['chi_vs'])*s['phi_l_g'] + s['chi_vs']) +
            (1 - s['chi_vs'])*s['k_g_g'] / (s['chi_vs']*s['phi_g_g'] + 1 - s['chi_vs'])
        )
        s['beta'] = -s['md']*s['cp_m_g'] / (2*np.pi*s['k_m_g']*1e-3*d_d)
        s['G'] = s['beta'] / (np.exp(s['beta'] - 1))
        s['Prm'] = s['mu_m_g']*s['cp_m_g'] / (1.0e-3*s['k_m_g']) * 1.0e3
        s['Num'] = 2 + 0.6*s['Red']**0.5 * s['Prm']**(1/3)
        s['Lem'] = s['k_m_g']*1.0e-3 / (s['cp_m_g']*s['diff_vm']*1.0e-4*s['rho_m_g'])
        s['PHI'] = (s['cp_l_g']/s['cp_g_g'] * s['Shm']/s['Num']) / s['Lem']
        s['Bt'] = (1 + s['Bm'])**s['PHI'] - 1
        s['FT'] = (1 + s['Bt'])**0.7 * np.log(1 + s['Bt']) / s['Bt']
        s['Num_star'] = 2 + (s['Num'] - 2) / s['FT']
        s['G_star'] = np.log(1 + s['Bt']) / s['Bt']
        s['Qs'] = (
            s['G_star']*np.pi*d_d*s['Num_star']*1.0e-3*s['k_m_g']*(t_g - t_d)
            - s['md']*s['vap_heat']
        )
        s['dT_dt'] = 6*s['Qs'] / (d_d**3 * np.pi * rho_d * cp_d)
        s['dppm_dt'] = s['md'] * 1.0e6 / self.mass_g
        s['psat_g'] = sat_pressure_g(t_g)
        s['rh'] = pvap_g / s['psat_g']
        return s


    def _ode_rhs(self, t: float, y: np.ndarray) -> list:
        """ODE right-hand side for solve_ivp: returns [d(d_d)/dt, d(t_d)/dt, d(ppm)/dt]."""
        d_d, t_d, ppm = y
        if d_d <= 0 or t_d <= 0:
            return [0.0, 0.0, 0.0]
        t_g = self.init_params['t_g']
        rho_d = self.init_params['rho_d']
        s = self._compute_intermediates(d_d, t_d, t_g, ppm, rho_d)
        return [s['d_dt'], s['dT_dt'], s['dppm_dt']]


    def run(self):
        max_it = int(self.model.simulation_config.modelling.max_iteration)
        timestep = self.model.simulation_config.modelling.timestep
        t_span = (0.0, max_it * timestep)
        t_g = self.init_params['t_g']
        rho_d = self.init_params['rho_d']

        y0 = [self.init_params['d_d'], self.init_params['t_d'], self.init_params['ppm']]

        # Stop when d_d reaches this threshold; 1e-10 m is physically "fully evaporated"
        # and avoids the singularity in the derivatives as d_d → 0 (d_dt ∝ 1/d_d², dT_dt ∝ 1/d_d³).
        _D_THRESHOLD = 1e-10

        def event_stop(t, y):
            return y[0] - _D_THRESHOLD
        event_stop.terminal = True
        event_stop.direction = -1

        logger.info('Starting modelling with solve_ivp (Radau)')
        sol = solve_ivp(
            self._ode_rhs,
            t_span,
            y0,
            method='Radau',
            events=event_stop,
            rtol=1e-6,
            atol=1e-10,
        )

        if not sol.success:
            logger.warning(f'solve_ivp did not converge: {sol.message}')

        if sol.t_events[0].size > 0:
            logger.info(f'Estimated evaporation time: {sol.t_events[0][0]:.6f} s')
        else:
            logger.warning('Maximum simulation time reached without complete evaporation')

        # Reconstruct full output dict with all intermediate quantities
        res = {}
        for i, t in enumerate(sol.t):
            d_d, t_d, ppm = sol.y[:, i]
            state = self._compute_intermediates(d_d, t_d, t_g, ppm, rho_d)
            state['d_d'] = d_d
            state['t_d'] = t_d
            state['ppm'] = ppm
            state['t_g'] = t_g
            state['rho_d'] = rho_d
            state['time'] = t
            for k, v in state.items():
                if isinstance(v, complex):
                    logger.warning(f'Found complex result for {k} at t={t:.6f}')
            res[i] = state

        return res
