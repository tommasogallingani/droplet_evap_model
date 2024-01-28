import tqdm
import numpy as np
from typing import Dict

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
        """
        Initialize the model

        :param model: Model
        :type model: DropletEvapModel

        """
        self.model = model

        self.sigma12, self.epsilon12, self.mass_g = self.eval_const()
        logger.info('Completed evaluation of constants')
        self.init_params = self.initialize_params()
        self.plotting = self.model.simulation_config.output.plotting
        self.csv = self.model.simulation_config.output.csv
        logger.info('Completed initialization of t_zero state')


    def eval_const(self):
        """
        Evaluate the constants of the model
        """
        sigma_d = self.model.fluid_properties.sigma_d
        sigma_g = self.model.gas_properties.sigma_g
        kb = self.model.simulation_config.environment.kb
        r = self.model.gas_properties.r
        eps_gas = self.model.gas_properties.eps_fact
        eps_liq = self.model.fluid_properties.eps_fact
        mm_g = self.model.gas_properties.mm_g
        pressure = self.model.simulation_config.environment.pressure
        # Zero conditions
        t_g_zero = self.model.simulation_config.t_zero.t_g_zero
        drop_conc = self.model.simulation_config.environment.drop_conc
        
        sigma12 = (sigma_d + sigma_g)/2
        epsilon_g = kb*eps_gas
        epsilon_d = kb*eps_liq
        epsilon12 = np.sqrt(epsilon_g*epsilon_d)
        mass_g = pressure*mm_g / \
            (drop_conc *
            r*t_g_zero)
        return sigma12, epsilon12, mass_g

    def initialize_params(self):
        """
        Initialize the parameters of the model for the zero state
        """
        params = {}
        params['d_d'] = self.model.simulation_config.t_zero.d_zero # droplet diameter m
        params['t_d'] = self.model.simulation_config.t_zero.t_d_zero  # droplet temperature K
        params['t_g'] = self.model.simulation_config.t_zero.t_g_zero  # gas temperature K
        # water content in gas phase ppm
        params['ppm'] = self.model.simulation_config.t_zero.water_content_g_zero
        params['rho_d'] = self.model.fluid_properties.rho_d  # liquid density g/m3
        params['t_layer'] = layer_temp(
            params['t_d'], params['t_g'])  # evaluate layer temperature
        # evaluate pressure saturation gas
        params['psat_g'] = sat_pressure_g(params['t_g'])
        params['pvap_g'] = vap_pressure_g(
            params['ppm'], self.model.simulation_config.environment.pressure, self.model.gas_properties.mm_g, self.model.fluid_properties.mm_d)
        params['rh'] = params['pvap_g']/params['psat_g']
        return params


    def evaluate_state(
            self,
            res_prec:dict,
            sigma12:float,
            epsilon12:float,
            mass_g:float
        ) -> Dict:
        """
        Evaluate the state of the model at the current timestamp

        :param res_prec: Previous state
        :param sigma12: Sigma12 constant
        :param epsilon12: Epsilon12 constant
        :param mass_g: Mass of the gas
        :return: Current state

        """
        kb = self.model.simulation_config.environment.kb
        r = self.model.gas_properties.r
        mm_g =  self.model.gas_properties.mm_g
        mm_d = self.model.fluid_properties.mm_d
        pressure = self.model.simulation_config.environment.pressure
        patm = self.model.simulation_config.environment.pressure_atm
        boiling_t = self.model.fluid_properties.boiling_t
        u_g = self.model.gas_properties.speed_g
        u_d = self.model.fluid_properties.speed_d
        timestep = self.model.simulation_config.modelling.timestep
        cp_d = self.model.fluid_properties.cp_l

        # Parameter for viscosity calculation
        viscosity_a_g = self.model.gas_properties.viscosity_a
        viscosity_b_g = self.model.gas_properties.viscosity_b
        viscosity_a_l = self.model.fluid_properties.viscosity_a
        viscosity_b_l = self.model.fluid_properties.viscosity_b
        # Parameter for heat capacity calculation gas
        cp_a_g = self.model.gas_properties.cp_a
        cp_b_g = self.model.gas_properties.cp_b
        cp_c_g = self.model.gas_properties.cp_c
        cp_d_g = self.model.gas_properties.cp_d
        cp_e_g = self.model.gas_properties.cp_e
        # Parameter for heat capacity calculation liquid
        cp_a_l = self.model.fluid_properties.cp_a
        cp_b_l = self.model.fluid_properties.cp_b
        cp_c_l = self.model.fluid_properties.cp_c
        cp_d_l = self.model.fluid_properties.cp_d
        cp_e_l = self.model.fluid_properties.cp_e
        # Parameter for thermal conductivity calculation gas
        k_a_g = self.model.gas_properties.k_a
        k_b_g = self.model.gas_properties.k_b
        k_c_g = self.model.gas_properties.k_c
        # Parameter for thermal conductivity calculation liquid
        k_a_l = self.model.fluid_properties.k_a
        k_b_l = self.model.fluid_properties.k_b
        k_c_l = self.model.fluid_properties.k_c


                                    
        res = res_prec.copy()
        # TODO Update liquid density as a function of temperature
        # eval layer temperature
        res['t_layer'] = layer_temp(res_prec['t_d'], res_prec['t_g'])
        res['kbtepsilon12'] = kb*res['t_layer'] / epsilon12
        res['omega'] = eval_omega(res['kbtepsilon12'])
        res['diff_vm'] = eval_diff_coeff(
            res['t_layer'], res['omega'], sigma12, pressure, mm_g, mm_d)
        res['vap_heat'] = eval_vap_heat(res_prec['t_d'])
        # surface vapour molar fraction at surface
        res['chi_vs'] = pressure/patm * \
            np.exp(res['vap_heat']*mm_d/r*(1/boiling_t-1/res_prec['t_d']))
        # surface vapour molar fraction in the gas
        res['chi_gs'] = 1 - res['chi_vs']
        res['chi_vg'] = res_prec['pvap_g']/pressure
        # vapour mass fraction in the gas
        res['y_vg'] = res['chi_vg']*mm_d / \
            (res['chi_vg']*mm_d+(1-res['chi_vg'])*mm_g)
        # vapour mass fraction at the surface
        res['y_vs'] = res['chi_vs']*mm_d/(res['chi_vs']*mm_d+res['chi_gs']*mm_g)
        res['Bm'] = (res['y_vs'] - res['y_vg']) / \
            (1-res['y_vs'])  # Spalding mass transfer number
        res['rho_g_g'] = pressure*mm_g/(r*res['t_layer'])
        res['rho_l_g'] = pressure*mm_d/(r*res['t_layer'])
        # average density
        res['rho_m_g'] = res['rho_g_g'] * \
            (1-res['chi_vs']) + res['chi_vs']*res['rho_l_g']

        # model convection
        res['mu_g_g'] = eval_viscosity(res['t_layer'], a=viscosity_a_g, b=viscosity_b_g)
        res['mu_l_g'] = eval_viscosity(res['t_layer'], a=viscosity_a_l, b=viscosity_b_l)
        # average viscosity
        res['mu_m_g'] = eval_m_viscosity(
            res['mu_g_g'], res['mu_l_g'], res['chi_vs'], mm_g, mm_d)
        # Reynold number
        res['Red'] = pressure/(r*res_prec['t_g'])*mm_g * \
            res_prec['d_d']*np.abs(u_g-u_d)/res['mu_m_g']*1.0e-3
        # Schmidt number at interface
        res['Scm'] = res['mu_m_g']/(res['rho_m_g']*res['diff_vm'])*1.0e7
        # Sherwood number
        res['Shm'] = 2+0.6*res['Red']**0.5*res['Scm']**(1/3)
        res['FM'] = (1+res['Bm'])**0.7*np.log(1+res['Bm'])/res['Bm']
        res['Shm_star'] = 2+(res['Shm']-2)/res['FM']
        # Mass losses
        res['md'] = np.pi*res_prec['d_d']*1.0e-4*res['diff_vm'] * \
            res['rho_m_g']*res['Shm_star']*np.log(1+res['Bm'])
        res['d_dt'] = -2*res['md']/(np.pi*res_prec['rho_d']*(res_prec['d_d'])**2)
        res['d_d'] = res_prec['d_d']+res['d_dt']*timestep
        # Evaluate humidity and droplet temperature
        res['cp_g_g'] = eval_cp(
            res['t_layer'],
            cp_a_g,
            cp_b_g,
            cp_c_g,
            cp_d_g,
            cp_e_g,
            mm_g
            )
        res['cp_l_g'] = eval_cp(
            res['t_layer'],
            cp_a_l,
            cp_b_l,
            cp_c_l,
            cp_d_l,
            cp_e_l,
            mm_d
            )
        # Jg/K Heat capacity mix ## TODO check y_vs
        res['cp_m_g'] = res['cp_l_g']*res['y_vs']+(1-res['y_vs'])*res['cp_g_g']
        # mW/mK thermal conductivity argon NIST with fitting excel
        res['k_g_g'] = eval_k(
            res['t_layer'],
            k_a_g,
            k_b_g,
            k_c_g
            )
        # mW/mK thermal conductivity vapour NIST with fitting excel
        res['k_l_g'] = eval_k(
            res['t_layer'],
            k_a_l,
            k_b_l,
            k_c_l
            )
        # argon to vapour wilke mixture calculation
        res['phi_g_g'] = eval_phi(res['mu_g_g'], res['mu_l_g'], mm_g, mm_d)
        # vapour to argon wilke mixture calculatio
        res['phi_l_g'] = eval_phi(res['mu_l_g'], res['mu_g_g'], mm_d, mm_g)
        res['k_m_g'] = res['chi_vs']*res['k_l_g']/((1-res['chi_vs'])*res['phi_l_g']+res['chi_vs'])+(
            1-res['chi_vs'])*res['k_g_g']/((res['chi_vs'])*res['phi_g_g']+1-res['chi_vs'])
        # Thermal balance
        res['beta'] = -res['md']*res['cp_m_g'] / \
            (2*np.pi*res['k_m_g']*1e-3*res['d_d'])
        # correction factor evaporation
        res['G'] = res['beta']/(np.exp(res['beta']-1))
        # Prandtl number
        res['Prm'] = res['mu_m_g']*res['cp_m_g']/(1.0e-3*res['k_m_g'])*1.0e3
        # Nusselt number
        res['Num'] = 2+0.6*res['Red']**0.5*res['Prm']**(1/3)
        # Abramzon-Sirignano model
        # Lewis number
        res['Lem'] = res['k_m_g']*1.0e-3 / \
            (res['cp_m_g']*res['diff_vm']*1.0e-4*res['rho_m_g'])
        res['PHI'] = (res['cp_l_g']/res['cp_g_g']*res['Shm']/res['Num'])/res['Lem']
        res['Bt'] = (1+res['Bm'])**res['PHI']-1
        res['FT'] = (1+res['Bt'])**0.7*(np.log(1+res['Bt']))/(res['Bt'])
        res['Num_star'] = 2+(res['Num']-2)/res['FT']
        res['G_star'] = np.log(1+res['Bt'])/res['Bt']
        # sensible heat transfer for droplet temperature change
        res['Qs'] = res['G_star']*np.pi*res['d_d']*res['Num_star']*1.0e-3 * \
            res['k_m_g']*(res['t_g']-res['t_d'])-res['md']*res['vap_heat']  # W
        # Temperature delta in droplet
        res['dT_dt'] = 6*res['Qs']/((res['d_d']**3)*np.pi*res['rho_d']*cp_d)  # K/s
        res['deltaT'] = res['dT_dt']*timestep  # K/step
        res['t_d'] = res_prec['t_d']+res['deltaT']

        # Taking into account the change of humidity due to evaporation
        res['ppm'] = (res_prec['ppm']*mass_g+res['md']*1.0e6*timestep)/mass_g
        res['pvap_g'] = vap_pressure_g(res['ppm'], pressure, mm_g, mm_d)
        res['psat_g'] = sat_pressure_g(res['t_g'])
        res['rh'] = res['pvap_g']/res['psat_g']
        return res
    

    def run(self):
        """
        Run the model
        """
        # State initialization
        res = {}
        res[0] = self.init_params
        # Parameters initialization
        max_it = int(self.model.simulation_config.modelling.max_iteration)
        timestep = self.model.simulation_config.modelling.timestep
        logger.info('Starting modelling')
        res[0]['time'] = 0
        for t in tqdm.tqdm(range(1, max_it)):
            # Calculate the state at the current timestamp
            res[t] = self.evaluate_state(res[t-1], self.sigma12, self.epsilon12, self.mass_g)
            # Add time info
            res[t]['time'] = t*timestep
            # Check if the result is complex
            for k, v in res[t].items():
                if isinstance(v, complex):
                    logger.warning(f"Found complex res for {k} variable at iteration {t} and timestamp {t*self.model['modelling']['timestep']}")
            # Break if the droplet diameter is negative or zero
            if res[t]['d_d'] <= 0 or res[t]['d_d'] == np.nan or res[t]['t_d'] <= 0:
                logger.info(f'Estimated evaporation time of {t*timestep} seconds')
                break
            # Max iteration reached
            if t==max_it-1:
                logger.warning('Maximum number of iteration reached')
        return res





