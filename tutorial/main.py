import pandas as pd
import numpy as np
import os
from datetime import datetime
from path import Path
import logging
import logging.config
import yaml
import argparse
import matplotlib.pyplot as plt
import tqdm
from dropletevapmodel import load_config, setup_logger


def eval_const(config):
    sigma12 = (config['liquid']['sigma_d']+config['gas_phase']['sigma_g'])/2
    epsilon_g = config['environment']['kb']*config['gas_phase']['eps_fact']
    epsilon_d = config['environment']['kb']*config['liquid']['eps_fact']
    epsilon12 = np.sqrt(epsilon_g*epsilon_d)
    mass_g = config['environment']['pressure']*config['gas_phase']['mm_g'] / \
        (config['environment']['drop_conc'] *
         config['gas_phase']['r']*config['t_zero']['t_g_zero'])
    return sigma12, epsilon12, mass_g


def layer_temp(t_drop, t_gas, alpha=1/3):
    return t_drop + alpha*(t_gas - t_drop)


def sat_pressure_g(t_gas, a=77.3450, b=0.0057, c=7235, d=8.2):
    return (np.exp(a+b*t_gas-c/t_gas))/(t_gas**d)


def vap_pressure_g(ppm, press, mm_gas, mm_drop):
    return (ppm*mm_gas*press) / (mm_drop*1.0e6 + mm_gas*ppm)


def eval_omega(kbtepsilon12, a=0.49924, b=574.81152, c=1, d=11267.71941, e=1, f=1.45781):
    return a+b/(c+kbtepsilon12*d)**(e/f)


def eval_diff_coeff(t_layer, omega, sigma12, p, mm_g, mm_d, a=1.86e-3, b=3/2, c=1/2, d=9.86923e-6, e=2):
    return a*(t_layer**(b))*(1/mm_d+1/mm_g)**(c)/(p*d*(sigma12**e)*omega)


def eval_vap_heat(t_d, a=-1.66282e-7, b=2.53533e-4, c=-0.14657, d=35.23021, e=-434.96806):
    return a*(t_d)**4+b*(t_d)**3+c*(t_d)**2+d*(t_d)+e


def eval_viscosity(t_layer, a, b):
    return t_layer*a+b


def eval_term(mu_g_g, mu_l_g, chi_vs, mm_g, mm_d):
    return mu_l_g/(1+(1-chi_vs)/chi_vs*((1+((mu_l_g/mu_g_g)**0.5)*((mm_d/mm_g)**0.25))**2)/(2.83*(1+mm_d/mm_g)**0.5)**0.5)


def eval_m_viscosity(mu_g_g, mu_l_g, chi_vs, mm_g, mm_d):
    term_a = eval_term(mu_g_g, mu_l_g, chi_vs, mm_g, mm_d)
    term_b = eval_term(mu_l_g, mu_g_g, 1-chi_vs, mm_d, mm_g)
    return term_a+term_b


def eval_cp(t_m, a, b, c, d, e, mm_x):
    return (a+b*(t_m)/1000+c*((t_m)/1000)**2+d*(t_m/1000)**3+e/((t_m/1000)**2))/mm_x


def eval_k(t_m, a, b, c):
    return a*t_m**2+b*t_m+c


def eval_phi(mu_g_g, mu_l_g, mm_g, mm_d):
    return (1+(mu_g_g/mu_l_g)**0.5*(mm_d/mm_g)**0.25)**2/(np.sqrt(8)*(1+(mm_g/mm_d))**0.5)


def initialize_params(sigma12, epsilon12, mass_g, config):
    arguments = {}
    arguments['d_d'] = config['t_zero']['d_zero']  # droplet diameter m
    arguments['t_d'] = config['t_zero']['t_d_zero']  # droplet temperature K
    arguments['t_g'] = config['t_zero']['t_g_zero']  # gas temperature K
    # water content in gas phase ppm
    arguments['ppm'] = config['t_zero']['water_content_g_zero']
    arguments['rho_d'] = config['liquid']['rho_d']  # liquid density g/m3
    arguments['t_layer'] = layer_temp(
        arguments['t_d'], arguments['t_g'])  # evaluate layer temperature
    # evaluate pressure saturation gas
    arguments['psat_g'] = sat_pressure_g(arguments['t_g'])
    arguments['pvap_g'] = vap_pressure_g(
        arguments['ppm'], config['environment']['pressure'], config['gas_phase']['mm_g'], config['liquid']['mm_d'])
    arguments['rh'] = arguments['pvap_g']/arguments['psat_g']
    return arguments


def evaluate_state(res_prec, sigma12, epsilon12, mass_g, config):
    kb = config['environment']['kb']
    r = config['gas_phase']['r']
    mm_g =  config['gas_phase']['mm_g']
    mm_d = config['liquid']['mm_d']
    pressure = config['environment']['pressure']
    patm = config['environment']['pressure_atm']
    boiling_t = config['liquid']['boiling_t']
    u_g = config['gas_phase']['speed_g']
    u_d = config['liquid']['speed_d']
    timestep = config['modelling']['timestep']
    cp_d = config['liquid']['cp_l']
                                
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
    res['mu_g_g'] = eval_viscosity(res['t_layer'], a=config['gas_phase']['viscosity_a'], b=config['gas_phase']['viscosity_b'])
    res['mu_l_g'] = eval_viscosity(res['t_layer'], a=config['liquid']['viscosity_a'], b=config['liquid']['viscosity_b'])
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
        res['t_layer'], config['gas_phase']['cp_a'], config['gas_phase']['cp_b'], config['gas_phase']['cp_c'], config['gas_phase']['cp_d'], config['gas_phase']['cp_e'], mm_g)
    res['cp_l_g'] = eval_cp(res['t_layer'], config['liquid']['cp_a'], config['liquid']['cp_b'], config['liquid']['cp_c'], config['liquid']['cp_d'], config['liquid']['cp_e'], mm_d)
    # Jg/K Heat capacity mix ## TODO check y_vs
    res['cp_m_g'] = res['cp_l_g']*res['y_vs']+(1-res['y_vs'])*res['cp_g_g']
    # mW/mK thermal conductivity argon NIST with fitting excel
    res['k_g_g'] = eval_k(res['t_layer'], config['gas_phase']['k_a'], config['gas_phase']['k_b'],  config['gas_phase']['k_c'])
    # mW/mK thermal conductivity vapour NIST with fitting excel
    res['k_l_g'] = eval_k(res['t_layer'], config['liquid']['k_a'], config['liquid']['k_b'],  config['liquid']['k_c'])
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


def model_evap(config, logger):
    sigma12, epsilon12, mass_g = eval_const(config)
    logger.info('Completed evaluation of constants')
    res = {}
    arguments = initialize_params(sigma12, epsilon12, mass_g, config)
    logger.info('Completed initialization of t_zero state')
    res[0] = arguments
    max_it = int(config['modelling']['max_iteration'])
    logger.info('Starting modelling')
    for t in tqdm.tqdm(range(1, max_it)):
        res[t] = evaluate_state(res[t-1], sigma12, epsilon12, mass_g, config)
        for k, v in res[t].items():
            if isinstance(v, complex):
                logger.warning(f"Found complex res for {k} variable at iteration {t} and timestamp {t*config['modelling']['timestep']}")
        if res[t]['d_d'] <= 0 or res[t]['d_d'] == np.nan or res[t]['t_d'] <= 0:
            logger.info(f'Estimated evaporation time of {t*config["modelling"]["timestep"]} seconds')
            break
        if t==max_it-1:
            logger.warning('Maximum number of iteration reached')
    if config['output']['plotting'] or config['output']['csv']:
        now = datetime.now()
        dt_string = now.strftime("%Y%m%d_%H%M%S")
        directory = 'res'+dt_string
        if not os.path.exists(os.path.join('output',directory)):
            os.makedirs(os.path.join('output',directory))
        with open(os.path.join('output',directory,'config.yml'), 'w') as outfile:
            yaml.dump(config, outfile, default_flow_style=False)
    if config['output']['csv']:
        logger.info(f'Saving results to csv')
        df = pd.DataFrame.from_dict(res, orient='index')
        df['time'] = df.index*config['modelling']['timestep']
        df.to_csv(os.path.join('output',directory,'res.csv'))
    if config['output']['plotting']:
        logger.info(f'Saving plots')
        time = []
        diameter = []
        temp = []
        ppm = []
        for key, val in res.items():
            time.append(config['modelling']['timestep']*key)
            diameter.append(val['d_d'])
            temp.append(val['t_d'])
            ppm.append(val['ppm'])
        fig, ax = plt.subplots(3, 1, figsize=(25,10))
        ax[0].scatter(x=time, y=diameter, label='Diameter')
        ax[0].set_ylim((0, max(diameter)))
        ax[0].set_ylabel('Diameter (m)')
        ax[0].set_xlabel('Time (s)')
        ax1 = ax[0].twinx()
        ax1.scatter(x=time, y=temp, c='k', label='Droplet temperature')
        ax1.set_ylabel('Temperature (K)')
        ev_c = -np.diff(np.power(np.array(diameter), 2))/config['modelling']['timestep']
        ax[1].scatter(x=time[1:], y=ev_c, c='r', label='Evaporation constant')
        ax[1].set_ylim((0, max(ev_c)))
        ax[1].set_ylabel('Evaporation constant (m^2/s)')
        ax[1].set_xlabel('Time (s)')
        ax2 = ax[1].twinx()
        ax2.scatter(x=time, y=ppm, c='g', label='Water content ppm')
        ax2.set_ylabel('Water content (ppm)')
        ax[2].scatter(x=time, y=np.power(np.array(diameter), 2), c='y', label='D-square')
        ax[2].set_ylim((0, max(np.power(np.array(diameter), 2))))
        ax[2].set_ylabel('Squared diameter (m^2)')
        ax[2].set_xlabel('Time (s)')
        fig.legend()
        fig.savefig(os.path.join('output',directory,'plots.png'), dpi=300)
        plt.draw()
        plt.show()
        



def main(args):

    logger = setup_logger(
        'modelling_droplet_evaporation',
        level=args.log_level,
        config_path=args.log_config
    )
    logger.info('Loading config data')
    config = load_config(args.config)
    if args.gas_type is not None:
        logger.info('Gas parameters will be taken from csv database file')
        config['gas_phase'] = pd.read_csv(os.path.join(args.data_path,'gas.csv'), index_col='property')[args.gas_type].to_dict()
    if args.liq_type is not None:
        logger.info('Liquid parameters will be taken from csv database file')
        config['liquid'] = pd.read_csv(os.path.join(args.data_path,'liq.csv'), index_col='property')[args.liq_type].to_dict()
    model_evap(config, logger)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Modelling'
    )
    parser.add_argument(
        '--config',
        help='config path',
        default='configs/config.yml'
    )
    parser.add_argument(
        '--gas_type',
        help='gas type',
        default='Ar'
    )
    parser.add_argument(
        '--liq_type',
        help='liquid type',
        default='H2O'
    )
    parser.add_argument(
        '--data_path',
        help='data path',
        default='configs'
    )
    parser.add_argument(
        '--log_level',
        help='logging level',
        default='INFO'
    )
    parser.add_argument(
        '--log_config',
        help='config log',
        default='configs/log.yml'
    )

    args = parser.parse_args()
    main(args)
