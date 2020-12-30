import pandas as pd
import numpy as np
from path import Path
import logging, logging.config
import yaml
import argparse


def load_config(path: str) -> dict:
    path_config = Path(path)
    with open(path_config, 'r') as fp:
        config = yaml.load(fp, Loader=yaml.FullLoader)
    return config


def setup_logger(name: str, level="DEBUG", config_path=None) -> logging.Logger:
    """
    Create a logger with a given name
    Args:
        name: (str) logger name
        level: (str) logger level (default DEBUG)
        config_path: (str) Default None. YAML logger configuration file. If no configuration is provided,
                    a basic default is used.
    Returns:
        logger
    """

    logger = logging.getLogger(name)

    if config_path is not None:
        logging_config = load_config(config_path)
        logging.config.dictConfig(logging_config)
    else:
        kwargs = {
            "format": "%(asctime)s - %(name)-7s - %(levelname)s - %(message)s",
            "datefmt": "%Y-%m-%d %H:%M:%S"
        }

        logging.basicConfig(**kwargs)

    logger.setLevel(level)
    return logger
def eval_const(config):
    sigma12 = (config['liquid']['sigma_d']+config['gas_phase']['sigma_g'])/2
    epsilon_g = config['environment']['kb']*config['gas_phase']['eps_fact']
    epsilon_d = config['environment']['kb']*config['liquid']['eps_fact']
    epsilon12 = np.sqrt(epsilon_g*epsilon_d)
    mass_g = config['environment']['pressure']*config['gas_phase']['mm_g']/(config['environment']['drop_conc']*config['gas_phase']['r']*config['t_zero']['t_g_zero'])
    return sigma12, epsilon12, mass_g

def layer_temp(t_drop, t_gas, alpha=1/3):
    return t_drop + alpha*(t_gas - t_drop)

def sat_pressure_g(t_gas, a=77.3450, b=0.0057, c = 7235, d = 8.2 ):
    return (np.exp(a+b*t_gas-c/t_gas))/(t_gas**d)

def vap_pressure_g(ppm, press, mm_gas, mm_drop):
    return (ppm*mm_gas*press) / (mm_drop*1.0e6 + mm_gas*ppm)

def eval_omega(kbtepsilon12, a =0.49924, b=574.81152,c=1,d=11267.71941, e=1 ,f=1.45781):
    return a+b/(c+kbtepsilon12*d)**(e/f)

def eval_diff_coeff(t_layer, omega, sigma12, p, mm_g, mm_d, a=1.86e-3,b=3/2,c=1/2,d=9.86923e-6, e=2):
    return a*(t_layer**(b))*(1/mm_d+1/mm_g)**(c)/(p*d*(sigma12**e)*omega)

def eval_vap_heat(t_d, a=-1.66282e-7, b= 2.53533e-4, c=-0.14657, d=35.23021, e=-434.96806):
    return a*(t_d)**4+b*(t_d)**3+c*(t_d)**2+d*(t_d)+e



def initialize_params(sigma12, epsilon12, mass_g, config):
    arguments={}
    arguments['d_d'] = config['t_zero']['d_zero'] # droplet diameter m
    arguments['t_d'] = config['t_zero']['t_d_zero'] # droplet temperature K
    arguments['t_g'] = config['t_zero']['t_g_zero'] # gas temperature K
    arguments['ppm'] = config['t_zero']['water_content_g_zero'] # water content in gas phase ppm
    arguments['rho_d'] = config['liquid']['rho_d'] # liquid density g/m3
    arguments['t_layer'] = layer_temp(arguments['t_d'], arguments['t_g']) # evaluate layer temperature
    arguments['psat_g'] = sat_pressure_g(arguments['t_g']) # evaluate pressure saturation gas
    arguments['pvap_g'] = vap_pressure_g(arguments['ppm'], config['environment']['pressure'], config['gas_phase']['mm_g'], config['liquid']['mm_d'])
    arguments['rh'] = arguments['pvap_g']/arguments['psat_g']
    return arguments

def evaluate_state(res_prec, sigma12, epsilon12, mass_g, kb, mm_g, mm_d, pressure):
    res=dict.fromkeys(res_prec.keys())
    # eval layer temperature
    res['t_layer'] = layer_temp(res_prec['t_d'], res_prec['t_g'])
    res['kbtepsilon12'] = kb*res['t_layer'] /epsilon12
    res['omega'] = eval_omega(res['kbtepsilon12'])
    res['diff_vm'] = eval_diff_coeff(res['t_layer'], res['omega'], sigma12, pressure, mm_g, mm_d)
    res['vap_heat'] = eval_vap_heat(res_prec['t_d'])

    return res

def model_evap(config):
    sigma12, epsilon12, mass_g = eval_const(config)
    res = {}
    arguments = initialize_params(sigma12, epsilon12, mass_g, config)
    time=int(np.floor(config['modelling']['duration']/config['modelling']['timestep']))
    res[0] = arguments
    for t in range(1, time):
        res[t] = evaluate_state(res[t-1], sigma12, epsilon12, mass_g, config['environment']['kb'], config['gas_phase']['mm_g'], config['liquid']['mm_d'], config['environment']['pressure'])

def main(args):

    logger = setup_logger(
        'modelling_droplet_evaporation',
        level=args.log_level,
        config_path=args.log_config
    )
    logger.info('Loading config data')
    config = load_config(args.config)
    model_evap(config)



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