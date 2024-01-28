from pathlib import Path
import pandas as pd
import numpy as np
import os
import logging
import logging.config
import yaml
import matplotlib.pyplot as plt


def dump_csv(
        res:dict,
        directory:str
    ) -> None:
    """
    Dump results to csv

    :param res: (dict) results
    :param directory: (str) directory name

    """
    df = pd.DataFrame.from_dict(res).T
    df.to_csv(os.path.join('output',directory,'res.csv'), index=True)

def plot_curves(
        res:dict,
        timestep:float,
        directory:str
    ) -> None:
    """
    Plot curves

    :param res: (dict) results
    :param timestep: (float) timestep
    :param directory: (str) directory name

    """
    time = []
    diameter = []
    temp = []
    ppm = []
    for key, val in res.items():
        time.append(timestep*key)
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
    ev_c = -np.diff(np.power(np.array(diameter), 2))/timestep
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


def load_config(
        path: str
        ) -> dict:
    """
    Load YAML configuration file
    
    :param path: (str) path to YAML configuration file
    :return: (dict)
    """
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


def get_query_folder() -> Path:
    cf = Path(__file__).parent.absolute()
    qdir = cf / 'queries'
    return qdir


def get_query_from_sql_file(
        file_name: str,
        kwargs: dict = None
) -> str:
    with open(get_query_folder() / file_name, 'r') as f:
        query = f.read()
    # injecting kwargs if not None
    if kwargs is not None:
        query = query.format(**kwargs)
    return query
