import pandas as pd
import numpy as np
import os
from datetime import datetime
from pathlib import Path
import logging
import logging.config
import yaml
import argparse
import matplotlib.pyplot as plt
import tqdm
from dropletevapmodel import setup_logger, load_config, post_processing
from dropletevapmodel import DropletEvapModel,EvapModel

logger = setup_logger(
    'Main',
    level='INFO',
    config_path='configs/log.yml'
)

def main(args):
    logger.info('Loading config data')
    config = load_config(args.config)
    if args.gas_type is not None:
        logger.info('Gas parameters will be taken from csv database file')
        config['gas_properties'] = pd.read_csv(os.path.join(args.data_path,'gas.csv'), index_col='property')[args.gas_type].to_dict()
    if args.liq_type is not None:
        logger.info('Liquid parameters will be taken from csv database file')
        config['fluid_properties'] = pd.read_csv(os.path.join(args.data_path,'liq.csv'), index_col='property')[args.liq_type].to_dict()
    
    logger.info('Running model')
    model = DropletEvapModel(**config)
    res = EvapModel(model=model).run()
    logger.info('Saving results')
    # Post processing
    post_processing(model=model, data=res, directory='res')


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
