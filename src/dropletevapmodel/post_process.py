import os
import yaml
import logging
import matplotlib.pyplot as plt
from datetime import datetime
from pathlib import Path
from dropletevapmodel._utils import dump_csv, plot_curves
from dropletevapmodel.model import DropletEvapModel

logger = logging.getLogger(__name__)

def post_processing(
        model: DropletEvapModel,
        data:dict,
        directory:str
        ) -> None:
    """
    Post processing

    :param model: (DropletEvapModel) model
    :param data: (dict) data
    :param directory: (str) directory name

    """

    plotting = model.simulation_config.output.plotting
    csv = model.simulation_config.output.csv
    timestep = model.simulation_config.modelling.timestep
    if plotting or csv:
        now = datetime.now()
        dt_string = now.strftime("%Y%m%d_%H%M%S")
        directory = 'res'+dt_string
        if not os.path.exists(os.path.join('output',directory)):
            os.makedirs(os.path.join('output',directory))
        with open(os.path.join('output',directory,'config.yml'), 'w') as outfile:
            yaml.dump(model.dict(), outfile, default_flow_style=False)
    if csv:
        logger.info('Dumping csv')
        dump_csv(res=data, directory=directory)
    if plotting:
        logger.info('Plotting curves')
        plot_curves(res=data, timestep=timestep, directory=directory)