import unittest
import numpy as np
import os
import pandas as pd
from dropletevapmodel import EvapModel, DropletEvapModel
from dropletevapmodel._utils import load_config

class TestEvapModel(unittest.TestCase):

    def setUp(self):
        # Set up the test environment, create a DropletEvapModel instance for testing
        # Ensure that you pass valid arguments to initialize your DropletEvapModel instance
        config = load_config('tests/assets/config.yml')
        config['gas_properties'] = pd.read_csv(os.path.join('tests/assets/gas.csv'), index_col='property')['Ar'].to_dict()
        config['fluid_properties'] = pd.read_csv(os.path.join('tests/assets/liq.csv'), index_col='property')['H2O'].to_dict()
    
        self.droplet_model = DropletEvapModel(**config)
        self.evap_model = EvapModel(self.droplet_model)

    def test_eval_const(self):
        # Test the eval_const method
        sigma12, epsilon12, mass_g = self.evap_model.eval_const()
        self.assertIsInstance(sigma12, (int, float))
        self.assertIsInstance(epsilon12, (int, float))
        self.assertIsInstance(mass_g, (int, float))

    def test_initialize_params(self):
        # Test the initialize_params method
        init_params = self.evap_model.initialize_params()
        self.assertIsInstance(init_params, dict)
        self.assertEqual(len(init_params), 9)  # Update this based on the actual number of expected parameters

    def test_run(self):
        # Test the run method
        result = self.evap_model.run()
        self.assertIsInstance(result, dict)
        self.assertTrue('time' in result[0])  # Check if 'time' key is present in the result

        self.assertEqual(np.round(result[max(result.keys())]['time'], 9), np.round(0.9992800000000001,9))  # Check if the last time step is equal to expected value

