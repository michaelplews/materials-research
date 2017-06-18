"""Unit tests for XAS.py"""

import unittest, sys, os, numpy as np, pandas as pd

wdir = os.path.dirname(__file__) # Find the current working directory
sys.path.append("..")
import XAS as xas

class XAS_init_tests(unittest.TestCase):

    def test__MDAdatafile_init(self):
        sample_data = [
            (wdir + "/test_data/htxs_scan.881", "SSRL", # Sample Data for SSRL
                ['Rounded Energy / eV',                 # index
                 'mono','Seconds',
                 'i0','tey','cmact','spare',                 'fy','ref','i1','SRS_i0','SRS_tey',
                 'Encoder','Grating_Encoder','Energy_encoder']
             )
            ]
        for filename, flavour, columns in sample_data:
            test = xas._MDAdatafile(filename, flavour=flavour)
            # Assert the output is of pandas dataframe type
            self.assertEqual(type(test.dataframe), pd.core.frame.DataFrame)
            # Assert index is as expected
            self.assertEqual(test.dataframe.index.name, columns[0])
            # Assert columns are as expected
            for i in range(1, len(columns)):
                self.assertEqual(test.dataframe.columns.values.tolist()[i-1], columns[i]) 
            # Assert number of rows is as expected
            self.assertEqual(len(test.dataframe), 288)

    def test_SSRL82_init(self):
        sample_data = [
            (wdir + "/test_data/", "htxs_scan", 881, 881) # Sample Data for SSRL82
            ]
        expected_columns = [
            'Rounded Energy / eV', # index
            'Energy / eV',
            'REF',
            'TFY',
            'sTFY',
            'TEY',
            'AEY'
            ]
        for d, b, s, e in sample_data:
            test = xas.SSRL82(directory=d, basename=b, start=s, end=e)
            # Assert the output is of pandas dataframe type
            self.assertEqual(type(test.processed_dataframe), pd.core.frame.DataFrame)
            # Assert index is as expected
            self.assertEqual(test.processed_dataframe.index.name, expected_columns[0])
            # Assert columns are as expected
            for i in range(1, len(expected_columns)):
                self.assertEqual(test.processed_dataframe.columns.values.tolist()[i-1], expected_columns[i]) 
            # Assert number of rows is as expected
            self.assertEqual(len(test.processed_dataframe), 268)
