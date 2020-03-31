"""Unit tests for XAS.py"""

import unittest, sys, os, numpy as np, pandas as pd

wdir = os.path.dirname(__file__) # Find the current working directory
sys.path.append("..")
sys.path.append(".")
import cabanapy.XAS as xas

class XAS_init_tests(unittest.TestCase):

    def test__MDAdatafile_init(self):
        sample_data = [
            (wdir + "/test_data/Blank_C_tape.886", "SSRL", # Sample Data for SSRL
            [
                'Rounded Energy / eV',  # index
                'mono','tey','aey', 'pey', 'tfy','refy',
            ],
            288),
            (wdir + "/test_data/JLApr16.0001", "IDC4", # Sample Data or 4-ID-C
            [
                'Rounded Energy / eV', # index
                '[1-D Positioner 1]  4idc1:m13.VAL\t 7T Sample Z\t LINEAR\t mm\t 4idc1:m13.RBV\t 7T Sample Z\t mm',
                '[1-D Detector   9]  4idc1:scaler1_calc4.VAL\t \t ',
                '[1-D Detector  10]  4idc1:scaler1_calc5.VAL\t \t ',
                '[1-D Detector  11]  4idc1:scaler1_calc6.VAL\t \t ',
            ],
             401),
            (wdir + "/test_data/SigScan.25702", "ALS", # Sample Data for ALS6312
            [
                'Rounded Energy / eV', # index
                'Energy','Channeltron','Izero','TEY_up','TEY_dn',
            ],
             258)]
        for filename, flavour, columns, data_length in sample_data:
            test = xas._MDAdatafile(filename, flavour=flavour)
            # Assert the output is of pandas dataframe type
            self.assertEqual(type(test.dataframe), pd.core.frame.DataFrame)
            # Assert index is as expected
            self.assertEqual(test.dataframe.index.name, columns[0])
            # Assert columns are as expected
            for i in range(1, len(columns)):
                self.assertIn(columns[i], test.dataframe.columns.values.tolist()) 
            # Assert number of rows is as expected
            self.assertEqual(len(test.dataframe), data_length)

    def test_SSRL82_init(self):
        sample_data = [
            (wdir + "/test_data/", "Blank_C_tape", 886, 886) # Sample Data for SSRL82
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

    def test_athena_init(self):
        sample_data = [
            (wdir + "/test_data/FeF2.nor") # Sample .nor data
        ]
        expected_columns = [
            'Rounded Energy / eV', # index
            'energy',
            'norm',
            'bkg_norm',
            'der_norm',
            'stddev'
        ]
        for f in sample_data:
            test = xas.athena(filename=f)
            # Assert the output is of pandas dataframe type
            self.assertEqual(type(test.dataframe), pd.core.frame.DataFrame)
            # Assert index is as expected
            self.assertEqual(test.dataframe.index.name, expected_columns[0])
            # Assert columns are as expected
            for i in range(1, len(expected_columns)):
                self.assertIn(expected_columns[i], test.dataframe.columns.values.tolist()) 
            # Assert number of rows is as expected
            self.assertEqual(len(test.dataframe), 451)
