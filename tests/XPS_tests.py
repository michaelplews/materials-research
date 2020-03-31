"""Unit tests for XPS.py"""

import unittest, sys, os, numpy as np, pandas as pd

wdir = os.path.dirname(__file__) # Find the current working directory
sys.path.append("..")
sys.path.append(".")
import cabanapy.XPS as xps

class XPS_init_tests(unittest.TestCase):

    def test_KratosAsciiFile_init(self):
        sample_data = [
            (wdir + "/test_data/KratosAsciiFile_sample.txt", "Iron Acetate", "Survey") # Sample Data for KratosAsciiFile
            ]
        expected_columns = [
                 'Kinetic Energy(eV)', # index
                 'Binding Energy(eV)',
                 'Intensity(Counts)',
                 'Intensity(Counts/sec)',
                 'Transmission Value'
            ]
        for filename, shortname, scan_type in sample_data:
            test = xps.KratosAsciiFile(filename, shortname, scan_type)
            # Assert the output is of pandas dataframe type
            self.assertEqual(type(test.dataframe), pd.core.frame.DataFrame)
            # Assert header_lines has been established correctly
            self.assertEqual(test.header_lines, 1)
            # Assert index is as expected
            self.assertEqual(test.dataframe.index.name, expected_columns[0])
            # Assert columns are as expected
            for i in range(1, len(expected_columns)):
                self.assertEqual(test.dataframe.columns.values.tolist()[i-1], expected_columns[i]) 
            # Assert number of rows is as expected
            self.assertEqual(len(test.dataframe), 1101)
