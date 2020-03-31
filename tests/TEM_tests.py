"""Unit tests for TEM.py"""

import unittest, sys, os, numpy as np, pandas as pd

wdir = os.path.dirname(__file__) # Find the current working directory
sys.path.append("..")
sys.path.append(".")
import cabanapy.TEM as tem

class TEM_init_tests(unittest.TestCase):
    """Unit test to assert filetypes are loaded correctly"""

    def test_VantageEmsaFile_init(self):
        sample_data = [
            (wdir + "/test_data/VantageEmsaFile_sample.emsa", "Sample") # Sample Data for VantageEmsaFile
        ]
        expected_columns = [
            'Energy (keV)',     # index
            'Counts'
        ]
        for filename, shortname in sample_data:
            test = tem.VantageEmsaFile(filename, shortname)
            # Assert the output is of pandas dataframe type
            self.assertEqual(type(test.dataframe), pd.core.frame.DataFrame)
            # Assert index is as expected
            self.assertEqual(test.dataframe.index.name, expected_columns[0])
            # Assert columns are as expected
            for i in range(1, len(expected_columns)):
                self.assertEqual(test.dataframe.columns.values.tolist()[i-1], expected_columns[i])
            # Assert number of rows is as expected
            self.assertEqual(len(test.dataframe), 2047)
            # Assert the metadata is as expected
            self.assertEqual(test.metadata['DATE'].strip(), '12-MAY-2017') # Date measurement was taken
            self.assertEqual(test.metadata['TIME'].strip(), '15:42: 0') # Time measurement was taken
            self.assertEqual(test.metadata['PEAKLAB'], [['3.310000', 'K', 'Ka'], # PeakLab peak identification data (if assigned during experiment)
                                                        ['6.390000', 'Fe', 'Ka'],
                                                        ['8.020000', 'Cu', 'Ka2'],
                                                        ['0.280000', 'C', 'Ka'],
                                                        ['0.670000', 'F', 'Ka'],
                                                        ['0.930000', 'Cu', 'La1'],
                                                        ['8.900000', 'Cu', 'Kb1'],
                                                        ['7.050000', 'Fe', 'Kb1'],
                                                        ['3.580000', 'K', 'Kb1'],
                                                        ['1.740000', 'Si', 'Ka']])
            self.assertEqual(test.metadata['RESULT'], [' Quantification Results', # Element weight% results (if calculated during experiment)
                                                       ' Element Weight%',
                                                       ' K   41.80 %',
                                                       ' Fe   58.20 %'])
