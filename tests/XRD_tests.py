"""Unit tests for XRD.py"""

import unittest, sys, os, numpy as np, pandas as pd

wdir = os.path.dirname(__file__) # Find the current working directory
sys.path.append("..")
sys.path.append(".")
import cabanapy.XRD as xrd

class XRD_functions_test(unittest.TestCase):
    """Tests for the generic functions in XRD.py. These functions are either used by other functions or used directly"""

    def test_export_csv(self):
        sample_data = [
            (np.arange(0, 99, 1), np.arange(200, 101, -1), './test_data.csv') # Sample data creates sample file
        ]
        for x, y, export_to in sample_data:
            xrd._export_csv(x, y, export_to)
            get_data = np.genfromtxt('./test_data.csv', delimiter=',', names=['x', 'y'], autostrip=True)
            # Pick a data point
            i = 25
            # Assert the x data is as expected
            self.assertEqual(25, get_data['x'][i])
            # Assert the y data is as expected
            self.assertEqual(175, get_data['y'][i])
            # Destroy the test file
            os.remove('./test_data.csv')
            
# class XRD_init_tests(unittest.TestCase):
#     """Tests to assert filetypes are initialized loaded correctly"""

#     def test_ICDDXmlFile_init(self):
#         sample_data = [
#             wdir + "/test_data/ICDDXmlFile_sample.xml", flavour="thousand")
#         ]
#         expected_columns 


# class XPS_init_tests(unittest.TestCase):

        
        
#     def test_Kratos_init(self):
#         sample_data = [
#             (wdir + "/test_data/Feacetate_survey.txt", "Iron Acetate", "Survey") # Sample Data for KratosAsciiFile
#             ]
#         expected_columns = [
#                  'Kinetic Energy(eV)',
#                  'Binding Energy(eV)',
#                  'Intensity(Counts)',
#                  'Intensity(Counts/sec)',
#                  'Transmission Value'
#             ]
#         for filename, shortname, scan_type in sample_data:
#             test = xps.KratosAsciiFile(filename, shortname, scan_type)
#             # Assert the output is of pandas dataframe type
#             self.assertEqual(type(test.dataframe), pd.core.frame.DataFrame)
#             # Assert header_lines has been established correctly
#             self.assertEqual(test.header_lines, 1)
#             # Assert index is as expected
#             self.assertEqual(test.dataframe.index.name, expected_columns[0])
#             # Assert columns are as expected
#             for i in range(1, len(expected_columns)):
#                 self.assertEqual(test.dataframe.columns.values.tolist()[i-1], expected_columns[i]) 
#             # Assert number of rows is as expected
#             self.assertEqual(len(test.dataframe), 1101)
