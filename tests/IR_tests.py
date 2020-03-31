"""Unit tests for IR.py"""

import unittest
import sys
import os
import numpy as np
import pandas as pd

wdir = os.path.dirname(__file__) # Find the current working directory
sys.path.append("..")
sys.path.append(".")
import cabanapy.IR as ir

class IR_init_test(unittest.TestCase):
    """
    Tests for loader objects in IR.py 
    """

    def test_OPUSFile_init(self):
        sample_data = [
            (wdir + "/test_data/OPUSFile_sample.0", "Test Data") # Sample OPUS file
        ]
        expected_columns = [
            'wavelength',
            'transmission'
        ]
        for filename, shortname in sample_data:
            test = ir.OPUSFile(filename, shortname)
            # Assert output is of np.array
            self.assertEqual(type(test.dataframe), np.ndarray)
            # Assert number of rows is as expected
            self.assertEqual(len(test.dataframe), 2542)
