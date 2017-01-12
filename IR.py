#-*- coding: utf-8 -*-

import numpy as np, pandas as pd, matplotlib.pyplot as plt

# Parent Class
class _DataFile():
    """Parent class for all IR data files"""
    filename = ""
    shortname = ""
    def plot(self, color='red', legend="", invert_xaxis=True):
        wavenumber = self.dataframe()['wavenumber']
        transmission = self.dataframe()['transmission']
        if legend:
            label=legend
        elif self.shortname:
            label=self.shortname
        else:
            label="no_label"
        plt.plot(wavenumber, transmission, linewidth = 1, label=label, color=color)
        if legend is not None:
            plt.legend(loc = 1, frameon=False)
        plt.subplots_adjust(hspace = 0,wspace = 0)
        ax = plt.gca()
        ax.invert_xaxis()
        
#Individual Functions
class CSVFile(_DataFile):
    """Class importing data from CSV files containing IR data"""

    def __init__(self, filename, shortname=""):
        self.filename = filename
        self.shortname = shortname

    def dataframe(self):
        dataframe = np.genfromtxt(self.filename, delimiter=',', names=['wavenumber', 'transmission'], skip_header=0, autostrip=True)
        return dataframe
        # return dataframe['wavenumber'], dataframe['transmission']
