#-*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xml
from bruker_opus_filereader import OpusReader


# Parent Class
class _DataFile():
    """Parent class for all IR data files"""
    filename = ""
    shortname = ""
    def plot(self, color='red', legend="", invert_xaxis=True):
        wavenumber = self.dataframe['wavenumber']
        transmission = self.dataframe['transmission']
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
        
#Individual Classes
class CSVFile(_DataFile):
    """Class importing data from CSV files containing IR data"""

    def __init__(self, filename, shortname=""):
        self.filename = filename
        self.shortname = shortname

        self.dataframe = np.genfromtxt(self.filename, delimiter=',',
                                  names=['wavenumber', 'transmission'],
                                  skip_header=0, autostrip=True)

class OPUSFile(_DataFile):
    """Imports data from an OPUS Data File (.0) to dataframe """

    info_dict = {
        'DAT': 'Date',
        'FXV': 'Upper X Value',
        'LXV': 'Lower X Value',
        'DXU': 'X Units',
        'MXY': 'Max Y Value',
        'MNY': 'Min Y Value',
        'NPT': 'Number of Points',
        'TIM': 'Time',
    }
    
    def __init__(self, filename, shortname=""):
        self.filename = filename
        self.shortname = shortname
        self.metadata = {}

        sample = OpusReader(filename)
        sample.readDataBlocks()

        x_range = sample['AB Data Parameter']['LXV']-sample['AB Data Parameter']['FXV']
        step = x_range/(sample['AB Data Parameter']['NPT']-1)
        wavenumber = np.arange(sample['AB Data Parameter']['FXV'],sample['AB Data Parameter']['LXV']+step, step=step)

        dt = {'names':['wavenumber', 'transmission'], 'formats':[np.float, np.float]}
        self.dataframe = np.zeros(len(wavenumber), dtype=dt)
        self.dataframe['wavenumber'] = wavenumber
        self.dataframe['transmission'] = sample['AB']*100

        # Add metadata
        for key in sample['AB Data Parameter']:
            if key in self.info_dict:
                self.metadata[self.info_dict[key]] = sample['AB Data Parameter'][key]
