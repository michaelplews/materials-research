# -*- coding: utf-8 -*-
#Classes and functions of XPS experiments
"""XPS data analysis for materials synthesized"""

import datetime, pandas as pd, numpy as np, matplotlib.pyplot as plt

# Parent Classes
class _DataFile():
    """Parent class for all XPS data files.

    Arguments
    ---------
    filename : str
        Name of the file you wish to import.
    shortname : str
        The name used in legend plotting and other identifying information (e.g. log).
    scan_type : str
        Type of scan (e.g. survey, O K-edge, etc.).
    header_lines : int
        Number of lines to skip when importing data from a file.
    """
    
    filename = ""
    shortname = ""
    scan_type = ""
    header_lines = 0

    def plot(self, signal='Intensity(Counts/sec)', color="black", legend='', scan_type=True):
        signal_data = self.dataframe[signal].values
        energy_data = self.dataframe['Binding Energy(eV)'].values
        if legend:
            plt.plot(energy_data, signal_data, linewidth = 1, label=legend, color=color)
            plt.legend(loc = 1, frameon = False).draggable(True)
        elif self.shortname: 
            if legend is None:
                plt.plot(energy_data, signal_data, linewidth = 1, color=color)
            else:
                plt.plot(energy_data, signal_data, linewidth = 1, label=self.shortname, color=color)
                plt.legend(loc = 1, frameon = False).draggable(True)
        else:
            plt.plot(energy_data, signal_data, linewidth = 1, label='no_label', color=color)
            plt.legend(loc = 1, frameon = False).draggable(True)
        if scan_type:
            plt.text(x=0.02, y=0.98, s=self.scan_type, horizontalalignment='left', verticalalignment='top', transform=plt.gca().transAxes)
        plt.ylabel(signal)
        plt.xlabel(self.dataframe['Binding Energy(eV)'].name)
        plt.axis([np.amin(energy_data), np.amax(energy_data), 0, np.amax(signal_data)*1.1])
        plt.gca().invert_xaxis()
        
class KratosAsciiFile(_DataFile):
    """Kratos Ascii file class opening .ascii files from Kratos AXIS-165 Surface Analysis System at UIC RRC facilities. The technique is surface sensitive (less that 8nm for XPS) with a spatial resolution in X and Y of down to 30µm.

    Arguments
    ---------
    filename : str
        Name of the file you wish to import.
    shortname : str
        The name used in legend plotting and other identifying information (e.g. log).
    scan_type : str
        Type of scan (e.g. survey, O 1s, etc.).
    header_lines : int
        Number of lines to skip when importing data from a file.
    """

    def __init__(self, filename, shortname, scan_type):
        self._log = []          # Empty the log at object initialization
        self.filename = filename
        self.shortname = shortname
        self.scan_type = scan_type
        self._AddLog('Object Initialized')

        # Open the file and read its contents
        file = open(filename, 'r')
        i=0
        for line in file:
            i+=1
            if 'Spectra ASCII data for data set...' in line:
                self.header_lines = i
        self.dataframe = pd.read_csv(filename, sep='\t', skiprows=self.header_lines, header=self.header_lines-1, index_col=0)
        file.close()
        self._AddLog('Dataframe created')

    def _AddLog(self, message):
        """Writes to object.log property to allow the user to observe and changes that have been made to the data since the object was initialised."""
        self._log.append('{time}:\t{message}'.format(time=datetime.datetime.now(), message=message))

    @property
    def log(self):
        """Prints formatted log output"""
        for line in self._log:
            print (line)
