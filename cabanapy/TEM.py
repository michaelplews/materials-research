# -*- coding: utf-8 -*-
#TEM technique
import csv, os, math, xml.etree.ElementTree as ET, numpy as np, pandas as pd
from xml.etree import ElementTree
import zipfile
import numpy as np, pandas as pd, matplotlib.pyplot as plt, dm3_lib as dm3
from matplotlib.colors import Normalize
import matplotlib.patches as patches
from scipy import interpolate


#Adapted from the dm3_lib package by Greg Jefferis from https://bitbucket.org/piraynal/pydm3reader/get/b7500989b83a.zip
class DM3File(object):
    """With help from the dm3_lib package, this object allows plotting of TEM images with appropriate nm scale"""
    filename = ""
    shortname = ""
    image = ""
    
    def __init__(self, filename, shortname=""):
        self.filename = filename
        self.shortname = shortname
    
    @property
    def outputcharset(self):
        return dm3.DM3(self.filename).outputcharset

    @property
    def tags(self):
        """Returns all image Tags."""
        return dm3.DM3(self.filename).tags

    @property
    def info(self):
        return dm3.DM3(self.filename).info


    @property
    def thumbnail(self):
        """Returns thumbnail as PIL Image."""
        return dm3.DM3(self.filename).thumbnail

    @property
    def thumbnaildata(self):
        """Returns thumbnail data as numpy.array"""
        return dm3.DM3(self.filename).thumbnaildata

    def makePNGThumbnail(self, tn_file=''):
        """Save thumbnail as PNG file."""
        return dm3.DM3(self.filename).makePNGThumbnail(tn_file=tn_file)


    @property
    def image(self):
        """Extracts image data as PIL Image"""
        return dm3.DM3(self.filename).image

    @property
    def imagedata(self):
        """Extracts image data as numpy.array"""
        return dm3.DM3(self.filename).imagedata

    @property
    def contrastlimits(self):
        """Returns display range (cuts)."""
        return dm3.DM3(self.filename).contrastlimits

    @property
    def cuts(self):
        """Returns display range (cuts)."""
        return dm3.DM3(self.filename).cuts


    @property
    def pxsize(self):
        """Returns pixel size and unit."""
        return dm3.DM3(self.filename).pxsize
        
    def plot(self, cmap="gray", scale_bar="", sb_loc = "bl", sb_color = "w", *args, **kwargs):
        """Plots the DM3 object as a plt.

        scale_bar argument accepts integer (in nanometers) e.g. scale_bar = 20 will add a white rectangle to the bottom left of the image that is 20nm wide"""
        ax=plt.gca()
        axis_size = 2048*self.pxsize[0]
        ax.imshow(self.imagedata, cmap=cmap, extent=[0, axis_size, axis_size, 0], norm=Normalize(self.contrastlimits[0],self.contrastlimits[1]), *args, **kwargs)
#           scale bar location dict
        if scale_bar:
            loc = {
                'bl': [axis_size*.05, axis_size*.90],
                'br': [(axis_size*.95)-scale_bar, axis_size*.90],
                'tr': [(axis_size*.95)-scale_bar, axis_size*.05],
                'tl': [axis_size*.05, axis_size*.05]
            }
            add_scale_bar = patches.Rectangle((50, 50), scale_bar, 5, linewidth=1, edgecolor=sb_color, facecolor=sb_color)
            add_scale_bar.set_bounds(loc[sb_loc][0],loc[sb_loc][1], scale_bar, axis_size*0.05)
            ax.add_patch(add_scale_bar)
        ax.set_ylabel(self.pxsize[1].decode('ascii'), fontsize=10)
        plt.subplots_adjust(hspace=0, wspace=0)


# EDX Data Classes
class VantageEmsaFile(object):
    """Vantage Emsa File class opens .emsa files created from VANTAGE 2.4 program used in conjunction with the JEOL 3010 EDX data collection

    Arguments
    ---------
    filename : str
        Name of the file you wish to import.
    shortname : str
        The name used in legend plotting and other identifying information (e.g. log).
    """
    
    def __init__(self, filename, shortname):
        self._log = []          # Empty the log file at initialization
        self.filename = filename
        self.shortname = shortname

        self.metadata = {}      # Empty the metadata
            
        # open the file and read its contents
        file = open(self.filename, 'r')
        i=0
        p=[]
        r=[]
        for line in file:
            if '#' in line:
                i+=1
                # Assign files metadata to a string
                (key, val) = line.split(':', 1)
                key=''.join(c for c in key if c not in '#\t\n ')
                val=''.join(c for c in val if c not in '#\t\n')
                if str(key) == 'PEAKLAB':
                    p.append(val.split())
                if str(key) == 'RESULT':
                    r.append(val)
                self.metadata[str(key)] = val
        self.metadata['RESULT'] = r
        self.metadata['PEAKLAB'] = p
        self.dataframe = pd.read_csv(filename, sep=',', skiprows=i, skipfooter=2, engine='python', index_col=False, header=None, names=['Energy (keV)', 'Counts'])
        self.dataframe = self.dataframe.set_index('Energy (keV)')
        file.close()
       
    def plot(self, color="grey", show_results=False, r_loc = 'tl', label_peaks=False, legend=''):
        signal_data = self.dataframe['Counts'].values
        energy_data = self.dataframe.index.values
        plt.vlines(x=energy_data, ymax=signal_data, ymin=[0 for x in signal_data], color=color)
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
        if show_results:
            loc = {             # x, y, horizontalalignment, verticalalignment
                'bl': [0, 0, 'left', 'bottom'],
                'br': [0.995, 0, 'right', 'bottom'],
                'tr': [0.995,0.98, 'right', 'top'],
                'tl': [0, 0.98, 'left', 'top']
            }
            plt.text(x=loc[r_loc][0], y=loc[r_loc][1], s='\n'.join(self.metadata['RESULT']), horizontalalignment=loc[r_loc][2], verticalalignment=loc[r_loc][3], transform=plt.gca().transAxes)
        # if label_peaks:
        #     peak = self.metadata['PEAKLAB']
        #     for i in range(0, len(self.metadata['PEAKLAB'])-1):
        #         plt.annotate(
        #             xy=[peak[i][0],label_peaks[i]],
        #             s=peak[i][1]
        #         )
        if label_peaks:
            for peak in self.metadata['PEAKLAB']:
                plt.annotate(
                    xy=[float(peak[0]),
                        max_in_range(
                            self.dataframe['Counts'].values,
                            self.dataframe.index.values,
                            low=float(peak[0])-0.05,
                            high=float(peak[0])+0.05,
                            plot=False,
                            do_return=True)][1],
                    s=str(peak[1]),
                    horizontalalignment='center',
                    verticalalignment='bottom'
                )
        plt.ylabel(self.dataframe['Counts'].name)
        plt.xlabel(self.dataframe.index.name)
        plt.axis([np.amin(energy_data), np.amax(energy_data), 0, np.amax(signal_data)*1.1])

# External Functions
def max_in_range(signal, index, low, high, plot=True, do_return=False):
    """Finds the maximum value of y in a given range of x"""
    get_signal = signal
    energy = index
    data = np.vstack((energy, get_signal))
    y_values = data[1][np.logical_and(low < data[0], data[0] < high)]
    x_values = data[0][np.logical_and(low < data[0], data[0] < high)]
    index_max_y = y_values.argmax()
    max_y = y_values[index_max_y]
    max_x = x_values[index_max_y]
    if plot:
        plt.plot(max_x,max_y,'om', color='red')
    if do_return:
        return max_x, max_y

def yforx(x, xdata, ydata): #Used by other functions
    """Sorts xdata into ascending order (required by splrep) and solves y (as fa) for a value of x. Also returns spl to allow for further calculations to take place quicker. Used by other functions."""
    #Sorts data into ascending order (required by splrep)
    data = xdata
    indices = np.argsort(data)
    data_index = data[indices]
    percent_index = ydata[indices]
    #interpolates data and solves y for a specific x
    spl = interpolate.splrep(data_index, percent_index, s=1) 
    #s=1 smoothing to eradicate duplicates
    fa = interpolate.splev(x,spl,der=0)     # f(a)
    return fa, spl

