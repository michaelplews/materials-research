#-*- coding: utf-8 -*-
#XRD technique
import csv, os, math, xml.etree.ElementTree as ET, numpy as np
import json
from pprint import pprint
import dm3_lib as tem
from xml.etree import ElementTree
import zipfile
import numpy as np, pandas as pd, matplotlib.pyplot as plt

# Parent Classes
class _DataFile():
        """Parent class for all XRD data files"""
        filename = ""
        shortname = ""
        
        def _new_source(self, old_wavelength, new_wavelength, between):
                two_theta, intensity = self.norm_dataframe(between)
                two_theta[:] = [2*math.degrees(np.arcsin((new_wavelength*np.sin(math.radians(x/2)))/old_wavelength)) for x in two_theta]
                return two_theta, intensity

        def _normalize(self):
                """
                Internal function. Normalize the data between 0 and 100.
                """
                self.dataframe['norm_intensity'] = self.dataframe['intensity']
                self.dataframe['norm_intensity'] -= self.dataframe['norm_intensity'].min()
                self.dataframe['norm_intensity'] /= self.dataframe['norm_intensity'].max() * 0.01

        def plot(self, normalized=True, color="red", legend="", new_source=False, old_wavelength="", new_wavelength="", between=[], style=""):
                if style:
                        style.use=(style)
                if new_source:
                        two_theta, intensity = self._new_source(old_wavelength, new_wavelength, between)               
                elif normalized:
                        two_theta, intensity = self.norm_dataframe(between)
                else:
                        two_theta, intensity = self.dataframe
                if legend:
                        plt.plot(two_theta, intensity, linewidth = 1, label=legend, color=color)
                elif self.shortname:
                        plt.plot(two_theta, intensity, linewidth = 1, label=self.shortname, color=color)
                        plt.legend(loc = 1).draggable(True)                             # fontsize and other parameters removed as we now pickle our pgf styles
                        plt.subplots_adjust(hspace=0, wspace=0)
                elif legend is None:
                        plt.plot(two_theta, intensity, linewidth = 1, label=None, color=color)                        
                else:
                        plt.plot(two_theta, intensity, linewidth = 1, label='no_label', color=color)
                        plt.legend(loc = 1).draggable(True)                             # fontsize and other parameters removed as we now pickle our pgf styles
                plt.subplots_adjust(hspace=0, wspace=0)
                        
        
class _ReferenceFile():
        """Parent Class for all XRD Reference Pattern Files"""
        filename = ""
        shortname = ""
        legend = ""
        xtal_system = ""
        flavour = ""
        reference = ""
        
        def plot(self, color="red", legend="", xtal=False, hkl=False, lim=80, style=""):
                if style:
                        plt.style.use(style)
                x, y, hkl_values, h, k, l, d = self.peak_data
                if legend:
                        plt.vlines(x,-10, y, label=legend,colors=color) 
                elif self.legend:
                        plt.vlines(x,-10, y, label=self.legend,colors=color) 
                elif self.shortname:
                        plt.vlines(x,-10, y, label=self.shortname,colors=color) 
                else:
                        plt.vlines(x,-10, y, label='no_label',colors=color) 
                if xtal and self.xtal_system:
                        plt.annotate(self.xtal_system, xy=(0.01, 0.95), xycoords='axes fraction', fontsize=10,
                                horizontalalignment='left', verticalalignment='top')
                plt.legend(loc = 1).draggable(True)                             # fontsize and other parameters removed as we now pickle our pgf styles
                if hkl and hkl_values:
                        for index in range (0, len(hkl_values)):
                                if float(x[index]) < lim:
                                        plt.annotate(hkl_values[index], xy=(x[index], y[index]), xycoords='data', fontsize=10, horizontalalignment='left', verticalalignment='bottom')
                plt.xlim(float(10),float(lim))
                plt.ylim(float(0),float(110))


# Object classes to open and process data files from different sources and formats
class BM11CSVfile(_DataFile):
        '''Loads CSV files exported from 11-BM .mda Mail-In Service'''
        run_no = 0
        scan_time = ""
        proposal = ""
        temp = ""
        barcode = ""
        header_lines = 0
        dataframe = ""
        column_index = []

        def __init__(self, filename, shortname=""):
                self.filename = filename
                self.shortname = shortname
                file = open(filename, 'r')
                self.column_index = []
                i = 0
                for line in file:
                        i += 1
                        if '2theta, intensity' in line:
                                self.header_lines = i
                                self.column_index = line.strip('##').split(',')                                         
                self.dataframe = pd.read_csv(filename, sep=',', skiprows=self.header_lines, header=None, names=self.column_index, index_col=0)

        def max_in_range(self, x, y, low, high):
                """Finds the maximum value of y in a given range of x"""
                data = np.vstack((x,y))         
                y_values = data[1][np.logical_and(low < data[0], data[0] < high)]
                x_values = data[0][np.logical_and(low < data[0], data[0] < high)]
                index_max_y = y_values.argmax()
                max_y = y_values[index_max_y]
                max_x = x_values[index_max_y]
                return max_x, max_y

                        
        def norm_dataframe(self, between=[]):
                intensity = self.dataframe[' intensity'].copy()
                two_theta = intensity.index.values.copy()
                intensity = intensity.as_matrix().copy()
                if between == []:
                        max = np.amax(intensity)
                else:
                        max_x, max = self.max_in_range(two_theta, intensity, between[0], between[1])
                min = np.amin(intensity)
                a = max-min
                intensity[:] = [((x-min)/a) * 100 for x in intensity]
                return two_theta, intensity

class BrukerBrmlFile(_DataFile):
        # Taken (with permission) from https://github.com/m3wolf/scimap and edited. Thanks Mark!
        '''Loads data from a Bruker .brml v4 file to an object'''

        def __init__(self, filename, shortname=""):
                self.filename = filename
                self.shortname = shortname
                with zipfile.ZipFile(filename) as zf:
                        dataFile = zf.open('Experiment0/RawData0.xml')
                        self._dataTree = ElementTree.parse(dataFile)
                        
        @property
        def sample_name(self):
                nameElement = self._dataTree.find('.//InfoItem[@Name="SampleName"]')
                name = nameElement.get('Value')
                return name

        @property
        def dataframe(self): #Used by other functions
                index = []
                countsList = []
                # Find all Datum entries in data tree
                data = self._dataTree.findall('.//Datum')
                for datum in data:
                        time, num, two_theta, theta, counts = datum.text.split(',')
                        index.append(float(two_theta))
                        countsList.append(int(counts))
                return index, countsList
                # # Build pandas DataFrame
                # df = pd.DataFrame(countsList, index=index, columns=['counts'])
                # df.index.name = 'two_theta'
                # return df 
        
        def norm_dataframe(self, between):
                two_theta, intensity = self.dataframe
                max = np.amax(intensity)
                min = np.amin(intensity)
                a = max-min
                intensity[:] = [((x-min)/a) * 100 for x in intensity]
                return two_theta, intensity
                
class XYFile(_DataFile):
        '''Class that imports data from ASCII .xy file to an object'''
        
        def __init__(self, filename, shortname="", sep='\s+', **kwargs):
                self.filename = filename
                self.shortname = shortname
                self.dataframe = pd.read_csv(self.filename, header=0, index_col=0, names=['two_theta', 'intensity'], sep=sep, **kwargs)
                self._normalize()

        # @property
        # def dataframe(self): #Used by other functions
        #         dataframe = np.genfromtxt(self.filename, delimiter=None, names=['two_theta', 'intensity'],skip_header=0, autostrip=True)
        #         return dataframe['two_theta'], dataframe['intensity']
                
        def plot(self, signal='norm_intensity', color='r', legend=''):
                signal = self.dataframe[signal].values
                two_theta = self.dataframe.index.values
                if legend:
                        plt.plot(two_theta, signal, linewidth = 1, label=legend, color=color)
                        plt.legend(loc = 2, frameon = False).draggable(True)
                elif self.shortname: 
                        if legend is None:
                                plt.plot(two_theta, signal, linewidth = 1, color=color)
                        else:
                                plt.plot(two_theta, signal, linewidth = 1, label=self.shortname, color=color)
                                plt.legend(loc = 2, frameon = False).draggable(True)
                else:
                        plt.plot(two_theta, signal, linewidth = 1, label='no_label', color=color)
                        plt.legend(loc = 2, frameon = False).draggable(True)
                        plt.axis([np.amin(two_theta), np.amax(two_theta), 0, np.amax(signal)*1.1])
                plt.subplots_adjust(hspace=0, wspace=0)
                        
        def norm_dataframe(self):
                two_theta = self.dataframe.index.values
                intensity = self.dataframe['intensity'].values
                max = np.amax(intensity)
                min = np.amin(intensity)
                a = max-min
                intensity[:] = [((x-min)/a) * 100 for x in intensity]
                return two_theta, intensity
                
class ICDDXmlFile(_ReferenceFile):
        """
        Class that imports ICDD xml file containing XRD reference data. 'flavour' ("thousand" or "hundred") optional argument gives the option to not divide all values by 10
        """
        pdf_number = ""
        dataframe = ""

        def __init__(self, filename, flavour="thousand"):
                self.filename = filename
                tree = ET.parse(filename)
                root = tree.getroot()
                formula = root.find('.//chemical_formula')
                self.shortname = formula.text
                try:
                        self.xtal_system = xtal.text
                except:
                        pass
                pdf = root.find('.//pdf_number')
                self.pdf_number = pdf.text
                legend = filename[:-4].split(r'/')
                self.legend = legend[-1]
                self.flavour = flavour
                try:
                        self.reference = root.find('.//references/reference_group/reference').text
                except:
                        pass
                # self.dataframe = 
                
        
        @property               # Now Legacy
        def peak_data(self):#Used by other functions
                tree = ET.parse(self.filename)
                root = tree.getroot()

                theta_list, intensity_list, h_list, k_list, l_list, hkl_list, d_list = [],[],[],[],[],[],[]

                for theta in root.findall('.//theta'):
                        theta_list.append(float(theta.text))

                for intensity in root.findall('.//intensity/intensity'):
                        intensity_list.append(intensity.text)
                
                for h in root.findall('.//intensity/h'):
                        h_list.append(h.text)
                
                for k in root.findall('.//intensity/k'):
                        k_list.append(k.text)
                        
                for l in root.findall('.//intensity/l'):
                        l_list.append(l.text)
                
                for h in range (0, len(h_list)):
                        try:
                                hkl_list.append(str(h_list[h] + k_list[h] + l_list[h]))
                        except:
                                pass

                for d in root.findall('.//intensity/da'):
                        d_list.append(float(d.text))
                
                intensity_list[:] = [x for x in intensity_list if x != '\n']
                intensity_list[:] = [x.rstrip('m') for x in intensity_list]
                intensity_list[:] = [x.lstrip('<') for x in intensity_list]
                intensity_list[:] = [x or 1 for x in intensity_list]
                if self.flavour == "thousand":
                        intensity_list[:] = [int(x)/10 for x in intensity_list]
                elif self.flavour == "hundred":
                        intensity_list[:] = [int(x) for x in intensity_list]
                intensity_list[:] = [x if x > 1 else 0 for x in intensity_list]
                
                return theta_list, intensity_list, hkl_list, h_list, k_list, l_list, d_list
        
        def bragg_law(self, d_list, wavelength):
                """Returns a list of new 2theta values given a list of d_values and a wavelength via Braggs law"""
                new_twotheta = []
                for d in d_list:
                        new_twotheta.append(2*math.degrees(np.arcsin(wavelength/(2*d))))
                return new_twotheta

        def plot_wavelength(self, wavelength, color="red", legend="", xtal=False, hkl=False, lim=80):
                x, y, hkl_values, h, k, l, d = self.peak_data

                x = self.bragg_law(d, wavelength)       

                if legend is "None":
                        plt.vlines(x,-10, y ,colors=color)
                elif legend:
                        plt.vlines(x,-10, y, label=legend,colors=color)
                elif self.legend:
                        plt.vlines(x,-10, y, label=self.legend,colors=color) 
                elif self.shortname:
                        plt.vlines(x,-10, y, label=self.shortname,colors=color) 
                else:
                        plt.vlines(x,-10, y, label='no_label',colors=color) 
                if xtal and self.xtal_system:
                        plt.annotate(self.xtal_system, xy=(0.01, 0.95), xycoords='axes fraction', fontsize=10,
                                horizontalalignment='left', verticalalignment='top')
                plt.legend(loc = 1).draggable(True)                             # fontsize and other parameters removed as we now pickle our pgf styles
                if hkl and hkl_values:
                        for index in range (0, len(hkl_values)):
                                if float(x[index]) < lim:
                                        plt.annotate(hkl_values[index], xy=(x[index], y[index]), xycoords='data', fontsize=10, horizontalalignment='left', verticalalignment='bottom')
                #plt.xlim(float(10),float(lim))
                #plt.ylim(float(0),float(110))

        def export_csv(self, export_to):
                # if export_to[-3:] != "csv":
                #         export_to = export_to + '.csv'
                x, y = self.peak_data[0:2]
                _export_csv(x, y, export_to)
                
class MaterProjJSON(_ReferenceFile):
        """Class that imports downloaded XRD JSON files from materialsproject.org for comparison against other XRD patterns"""
        mp_number = ""
        json_read = ""
        amplitude, hkl, two_theta, d_spacing = [], [], [], []
        dataframe = ""
        rad_source = ""

        def __init__(self, filename, shortname=""):
                self.filename = filename
                self.shortname = shortname
                with open(filename) as json_file:
                        self.json_read = json.load(json_file)
                        self.rad_source = self.json_read['wavelength']['element']
                self.amplitude, self.hkl, self.two_theta, self.d_spacing = [], [], [], []
                for lines in self.json_read['pattern'][:]:
                        self.amplitude.append(lines[0])
                        self.hkl.append(str(lines[1][0]) + str(lines[1][1]) + str(lines[1][2]))
                        self.two_theta.append(lines[2])
                        self.d_spacing.append(lines[3])
                json_file.close()
                self.legend = self.mp_number    

        def export_csv(self, export_to):
                _export_to(self.two_theta, self.amplitude, export_to)
                
# Internal functions
def _export_csv(x, y, export_to):
        """Exports x and y data to a comma-separated ascii file.

        Arguments
        ---------
        x : list
            list of x values to be written to the first column in the .csv file
        y : list
            list of y values to be written to the second column in the .csv file
        export_to : str
            destination path to file
        """

        with open(export_to, 'w', newline='') as e:
                writer = csv.writer(e, delimiter=',')
                for i in range (0, len(x)):
                        writer.writerow([x[i], y[i]])
        
# Legacy functions
def retro_plot_xrd(filename, color="red", legend="", number_of_stacked=0):
        #normalize_xrd(filename)
        data = np.genfromtxt(filename + "_normalized.txt", delimiter='\t', names=['x', 'y'],skip_header=1, autostrip=True)
        data['y'] += 10 * number_of_stacked
        if legend:
                plt.plot(data['x'], data['y'], linewidth = 1, label=legend, color=color)
        else:
                plt.plot(data['x'], data['y'], linewidth = 1, label='no_label', color=color)
        plt.axis([10, 80, -10, 110 + number_of_stacked*10])
        plt.legend(loc = 1).draggable(True)                             # fontsize and other parameters removed as we now pickle our pgf styles
        plt.subplots_adjust(hspace=0)   
        
def plot_xrd(object, color="red", legend="", number_of_stacked=0):
        #normalize_xrd(filename)
        data = np.genfromtxt(object.filename + "_normalized.txt", delimiter='\t', names=['x', 'y'],skip_header=1, autostrip=True)
        data['y'] += 10 * number_of_stacked
        if legend:
                plt.plot(data['x'], data['y'], linewidth = 1, label=legend, color=color)
        elif object.shortname:
                plt.plot(data['x'], data['y'], linewidth = 1, label=object.shortname, color=color)
        else:
                plt.plot(data['x'], data['y'], linewidth = 1, label='no_label', color=color)
        plt.axis([10, 80, -10, 110 + number_of_stacked*10])
        plt.legend(loc = 1).draggable(True)                             # fontsize and other parameters removed as we now pickle our pgf styles
        plt.subplots_adjust(hspace=0)

