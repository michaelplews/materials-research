#TGA Technique 
import numpy as np, matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd

class _DataFile(object):
    def plot(self, x='Temperature (oC)', y='Mass (mg)', color='blue'):
        """Plots x vs y based on self.signal columns. Further calculation based plotting can be done using other commands."""
        plt.plot(self.dataframe[x].values, self.dataframe[y].values, color=color)
        plt.xlabel(self.dataframe[x].name, fontsize=12)
        plt.ylabel(self.dataframe[y].name, fontsize=12)
    
class TGAFile(object):
    """Loads data to an object from a .txt file created in Universal Analysis 2000 by TA Instruments. Encoding options on export from UA software must be ANSI"""
    filename = ""
    original_filename = ""
    date = ""
    time = ""
    sample = ""
    mass = ""
    method = ""
    comment = ""
    signal = []
    data_array = ""
    color = "blue"
    
    def __init__(self, filename, shortname=""):
        self.filename = filename
        self.shortname = shortname
        file = open(filename, 'r', encoding='latin-1')
        i = 0
        skip_line = 0
        self.signal = []
        for line in file:
            i += 1
            if "OrgFile" in line:
                self.original_filename = line.strip("OrgFile" + "\t").strip('\n').replace('\t', ' ').replace('\\', "/")
                skip_line = i
            elif "Date" in line:
                self.date = line.strip("Date\t").strip('\n').replace('\t', ' ')
                #self.date = datetime.strptime(line.strip('Date\t'), '%d-%b-%y')
            elif "Time" in line and "Sig" not in line:
                self.time = line.strip("Time\t").strip('\n').replace('\t', ' ')
                #self.time = datetime.strptime(line.strip('Time\t'), '%I:%M:%S')
            elif "Sample" in line and "Sig" not in line:
                self.sample = line.strip("Sample\t").strip('\n').replace('\t', ' ')
            elif "Size" in line:
                self.mass = line.strip("Size\t").strip('\n').replace('\t', ' ')
            elif "Method" in line:
                self.method = line.strip("Method\t").strip('\n').replace('\t', ' ')
            elif "Comment" in line:
                self.comment = line.strip("Comment\t").strip('\n').replace('\t', ' ')
            elif "Sig" in line:
                signal_entry = line.split('\t')
                self.signal.append(signal_entry[-1].strip('\n'))
            elif "StartOfData" in line:
                skip_line = i
        self.data_array = np.genfromtxt(filename, delimiter=None, skip_header=skip_line, autostrip=True, unpack=True)
    
    @property
    def all(self):
        print ('Filename: ' + self.filename)
        print ('Original filename: ' + self.original_filename)
        print ('Date: ' + self.date)
        print ('Time: ' + self.time)
        print ('Sample: ' + self.sample)
        print ('Mass: ' + self.mass)
        print ('Method: ' + self.method)
        print ('Comment: ' + self.comment)
        av_sigs = " | ".join(self.signal)
        print ('Available Signals :\n' + av_sigs)
        print ('Data: ' + str(len(self.data_array[0])) + " data points")
        print ('Label colour: ' + self.color)
    
    def plot(self, x=0, y=1, color='blue'):
        """Plots x vs y based on self.signal columns. Further calculation based plotting can be done using other commands."""
        plt.plot(self.data_array[x], self.data_array[y], color=color)
        plt.xlabel(self.signal[x], fontsize=12)
        plt.ylabel(self.signal[y], fontsize=12)
    
    def yforx(self, x, xdata, ydata): #Used by other functions
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
        
    def weight_percent(self, y): #Used by other functions
        """Calculates weight percent. Used by other functions. Used by other functions"""
        percent = self.data_array[y]
        max = np.amax(self.data_array[y])
        percent[:] = [x/max*100 for x in percent]
        return percent
    
    def plot_percent(self, x=1, y=2, color=""):     
        """Calculates percent mass based on column indices provided as arguments. Defaults are x=1, y=2"""
        percent = self.weight_percent(y)
        if color:
            plt.plot(self.data_array[x], percent, color=color, linewidth=1.5, label=self.shortname)
        else:
            plt.plot(self.data_array[x], percent, color=self.color, linewidth=1.5, label=self.shortname)
        plt.ylabel(r'Weight Percent / %', fontsize=12)
        plt.xlabel(self.signal[x], fontsize=12)
        legend = plt.legend(loc=1, frameon = 1, fontsize=15, framealpha=1)
        if legend:
            frame = legend.get_frame()
            frame.set_color('white')
        
    def plot_step(self, a, x=1, y=2, color="", line=True):
        """Plots point on curve of object"""
        fa, spl = self.yforx(a, self.data_array[x], self.data_array[y])
        plt.plot(a,fa,'om', color='red')
        if line:
            if color:
                plt.hlines(y=fa,xmin=a,xmax=np.amax(self.data_array[x])*10,color=color, linestyle='dashed')
            else:
                plt.hlines(y=fa,xmin=a,xmax=np.amax(self.data_array[x])*10,color=self.color, linestyle='dashed')
            
    def calculate_step(self, a, b, x=1, y=2):
        """Calculates the y difference between two x points as a float"""
        fa, spl = self.yforx(a, self.data_array[x], self.data_array[y])
        fb = interpolate.splev(b,spl,der=0)
        diff = abs(fa - fb)
        return float(diff)
        
    def gas_change(self, location, gas1="", gas2="", fontsize=15):
        plt.axvline(x=location, color='black', linestyle='dashed')
        plt.annotate(gas1, xy=(location-2, -0), xycoords='data', fontsize=fontsize,
            horizontalalignment='right', verticalalignment='top')
        plt.annotate(gas2, xy=(location+2, -0), xycoords='data', fontsize=fontsize,
            horizontalalignment='left', verticalalignment='top')
    

class KPFile(_DataFile):
    """Loads datafiles from TGA data obtained at Poppelmeier facilities, Northwestern University"""
    filename = ""
    original_filename = ""
    date = ""
    time = ""
    sample = ""
    mass = ""
    method = ""
    comment = ""
    signal = []
    data_array = ""
    color = "blue"

    def __init__(self, filename, shortname=""):
        self.filename = filename
        self.shortname = shortname
        file = open(filename, 'r', encoding='latin-1')
        i = 0
        skip_line = 0
        self.signal = []
        # self.data_array = pd(filename, delimiter='\s', skip_header=1, autostrip=True, unpack=True)
        self.dataframe = pd.read_csv(filename, encoding='utf-16', sep='\t', skiprows=0, header=0, index_col=0)
