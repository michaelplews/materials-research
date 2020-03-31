# -*- coding: utf-8 -*-
#Classes and functions of XAS experiments

import datetime, numpy as np, operator, pandas as pd, matplotlib.pyplot as plt
from scipy import interpolate

# Parent Classes
class _DataFile():
    """Parent class for all XAS data files"""
    filename = ""
    shortname = ""
    dataframe = ""
    
    def plot(self, signal, color="red", legend=''):
        if self.processed_dataframe is not None:
            signal = self.processed_dataframe[signal].values
            energy = self.processed_dataframe['Energy / eV'].values
        elif self.dataframe is not None:
            signal = self.dataframe[signal].values
            energy = self.dataframe.index.values
        else:
            RaiseException('No dataframe found')
        if legend:
            plt.plot(energy, signal, linewidth = 1, label=legend, color=color)
            plt.legend(loc = 2, frameon = False).draggable(True)
        elif self.shortname: 
            if legend is None:
                plt.plot(energy, signal, linewidth = 1, color=color)
            else:
                plt.plot(energy, signal, linewidth = 1, label=self.shortname, color=color)
                plt.legend(loc = 2, frameon = False).draggable(True)
        else:
            plt.plot(energy, signal, linewidth = 1, label='no_label', color=color)
            plt.legend(loc = 2, frameon = False).draggable(True)
        # plt.axis([np.amin(energy), np.amax(energy), 0, np.amax(signal)*1.1])
        # plt.subplots_adjust(hspace=0, wspace=0)

    def _ScaleRef(self, column_id, name=None, divisor=None, trim="", smooth=None):
        """
        Internal function. Scaled data between 1 and 0 irrespective of wave shape. Used for STD and TEY data analysis
        """
        if name is None:
            name = column_id
        if trim:        
            new_frame = self.normalized_dataframe[column_id].iloc[trim[0]:trim[1]].copy()
        else:
            new_frame = self.normalized_dataframe[column_id].copy()
        if divisor:
            new_frame = new_frame/self.normalized_dataframe[divisor].copy()
        if smooth:
            new_frame = pd.Series.rolling(new_frame, window=smooth,center=True,min_periods=smooth).mean()       # Updated due to future warning
            #new_frame = pd.rolling_mean(new_frame, smooth, smooth, center=True)                  # Legacy
        new_frame -= new_frame.min()
        new_frame /= new_frame.max()
        new_frame.columns = ['Rounded Energy / eV', name]
        return new_frame

    def _SumData(self):
        """
        Internal function. Averages data across multiple _MDAFile objects
        """
        normalized_dataframe = []
        normalized_dataframe = pd.concat([scan.dataframe for scan in self._MDAlist]).groupby(level=0).mean()
        #         normalized_dataframe = pd.Dataframe([scan.dataframe for scan in self._MDAlist])
        return normalized_dataframe

    def _ScaleAbs(self, column_id, name, divisor=None, head=5, tail=5, smooth=None, trim=""):
        """
        Internal function. Scales data absolutely between 1 and 0 depending on the wave shape. Used for TFY data analysis.
        """
        if trim:
            new_frame = self.normalized_dataframe[column_id].iloc[trim[0]:trim[1]].copy()
        else:
            new_frame = self.normalized_dataframe[column_id].copy()
        if divisor:
            new_frame = new_frame/self.normalized_dataframe[divisor].copy()
        # if smooth is None:
        #     endloc = self.normalized_dataframe['energy'].iloc[-1]
        # else:
        endloc = self.normalized_dataframe['Rounded Energy / eV'].iloc[-1]
        v_minloc = new_frame.idxmin()
        if v_minloc < endloc/2+v_minloc/2:
            new_frame -= new_frame.head(head).mean()
            new_frame /= new_frame.tail(tail).mean()
        else:
            new_frame -= new_frame.tail(tail).mean()
            new_frame /= new_frame.head(head).mean()
        if smooth:
            new_frame = pd.Series.rolling(new_frame, window=smooth,center=True,min_periods=smooth).mean()       # Updated due to future warning
            #new_frame = pd.rolling_mean(new_frame, smooth, smooth, center=True)                  # Legacy
        new_frame.columns = ['Rounded Energy / eV', name]
        return new_frame

    def yforx(self, x, xdata, ydata): #Used by other functions
        """
        Sorts xdata into ascending order (required by splrep) and solves y (as fa) for a value of x. Also returns spl to allow for further calculations to take place quicker. Used by other functions.
        """
        #Sorts data into ascending order (required by splrep)
        data = xdata
        indices = np.argsort(data)
        data_index = data[indices]
        percent_index = ydata[indices]
        #interpolates data and solves y for a specific x
        spl = interpolate.splrep(data_index, percent_index, s=0) 
        #If errors, try s=1 smoothing to eradicate duplicates
        fa = interpolate.splev(x,spl,der=0)   # f(a)
        return fa, spl
  
    def align(self, x_from, x_to, dataframe='processed_dataframe', use_index=False):
        """
        Applies the appropriate multiplication/division to shift waves in order to account for beam drift and other sources of expermental error. Should be used on a similar peak where x_from=original x location and x_to=desired x location.
        """
        dict = {
            'dataframe': self.dataframe,
            'processed_dataframe': self.processed_dataframe
        }
        
        if use_index:
            dict[dataframe].reset_index()
        energy = dict[dataframe]['Energy / eV']
        if x_to > x_from:
            energy *= x_to/x_from
        else:
            energy /= x_from/x_to

        # if use_index:
        #     dict[dataframe].set_index('Energy / eV')

        self._AddLog(self.shortname + " aligned from " + str(x_from) + "eV to " + str(x_to) + "eV")

    def max_in_range(self, signal, low, high, plot=True, do_return=False, dataframe='processed_dataframe', use_index=False):
        """
        Finds the maximum value of y in a given range of x
        """
        dict = {
            'dataframe':self.dataframe,
            'processed_dataframe':self.processed_dataframe
        }
        
        get_signal = dict[dataframe][signal].values
        if use_index:
            energy = dict[dataframe].index.values
        else:
            energy = dict[dataframe]['Energy / eV'].values
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
        
    def _AddLog(self, message):
        """
        Writes to object.log property to allow the user to observe and changes that have been made to the data since the object was initialised.
        """
        self._log.append('{time}:\t{message}'.format(time=datetime.datetime.now(), message=message))


    @property
    def log(self):
        """Prints formatted log output"""
        for line in self._log:
            print (line)

    def fit_linear(self, signal, x1, x2):
        """
        Fits a linear regression to two defined points. Returns the equation of the straight line (m, c). 
        m = gradient, c = constant
        """
        energy = self.processed_dataframe['Energy / eV'].values
        y1, y1spl = self.yforx(x1, energy, self.processed_dataframe[signal].values)
        y2, y2spl = self.yforx(x2, energy, self.processed_dataframe[signal].values)
        # Taken from http://stackoverflow.com/a/21566184
        points = [(x1, y1),(x2, y2)]
        x_coords, y_coords = zip(*points)
        A = np.vstack([x_coords,np.ones(len(x_coords))]).T
        m, c = np.linalg.lstsq(A, y_coords)[0]
        return m, c

    def plot_linear(self, m, c, color='blue'):
        """Plots a linear equation of y=mx+c"""
        y = [m*x+c for x in self.processed_dataframe['Energy / eV'].values]
        plt.plot(self.processed_dataframe['Energy / eV'].values, y, color=color)


    def renormalize(self, signal, TFY_smooth=7, trim_tey="", trim_tfy=""):
        """Provides renormalization in case data needs to be reworked"""
        self.processed_dataframe[signal] -= self.processed_dataframe[signal].min()
        self.processed_dataframe[signal] /= self.processed_dataframe[signal].max()
                   
    def subtract_linear(self, signal, m, c):
        """Subtracts a linear equation of y=mx+c from a specified signal"""
        y = [m*x+c for x in self.processed_dataframe['Energy / eV'].values]
        self.processed_dataframe[signal] = self.processed_dataframe[signal].subtract(y)
        self._AddLog(self.shortname + " " + signal +" background subtraction of " + "y = {m}x + {c}".format(m=m, c=c))


# Sub-classes
class _MDAdatafile():
    """
    Loads .0001 data files produced at Argonne National Lab, or Advanced Light Source.
    """
    filename = ""
    basename = ""
    ext = ""
    header_lines = 0
    dataframe = ""
    column_index = []
    scan_datetime = ""
    
    def __init__(self, filename, header_lines=0, flavour="IDC4"):
        """
        filename: 
            Pretty self explanatory.
        header_lines: 
            How many lines of the file needs to be skipped before the data starts (titles will therefore be self.header_lines-1).
        flavour: 
            Just a way of identifying how the data should be extracted depending on which facility it was collected at.
        """
        self.filename = filename
        self.ext = filename[-4:]
        self.basename = filename[:-5].split(r'/')[-1]
        file = open(filename, 'r')
        self.column_index = []
        self.header_lines = header_lines
        i = 0
        if flavour == "IDC4":
            for line in file:
                i += 1
                if '# 1-D Scan Values' in line:
                    self.header_lines = i
                elif '# Scan time' in line:
                    datetime_string = line.strip('# Scan time = ').strip('\n').replace('\t', ' ')
                    self.scan_datetime = datetime.datetime.strptime(datetime_string, "%b %d, %Y %H:%M:%S.%f")
                elif '[' and ']' in line:
                    self.column_index.append('['+line.split('[',1)[-1].strip('#').strip('\n').replace(',','\t'))
            self.dataframe = pd.read_csv(filename, sep=' ', skiprows=self.header_lines, header=None, names=self.column_index, index_col=0)
            self.dataframe['Rounded Energy / eV'] = pd.Series(np.around(self.dataframe[self.column_index[1]], decimals=1), index=self.dataframe.index)
            self.dataframe.index = self.dataframe['Rounded Energy / eV']
        elif flavour == "ALS":
            for line in file:
                i += 1
                if "Time (s)" in line: # This really depends on the data titles in the datafile. Might have to find a better way to do this in the future.
                    self.header_lines = i
            self.dataframe = pd.read_csv(filename, sep='\t', skiprows=self.header_lines-1, header=0, index_col=0)
            try:
                self.dataframe['Rounded Energy / eV'] = pd.Series(np.around(self.dataframe['Energy'], decimals=1), index=self.dataframe.index)
            except:
                self.dataframe['Rounded Energy / eV'] = pd.Series(np.around(self.dataframe['Mono Energy'], decimals=1), index=self.dataframe.index)
            self.dataframe.index = self.dataframe['Rounded Energy / eV']
        elif flavour == "SSRL":
            self.dataframe = pd.read_csv(filename, sep=' ', skiprows=0, header=0, index_col=0)
            self.dataframe.columns = [s.strip(' ') for s in self.dataframe.columns.tolist()]
            self.dataframe['Rounded Energy / eV'] = pd.Series(np.around(self.dataframe.index, decimals=1), index=self.dataframe.index)
            self.dataframe['mono'] = self.dataframe.index.copy()
            self.dataframe.index = self.dataframe['Rounded Energy / eV']
             
    @property
    def available_signals(self):
        print ('Column\t|\tDetector')
        for line in range (0, len(self.column_index)):
            print (str(line) + '\t|\t' + self.column_index[line])


# Object classes to open and process datafiles from different sources and formats
class igorTXTFile(_DataFile):
    """Loads .txt wave files exported by IGOR Pro"""
    filename = ""
    dataframe = ""  

    def __init__(self, filename, shortname=""):
        self.filename = filename
        self.shortname = shortname
        self.dataframe = np.genfromtxt(filename, delimiter='\t', skip_header=1, autostrip=True, unpack=True)

class IDC4(_DataFile):
    """
    Loads .0001 data files produced by 4-ID-C at Argonne National Lab, Uses the _MDAFile class to average data scans and produce a normalized array.
    """
    directory = ""
    basename = ""
    start = ""
    end = ""
    exclude = []
    normalized_dataframe = ""
    processed_dataframe = ""
    dataframe = ""
    beamline = "4-ID-C, Advanced Photon Source"
    _log = []

    #Variables that will change if instrument set up is adjusted
    Energy_id = '[1-D Positioner 1]  4idc1:SGM1:Energy\t SGM1:Energy\t LINEAR\t eV\t 4idc1:SGM1:EnergyRBV\t Energy readback\t eV'   
    TFY_id = '[1-D Detector  10]  4idc1:scaler1_calc5.VAL\t \t '
    TEY_id = '[1-D Detector   9]  4idc1:scaler1_calc4.VAL\t \t '
    REF_id = ''
    STD_id = '[1-D Detector  11]  4idc1:scaler1_calc6.VAL\t \t '

    def __init__(self, directory, basename, start, end, exclude=None, shortname="", TFY_smooth=7, trim_tey="", trim_tfy=""):
        self._log = []  #reset the log to be empty
        self.directory = directory
        self.basename = basename
        self.start = int(start)
        self.end = int(end)
        self.shortname = shortname
        self._MDAlist = []
        if exclude and type(exclude) is list:
            self.exclude = list(map(int, exclude))
        elif exclude and type(exclude) is str:
            self.exclude = list(map(int, str(exclude)))    
        self._MDAlist = [_MDAdatafile(directory + basename + '.' + str(ext).zfill(4)) 
                        for ext in range(self.start, self.end+1)
                        if ext not in self.exclude]
        self.normalized_dataframe = self._SumData()
        self._AddLog('Object created (__init__)')
        self._AddLog('normalized_dataframe created')      

        #More dataframes can be appended if needed
        self.normalized_dataframe['TFY'] = self._ScaleAbs(self.TFY_id, 'TFY', trim=trim_tfy)
        self.normalized_dataframe['sTFY'] = self._ScaleAbs(self.TFY_id, 'TFY', smooth=TFY_smooth, trim=trim_tfy)
        self._AddLog('TFY smoothed: rolling mean window = ' + str(TFY_smooth))
        self.normalized_dataframe['TEY'] = self._ScaleRef(self.TEY_id, name='TEY', trim=trim_tey)   
        self.normalized_dataframe['STD'] = self._ScaleRef(self.STD_id, name='STD')
        self.normalized_dataframe['TFY2'] = self._ScaleRef(self.TFY_id, 'TFY', trim=trim_tfy)
        self.normalized_dataframe['sTFY2'] = self._ScaleRef(self.TFY_id, 'TFY', smooth=TFY_smooth, trim=trim_tfy)

        
        #Tidy into a seperate processed dataframe
        self.processed_dataframe = self.normalized_dataframe.ix[:,[
            self.Energy_id,
            'TFY',
            'sTFY',
            'TEY',
            'STD',
            'TFY2',
            'sTFY2',
        ]].copy()
        self.processed_dataframe.rename(columns={self.Energy_id:'Energy / eV'}, inplace=True)
        self._AddLog('processed_dataframe created')

    @property
    def info(self):
        print ('Filename:\t|\tAquisition Started:\t\t|\tNumber of Data Points:')
        for index in range(0, len(self._MDAlist)):
            try:
                print (self._MDAlist[index].basename + '.' + self._MDAlist[index].ext+'\t|\t'+str(self._MDAlist[index].scan_datetime)+'\t|\t'+str(len(self._MDAlist[index].dataframe)))
            except IndexError:
                print (self._MDAlist[index].basename + '.' + self._MDAlist[index].ext+'\t|\t'+str(self._MDAlist[index].scan_datetime)+'\t|\t'+'Empty File')
        
class DATFile(_DataFile):
    """
    So far only one case of this extension appearing. When JF sends reference data as a request.
    """
    dataframe = ""  

    def __init__(self, filename, shortname=""):
        self.filename = filename
        self.shortname = shortname
        self.dataframe = np.genfromtxt(filename, delimiter='\t', skip_header=1, autostrip=True, unpack=True)

class ALS6312(_DataFile):
    """
    Loads .txt data files produced by beamline 6312 at The Advanced Light Source, Lawrence Berkeley National Laboratory. The _MDAFile class is not needed as data was pre-averaged at the beamline before release. Before processing all files were renamed using bash scripts:

    $ for i in *-Avg.txt ; do mv "$i" "${i/-Avg.txt/.txt}" ; done
    $ for i in SigScan* ; do mv "$i" "${i%".txt"}" ; done
    $ for i in SigScan* ; do mv "$i" "${i/#"SigScan"/"SigScan."}" ; done

    This renamed the files from "SigScan22222-Avg.txt" to "SigScan.22222". This makes it easier to use the prewritten IDC4 object code.
    """
    directory = ""
    basename = ""
    start = ""
    end = ""
    exclude = []
    normalized_dataframe = ""
    processed_dataframe = ""
    beamline = "6321, Advanced Light Source"
    _log = []

    #Variables that will change if instrument set up is adjusted
    Energy_id = 'Energy'  
    TFY_id = 'Channeltron'
    TEY_id = 'TEY_up'
    I0_id = 'Izero'
    REF_id = None
    STD_id = None

    def __init__(self, directory, basename, start, end=0, exclude=None, shortname="", tey_detector="", TFY_smooth=7, trim_tey="", trim_tfy=""):
        self._log = []  #reset the log to be empty
        self.directory = directory
        self.basename = basename
        self.start = int(start)
        if end:
            self.end = int(end)
        else:
            self.end = int(start)
        self.shortname = shortname
        self._MDAlist = []
        if exclude and type(exclude) is list:
            self.exclude = list(map(int, exclude))
        elif exclude and type(exclude) is str:
            self.exclude = list(map(int, str(exclude)))
        self._MDAlist = [_MDAdatafile(directory + basename + '.' + str(ext), flavour='ALS') 
                        for ext in range(self.start, self.end+1)
                        if ext not in self.exclude]
        self.normalized_dataframe = self._SumData()
        self._AddLog('Object created (__init__)')
        self._AddLog('normalized_dataframe created')      

        #More dataframes can be appended if needed
        self.normalized_dataframe['TFY'] = self._ScaleAbs(self.TFY_id, 'TFY', divisor=self.I0_id, trim=trim_tfy)
        self.normalized_dataframe['sTFY'] = self._ScaleAbs(self.TFY_id, 'sTFY', divisor=self.I0_id, smooth=TFY_smooth, trim=trim_tfy)
        self._AddLog('TFY smoothed: rolling mean window = ' + str(TFY_smooth))
        if tey_detector:
            self.normalized_dataframe['TEY'] = self._ScaleRef(tey_detector, 'TEY', divisor=self.I0_id, trim=trim_tey)
        else:
            self.normalized_dataframe['TEY'] = self._ScaleRef(self.TEY_id, 'TEY', divisor=self.I0_id, trim=trim_tey)
        self._AddLog('TFY, sTFY and TEY divided by I0')

        #Tidy into a seperate processed dataframe
        self.processed_dataframe = self.normalized_dataframe.ix[:,[
            self.Energy_id,
            'TFY',
            'sTFY',
            'TEY'
        ]].copy()
        self.processed_dataframe.rename(columns={self.Energy_id:'Energy / eV'}, inplace=True)
        self._AddLog('processed_dataframe created')

    @property
    def info(self):
        print ('Filename:\t|\tAquisition Started:\t\t|\tNumber of Data Points:')
        for index in range(0, len(self._MDAlist)):
            try:
                print (self._MDAlist[index].basename + '.' + self._MDAlist[index].ext+'\t|\t'+str(self._MDAlist[index].scan_datetime)+'\t|\t'+str(len(self._MDAlist[index].dataframe)))
            except IndexError:
                print (self._MDAlist[index].basename + '.' + self._MDAlist[index].ext+'\t|\t'+str(self._MDAlist[index].scan_datetime)+'\t|\t'+'Empty File')
    
class ALS801(_DataFile):
    """
     Loads .txt data files produced by beamline 8.0.1 (also known as 8.0.3) at The Advanced Light Source, Lawrence Berkeley National Laboratory. The _MDAFile class is not needed as data was pre-averaged at the beamline before release. Before processing all files were renamed using bash scripts:


    $ for i in *-Avg.txt ; do mv "$i" "${i/-Avg.txt/.txt}" ; done
    $ for i in SigScan* ; do mv "$i" "${i%".txt"}" ; done
    $ for i in SigScan* ; do mv "$i" "${i/#"SigScan"/"SigScan."}" ; done

    This renamed the files from "SigScan22222-Avg.txt" to "SigScan.22222". This makes it easier to use the prewritten IDC4 object code.
    """
    directory = ""
    basename = ""
    start = ""
    end = ""
    exclude = []
    normalized_dataframe = ""
    processed_dataframe = ""
    beamline = "8.0.1, Advanced Light Source"
    _log = []

    #Variables that will change if instrument set up is adjusted
    Energy_id = 'Mono Energy'     
    TFY_id = 'Counter 2'
    TEY_id = 'Counter 1'
    I0_id = 'Counter 3'   # Can also use 'Counter 0' if the signal-to-noise ratio is better
    REF_id = None
    STD_id = None

    def __init__(self, directory, basename, start, end=0, exclude=None, shortname="", tey_detector="", TFY_smooth=7, trim_tey="", trim_tfy=""):
        self._log = []  #reset the log to be empty
        self.directory = directory
        self.basename = basename
        self.start = int(start)
        if end:
            self.end = int(end)
        else:
            self.end = int(start)
        self.shortname = shortname
        self._MDAlist = []
        if exclude and type(exclude) is list:
            self.exclude = list(map(int, exclude))
        elif exclude and type(exclude) is str:
            self.exclude = list(map(int, str(exclude)))
        self._MDAlist = [_MDAdatafile(directory + basename + '.' + str(ext), flavour='ALS') 
                        for ext in range(self.start, self.end+1)
                        if ext not in self.exclude]
        self.normalized_dataframe = self._SumData()
        self._AddLog('Object created (__init__)')
        self._AddLog('normalized_dataframe created')      

        #More dataframes can be appended if needed
        self.normalized_dataframe['TFY'] = self._ScaleAbs(self.TFY_id, 'TFY', divisor=self.I0_id, trim=trim_tfy)
        self.normalized_dataframe['sTFY'] = self._ScaleAbs(self.TFY_id, 'sTFY', divisor=self.I0_id, smooth=TFY_smooth, trim=trim_tfy)
        self._AddLog('TFY smoothed: rolling mean window = ' + str(TFY_smooth))
        if tey_detector:
            self.normalized_dataframe['TEY'] = self._ScaleRef(tey_detector, 'TEY', divisor=self.I0_id, trim=trim_tey)
        else:
            self.normalized_dataframe['TEY'] = self._ScaleRef(self.TEY_id, 'TEY', divisor=self.I0_id, trim=trim_tey)
        self._AddLog('TFY, sTFY and TEY divided by I0')

        #Tidy into a seperate processed dataframe
        self.processed_dataframe = self.normalized_dataframe.ix[:,[
            self.Energy_id,
            'TFY',
            'sTFY',
            'TEY'
        ]].copy()
        self.processed_dataframe.rename(columns={self.Energy_id:'Energy / eV'}, inplace=True)
        self._AddLog('processed_dataframe created')

    @property
    def info(self):
        print ('Filename:\t|\tAquisition Started:\t\t|\tNumber of Data Points:')
        for index in range(0, len(self._MDAlist)):
            try:
                print (self._MDAlist[index].basename + '.' + self._MDAlist[index].ext+'\t|\t'+str(self._MDAlist[index].scan_datetime)+'\t|\t'+str(len(self._MDAlist[index].dataframe)))
            except IndexError:
                print (self._MDAlist[index].basename + '.' + self._MDAlist[index].ext+'\t|\t'+str(self._MDAlist[index].scan_datetime)+'\t|\t'+'Empty File')

class SSRL82(_DataFile):
    """
    Imports data files from beamline 8.2 at SSRL, Stanford Linear Accelerator Center. 
    Resulting object contains a raw dataframe, and a processed dataframe.  
    Files should be renamed to form Name.001 where Name is a common name shared by similar files. Make copies of files before renaming, as angle and position data is stored int he filename. 
    Notes:
     - Files should be taken from spec_export/ (these have been pre-processed)
     - These signals are already divided by i0, so divisor is set to None
    
    Files should be renamed using bash:

    $ for i in *.dat ; do mv $i $(echo $i | sed -e 's/.*_//') ; done
    $ for i in *.dat ; do mv "$i" "Scan.$i" ; done
    $ for i in *.dat ; do mv "$i" "${i%".dat"}" ; done
    
    """
    directory = ""
    basename = ""
    start = ""
    end = ""
    exclude = []
    normalized_dataframe = ""
    processed_dataframe = ""
    beamline = "8.2, Stanford Linear Accelerator Center"
    _log = []
    default_trim = {
        'tey': [],
        'tfy': [],
        'aey': [],
        'pey': [],
    }

    #Variables that will change if instrument set up is adjusted
    Energy_id = 'mono'
    TFY_id = 'tfy'
    TEY_id = 'tey'
    AEY_id = 'aey'
    PEY_id = 'pey'
    I0_id = 'i0'
    REF_id = 'refy'
    STD_id = None
    
    def __init__(self, directory, basename, start, end=0, exclude=None, shortname="", tey_detector="", TFY_smooth=7, trim={}):
        self.default_trim = { # reset the trim dict
            'tey': [],
            'tfy': [],
            'aey': [],
            'pey': [],
        }
        self._log = []  #reset the log to be empty
        self.directory = directory
        self.basename = basename
        self.start = int(start)
        if end:
            self.end = int(end)
        else:
            self.end = int(start)
        self.shortname = shortname
        self._MDAlist = []
        if exclude and type(exclude) is list:
            self.exclude = list(map(int, exclude))
        elif exclude and type(exclude) is str:
            self.exclude = list(map(int, str(exclude)))
        self._MDAlist = [_MDAdatafile(directory + basename + '.' + str(ext), flavour='SSRL') 
                        for ext in range(self.start, self.end+1)
                        if ext not in self.exclude]
        self.normalized_dataframe = self._SumData()
        self._AddLog('Object created (__init__)')
        self._AddLog('normalized_dataframe created')      

        # Apply any update the trim dict with new data 
        self.default_trim.update(trim)
        
        # More dataframes can be appended if needed
        self.normalized_dataframe['REF'] = self._ScaleRef(self.REF_id, 'REF', divisor=None)
        self.normalized_dataframe['TFY'] = self._ScaleRef(self.TFY_id, 'TFY', divisor=None, trim=self.default_trim['tfy'])
        self.normalized_dataframe['sTFY'] = self._ScaleRef(self.TFY_id, 'sTFY', divisor=None, smooth=TFY_smooth, trim=self.default_trim['tfy'])
        self._AddLog('TFY smoothed: rolling mean window = ' + str(TFY_smooth))
        if tey_detector:
            self.normalized_dataframe['TEY'] = self._ScaleRef(tey_detector, 'TEY', divisor=None, trim=self.default_trim['tey'])
        else:
            self.normalized_dataframe['TEY'] = self._ScaleRef(self.TEY_id, 'TEY', divisor=None, trim=self.default_trim['tey'])
        self.normalized_dataframe['AEY'] = self._ScaleRef(self.AEY_id, 'AEY', divisor=None, trim=self.default_trim['aey'])
        self.normalized_dataframe['PEY'] = self._ScaleRef(self.PEY_id, 'PEY', divisor=None, trim=self.default_trim['pey'])
        self._AddLog('REF, TFY, sTFY, TEY ,PEY and AEY divided by I0')

        #Tidy into a seperate processed dataframe
        self.processed_dataframe = self.normalized_dataframe.ix[:,[
            self.Energy_id,
            'REF',
            'TFY',
            'sTFY',
            'TEY',
            'AEY',
            'PEY'
        ]].copy()
        self.processed_dataframe.rename(columns={self.Energy_id:'Energy / eV'}, inplace=True)
        self._AddLog('processed_dataframe created')

    @property
    def info(self):
        print ('Filename:\t|\tAquisition Started:\t\t|\tNumber of Data Points:')
        for index in range(0, len(self._MDAlist)):
            try:
                print (self._MDAlist[index].basename + '.' + self._MDAlist[index].ext+'\t|\t'+str(self._MDAlist[index].scan_datetime)+'\t|\t'+str(len(self._MDAlist[index].dataframe)))
            except IndexError:
                print (self._MDAlist[index].basename + '.' + self._MDAlist[index].ext+'\t|\t'+str(self._MDAlist[index].scan_datetime)+'\t|\t'+'Empty File')
        
class athena(_DataFile):
    """
    Files exported from Athena (part of the Horae package). Resulting object contains a dataframe that is unprocessed (as that was presumably already done in Athena). One file should be imported per object.
    """
    _log = []
    header_lines = 0
    column_list = []
    processed_dataframe = None
    
    #Variables that will change if instrument set up is adjusted
    Energy_id = 'energy'
    norm_id = 'norm'
    bkg_norm_id = 'bkg_norm'
    der_norm_id = 'der_norm'
    stddev_id = 'stddev'

    def __init__(self, filename, shortname="", signal_type='nor'):
        self._log = []  #reset the log to be empty
        self.filename = filename
        self.shortname = shortname
        f = open(filename)

        # Empty column_list
        column_list = []
        i=0
        for line in f:
            if '#---------' in line:
                self.header_lines = i+2
                self.column_list = [s for s in next(f).strip('#').strip('\n').split(' ') if s is not '']
            else:
                i+=1
        self.dataframe = pd.read_csv(filename, sep='\s+', skiprows=self.header_lines, names=self.column_list, index_col=0)
        self.dataframe['Rounded Energy / eV'] = pd.Series(np.around(self.dataframe.index, decimals=1), index=self.dataframe.index)
        self.dataframe['energy'] = self.dataframe.index.copy()
        self.dataframe.index = self.dataframe['Rounded Energy / eV']

        self._AddLog('dataframe created')

class datathief(_DataFile):
    """
    Loads .txt files exported by Datathief.jar
    """
    directory = ""
    basename = ""
    start = ""
    end = ""
    exclude = []
    normalized_dataframe = ""
    dataframe = ""
    processed_dataframe = None
    _log = []

    #Variables that will change if instrument set up is adjusted
    energy_id = 'energy'
    signal_id = 'signal'

    def __init__(self, filename, shortname="", trim="", method='ScaleRef'):
        self._log = []  #reset the log to be empty
        self.shortname = shortname
        self._MDAlist = []

        # Read the file to create the dataframe
        self.normalized_dataframe = pd.read_csv(filename, names=['energy', 'signal'], skiprows=1)

        self._AddLog('Object created (__init__)')
        self._AddLog('normalized_dataframe created')      

        dict={
            'ScaleRef':self._ScaleRef,
            'ScaleAbs':self._ScaleAbs,
        }
        
        #More dataframes can be appended if needed
        self.normalized_dataframe['signal'] = dict[method](self.signal_id, 'signal', trim=trim)
        
        #Tidy into a seperate processed dataframe
        #Tidy into a seperate processed dataframe
        self.processed_dataframe = self.normalized_dataframe.ix[:,[
            self.energy_id,
            'signal'
        ]].copy()
        self.processed_dataframe.rename(columns={self.energy_id:'Energy / eV'}, inplace=True)
        self._AddLog('processed_dataframe created')

        
# Extra Functions for Compatibility
def load_file_to_dataframe(filename, beamline):
    """Extracts data to a dataframe"""
    if beamline == "APS: 4-ID-C":
        try:
            dataframe = _MDAdatafile(filename)
        except IndexError:
            return 'This data was not collected at {beamline}. Please update the edit this sample to switch to the appropriate beamline.'.format(beamline=beamline)
        else:
            df = dataframe.dataframe.rename(columns={IDC4.TEY_id: 'TEY', IDC4.TFY_id: 'TFY'})
            return df
    elif beamline == "ALS: 6.1.2":
        dataframe = _MDAdatafile(filename, flavour="ALS")
        df = dataframe.dataframe.rename(columns={ALS6312.TEY_id: 'TEY', ALS6312.TFY_id: 'TFY'})
        return df
    elif beamline == "ALS: 8.0.1":
        dataframe = _MDAdatafile(filename, flavour="ALS")
        df = dataframe.dataframe.rename(columns={ALS801.TEY_id: 'TEY', ALS801.TFY_id: 'TFY'})
        return df
