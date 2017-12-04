# -*- coding: utf-8 -*-
#EChem technique

import numpy as np, matplotlib.pyplot as plt, pandas as pd

class processedMPTFile(object):
    """Loads data from a processed (IQxnE) .mpt file produced in EC-Lab or BT-Lab. Will also work for .txt files exported by the software."""
    filename = ""
    shortname = ""
    channel = ""
    date = ""
    device = ""
    material = ""
    initial_state = ""
    electrolyte = ""
    comments = []
    mass_am = ""
    mol_weight = ""
    surface_area = ""
    skip_line = 1
    dataframe = ""
    start_of_method = 3
    efficiency_column = 0
    capacity_column = 0
    cycle_number_column = 0
       
    def _lines_to_skip(self, filename):
        file = open(filename, 'r', encoding='latin-1')
        for line in file:
            if "Nb header lines : " in line:
                return int(line.strip("Nb header lines : "))
                break
        return 1
    
    def _get_columns(self):
        file = open(self.filename, 'r', encoding='latin-1')
        line = file.readlines()[self.skip_line-1]
        columns = line.split('\t')
        return columns
    
    def __init__(self, filename, shortname=""):
        columns = []
        self.comments = []
        self.filename = filename
        self.shortname = shortname
        self.skip_line = self._lines_to_skip(filename)
        file = open(filename, 'r', encoding='latin-1')
        columns = self._get_columns()
        for column in range(0, len(columns)):
            if 'cycle number' in columns[column]:             
                self.cycle_number_column = column
            if 'Q discharge' in columns[column]:
                self.capacity_column = column
            if 'Efficiency' in columns[column]:
                self.efficiency_column = column
        i = 0
        for line in file:
            i += 1
            if "Run on channel : " in line:
                self.channel = line.strip("Run on channel : ").strip('\n')
            elif "Acquisition started on : " in line:
                self.date = line.strip("Acquisition started on : ").strip('\n').replace('\t', ' ')
            elif "Device : " in line:
                self.device = line.strip("Device : ").strip('\n').replace('\t', ' ')
            elif "Electrode material : " in line:
                self.material = line.strip("Electrode material : ").strip('\n').replace('\t', ' ')
            elif "Initial state : " in line:
                self.initial_state = line.strip("Initial state : ").strip('\n').replace('\t', ' ')
            elif "Electrolyte : " in line:
                self.electrolyte = line.strip("Electrolyte : ").strip('\n').replace('\t', ' ')
            elif "Comments : " in line:
                self.comments.append(line.strip("Comments : ").replace('\t', ' '))
            elif "Mass of active material : " in line:
                self.mass_am = line.strip("Mass of active material : ").strip('\n').replace('\t', ' ')
            elif "Electrode surface area : " in line:
                self.surface_area = line.strip("Electrode surface area : ").strip('\n').replace('\t', ' ')
                self.start_of_method = i + 2
            elif "Molecular weight of active material (at x = 0) : " in line:
                self.mol_weight = line.strip("Molecular weight of active material (at x = 0) : ").strip('\n').replace('\t', ' ')
        self.dataframe = np.genfromtxt(filename, delimiter=None, skip_header=self.skip_line, autostrip=True, unpack=True)
                
    @property
    def all(self):
        print ('Filename: ' + self.filename)
        print ("Shortname: " + self.shortname)
        print ("Channel: " + self.channel)
        print ("Date Started: " + self.date)
        print ("Device: " + self.device)
        print ("Active Material: " + self.material)
        print ("Initial State: " + self.initial_state)
        print ("Electrolyte: " + self.electrolyte)
        print ("Comments: ")
        for line in self.comments:
            print ("\t" + line.strip('\n'))
        print ("Mass of active material: " + self.mass_am)
        print ('Molecular weight at x=0: ' + str(self.mol_weight))
        print ('Electrode Surfance Area: ' + self.surface_area)
    
    @property
    def method(self):
        with open(self.filename, encoding='latin-1') as file:
            for lines in file.readlines()[self.start_of_method:self.skip_line-1]:
                print (lines.strip('\n'))
    
    def _get_columns(self):
        file = open(self.filename, 'r', encoding='latin-1')
        line = file.readlines()[self.skip_line-1]
        columns = line.split('\t')
        return columns

    @property
    def show_columns(self):
        columns = self._get_columns()
        for column in range (0,len(columns)):
            print (str(column).zfill(2) +' | ' +columns[column])

#class to import ASCII .txt and .mpt files produced in EC-Lab and BT-Lab
class MPTFile(object):
    """Loads data from a .mpt file produced in EC-Lab or BT-Lab. Will also work for .txt files exported by the software."""
    filename = ""
    shortname = ""
    channel = ""
    date = ""
    device = ""
    material = ""
    initial_state = ""
    electrolyte = ""
    comments = []
    mass_am = ""
    mol_weight = ""
    surface_area = ""
    skip_line = 0
    dataframe = ""
    start_of_method = 3
    column_list = []            # Assigned in _lines_to_skip
       
    def _lines_to_skip(self, filename):
        file = open(filename, 'r', encoding='latin-1')
        for line in file:
            if "Nb header lines : " in line:
                self.column_list = line.strip("Nb header lines : ").split('\t')
                return int(line.strip("Nb header lines : "))
    
    def __init__(self, filename, shortname=""):
        self.comments = []
        self.filename = filename
        self.shortname = shortname
        self.skip_line = self._lines_to_skip(filename)
        file = open(filename, 'r', encoding='latin-1')
        i = 0
        for line in file:
            i += 1
            if "Run on channel : " in line:
                self.channel = line.strip("Run on channel : ").strip('\n')
            elif "Acquisition started on : " in line:
                self.date = line.strip("Acquisition started on : ").strip('\n').replace('\t', ' ')
            elif "Device : " in line:
                self.device = line.strip("Device : ").strip('\n').replace('\t', ' ')
            elif "Electrode material : " in line:
                self.material = line.strip("Electrode material : ").strip('\n').replace('\t', ' ')
            elif "Initial state : " in line:
                self.initial_state = line.strip("Initial state : ").strip('\n').replace('\t', ' ')
            elif "Electrolyte : " in line:
                self.electrolyte = line.strip("Electrolyte : ").strip('\n').replace('\t', ' ')
            elif "Comments : " in line:
                self.comments.append(line.strip("Comments : ").replace('\t', ' '))
            elif "Mass of active material : " in line:
                self.mass_am = line.strip("Mass of active material : ").strip('\n').replace('\t', ' ')
            elif "Electrode surface area : " in line:
                self.surface_area = line.strip("Electrode surface area : ").strip('\n').replace('\t', ' ')
                self.start_of_method = i + 2
            elif "Molecular weight of active material (at x = 0) : " in line:
                self.mol_weight = line.strip("Molecular weight of active material (at x = 0) : ").strip('\n').replace('\t', ' ')
        # self.dataframe = np.genfromtxt(filename, delimiter=None, skip_header=self.skip_line, autostrip=True, unpack=True)
        # line = file.readlines()[self.skip_line-1]
        # self.column_list = line.split('\t')
        self.dataframe = pd.read_csv(filename, sep='\t', skiprows=self.skip_line-1, header=0, index_col=0, encoding='latin-1')
        file.close()
        
    @property
    def all(self):
        print ('Filename: ' + self.filename)
        print ("Shortname: " + self.shortname)
        print ("Channel: " + self.channel)
        print ("Date Started: " + self.date)
        print ("Device: " + self.device)
        print ("Active Material: " + self.material)
        print ("Initial State: " + self.initial_state)
        print ("Electrolyte: " + self.electrolyte)
        print ("Comments: ")
        for line in self.comments:
            print ("\t" + line.strip('\n'))
        print ("Mass of active material: " + self.mass_am)
        print ('Molecular weight at x=0: ' + str(self.mol_weight))
        print ('Electrode Surfance Area: ' + self.surface_area)
    
    @property
    def method(self):
        with open(self.filename, encoding='latin-1') as file:
            for lines in file.readlines()[self.start_of_method:self.skip_line-1]:
                print (lines.strip('\n'))
    
    @property
    def show_columns(self):
        file = open(self.filename, 'r', encoding='latin-1')
        line = file.readlines()[self.skip_line-1]
        columns = line.split('\t')
        for column in range (0,len(columns)):
            print (str(column).zfill(2) +' | ' +columns[column])
    
    def _moving_average(self, a, n=3) :
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n
    
    def specific_capacity(self, capacity_column, mass_am):
        self.dataframe = self.dataframe.assign(specific_capacity=self.dataframe[capacity_column]/mass_am)
                        
    def plot(self, x, y, legend="", color="red", cycle='', cycle_column="", disconnect=False, end_point=False, end_point_color='', linewidth=1.5, *args, **kwargs):
        if cycle or cycle is 0:
            xdata, ydata=[],[]
            for data_point in range (0, len(self.dataframe[cycle_column])):
                if self.dataframe[cycle_column].values[data_point] == cycle:
                    xdata.append(self.dataframe[x].values[data_point])
                    ydata.append(self.dataframe[y].values[data_point])
                    #print (xdata, ydata)
                else:
                    pass 
            if disconnect:
                pos = np.where(np.abs(np.diff(xdata)) >= 0.001 )[0]+1
                xdata = np.insert(xdata, pos, np.nan)
                ydata = np.insert(ydata, pos, np.nan)
            else:
                pass            
            if legend is "None":
                plt.plot(xdata, ydata, color=color, linewidth=linewidth)
            elif legend:
                plt.plot(xdata, ydata, color=color, linewidth=linewidth, label=legend)
            elif self.shortname:
                plt.plot(xdata,ydata, color=color, linewidth=linewidth, label=self.shortname)
            else:
                plt.plot(xdata,ydata, color=color, linewidth=linewidth, label="no_label")
        
        else:
            xdata, ydata = self.dataframe[x], self.dataframe[y]
            if disconnect:
                pos = np.where(np.abs(np.diff(xdata)) >= 0.01 )[0]+1
                xdata = np.insert(xdata, pos, np.nan)
                ydata = np.insert(xdata, pos, np.nan)
            else:
                pass
            if legend is "None":
                plt.plot(xdata, ydata, color=color, linewidth=linewidth)
            elif legend:
                plt.plot(xdata, ydata, color=color, linewidth=linewidth, label=legend)
            elif self.shortname:
                plt.plot(xdata, ydata, color=color, linewidth=linewidth, label=self.shortname)
            else:
                plt.plot(xdata, ydata, color=color, linewidth=linewidth, label="no_label")
        if end_point:
            if end_point_color:
                plt.scatter(self.dataframe[x][-1], self.dataframe[y][-1], color=end_point_color, *args, **kwargs)
            else:
                plt.scatter(self.dataframe[x][-1], self.dataframe[y][-1], color=color, *args, **kwargs)
        plt.legend(loc = 1, frameon = False)
        plt.xlim(np.amin(self.dataframe[x])-0.1,np.amax(self.dataframe[x])+0.1)
        plt.ylim(np.amin(self.dataframe[y])-0.2,np.amax(self.dataframe[y])+0.2)


    def plot_diffcap(self, x=10, y=22, n=5, legend="", color="red"):
        """Plots differential capacity. Moving average of n=5 by default"""
        Ecell = self.dataframe[x]
        q = self.dataframe[y]
        qo = self.dataframe[x][0]
        q = [x-qo for x in q]
        dE = [x-y for x,y in zip(Ecell, Ecell[2:])]
        dq = [x-y for x,y in zip(q,q[2:])]
        dE = np.array(dE)
        dq = np.array(dq)
        dqdv = dq/dE
        dqdv[dqdv == 0] = None
        Ecell = Ecell[np.logical_not(np.isnan(dqdv))] #removes None data from corresponding data in the other axis
        dqdv = dqdv[np.logical_not(np.isnan(dqdv))] #removes None data from original data to make arrays the same size 
        ma = self._moving_average(abs(dqdv), n=n)
        if legend:
            plt.plot(Ecell[n-1:],ma, color=color, linewidth=1.5, label=legend)
        elif self.shortname:
            plt.plot(Ecell[n-1:],ma, color=color, linewidth=1.5, label=self.shortname)
        else:
            plt.plot(Ecell[n-1:],ma, color=color, linewidth=1.5, label="no_label")
