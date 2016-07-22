# -*- coding: utf-8 -*-
#Classes and functions of XAS experiments

import datetime, numpy as np, operator, pandas as pd, matplotlib.pyplot as plt
from scipy import interpolate

class igorTXTFile(object):
	"""Loads .txt wave files exported by IGOR Pro"""
	filename = ""
	shortname = ""
	dataframe = ""	

	def __init__(self, filename, shortname=""):
		self.filename = filename
		self.shortname = shortname
		self.dataframe = np.genfromtxt(filename, delimiter='\t', skip_header=1, autostrip=True, unpack=True)

	def plot(self, color="red", legend=""):
		if legend:
			plt.plot(self.dataframe[0],self.dataframe[1], color=color, linewidth=1.5, label=legend)
		elif self.shortname:
			plt.plot(self.dataframe[0],self.dataframe[1], color=color, linewidth=1.5, label=self.shortname)
		else:
			plt.plot(self.dataframe[0],self.dataframe[1], color=color, linewidth=1.5, label='no_label')
				
		plt.legend(loc = 1, frameon = False)
		plt.xlim(np.amin(self.dataframe[0]),np.amax(self.dataframe[0]))
		plt.ylim(np.amin(self.dataframe[1]),np.amax(self.dataframe[1])*1.1)
	
	def max_in_range(self, low, high, plot=True, do_return=False):
		"""Finds the maximum value of y in a given range of x"""
		y_values = self.dataframe[1][np.logical_and(low < self.dataframe[0], self.dataframe[0] < high)]
		x_values = self.dataframe[0][np.logical_and(low < self.dataframe[0], self.dataframe[0] < high)]
		index_max_y = y_values.argmax()
		max_y = y_values[index_max_y]
		max_x = x_values[index_max_y]
		if plot:
			plt.plot(max_x,max_y,'om', color='red')
		if do_return:
			return max_x, max_y

class _MDAdatafile(object):
	"""Loads .0001 data files produced by 4-ID-C at Argonne National Lab"""
	filename = ""
	basename = ""
	ext = ""
	header_lines = 0
	dataframe = ""
	column_index = []
	scan_datetime = ""
	
	def __init__(self, filename, header_lines=0, flavour="IDC4"):
		self.filename = filename
		self.ext = filename[-4:]
		self.basename = filename[:-5].split(r'/')[-1]
		file = open(filename, 'r')
		self.column_index = []
		self.header_lines = header_lines
		i = 0
		if flavour == "ICD4":
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
				if "Time of Day" in line:
					self.header_lines = i
			self.dataframe = pd.read_csv(filename, sep='\t', skiprows=self.header_lines-1, header=0, index_col=0)
			self.dataframe['Rounded Energy / eV'] = pd.Series(np.around(self.dataframe['Energy'], decimals=1), index=self.dataframe.index)
			self.dataframe.index = self.dataframe['Rounded Energy / eV']

	@property
	def available_signals(self):
		print ('Column\t|\tDetector')
		for line in range (0, len(self.column_index)):
			print (str(line) + '\t|\t' + self.column_index[line])
	
class IDC4(object):
	"""Loads .0001 data files produced by 4-ID-C at Argonne National Lab, Uses the _MDAFile class to average data scans and produce a normalized array."""
	directory = ""
	basename = ""
	start = ""
	end = ""
	exclude = []
	shortname= ""
	normalized_dataframe = ""
	processed_dataframe = ""
	beamline = "4-ID-C, Advanced Photon Source"
	_log = []

	#Variables that will change if instrument set up is adjusted
	Energy_id = '[1-D Positioner 1]  4idc1:SGM1:Energy\t SGM1:Energy\t LINEAR\t eV\t 4idc1:SGM1:EnergyRBV\t Energy readback\t eV'	
	TFY_id = '[1-D Detector  10]  4idc1:scaler1_calc5.VAL\t \t '
	TEY_id = '[1-D Detector   9]  4idc1:scaler1_calc4.VAL\t \t '
	REF_id = ''
	STD_id = '[1-D Detector  11]  4idc1:scaler1_calc6.VAL\t \t '

	def __init__(self, directory, basename, start, end, exclude=None, shortname="", TFY_smooth=7):
		self._log = []	#reset the log to be empty
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
		self.normalized_dataframe['TFY'] = self._ScaleAbs(self.TFY_id, 'TFY')
		self.normalized_dataframe['sTFY'] = self._ScaleAbs(self.TFY_id, 'TFY', smooth=TFY_smooth)
		self._AddLog('TFY smoothed: rolling mean window = 7')
		self.normalized_dataframe['TEY'] = self._ScaleRef(self.TEY_id, 'TEY')	
		self.normalized_dataframe['STD'] = self._ScaleRef(self.STD_id, 'STD')

		#Tidy into a seperate processed dataframe
		self.processed_dataframe = self.normalized_dataframe.ix[:,[
			self.Energy_id,
			'TFY',
			'sTFY',
			'TEY',
			'STD'
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
	
	@property
	def log(self):
		"""Prints formatted log output"""
		for line in self._log:
			print (line)

	def plot(self, signal, color="red", legend=''):
		signal = self.processed_dataframe[signal].values
		energy = self.processed_dataframe['Energy / eV'].values
		if legend:
			plt.plot(energy, signal, linewidth = 1, label=legend, color=color)
		elif self.shortname: 
			plt.plot(energy, signal, linewidth = 1, label=self.shortname, color=color)
		else:
			plt.plot(energy, signal, linewidth = 1, label='no_label', color=color)
		plt.axis([np.amin(energy), np.amax(energy), 0, np.amax(signal)*1.1])
		plt.legend(loc = 2, frameon = False).draggable(True)
		#plt.subplots_adjust(hspace=0, wspace=0)	

	def align(self, x_from, x_to):
		"""Applies the appropriate multiplication/division to shift waves in order to account for beam drift and other sources of expermental error. Should be used on a similar peak where x_from=original x location and x_to=desired x location."""		
		energy = self.processed_dataframe['Energy / eV']
		if x_to > x_from:
			energy *= x_to/x_from
		else:
			energy /= x_from/x_to
		self._AddLog(self.shortname + " aligned from " + str(x_from) + "eV to " + str(x_to) + "eV")

	def max_in_range(self, signal, low, high, plot=True, do_return=False):
		"""Finds the maximum value of y in a given range of x"""
		get_signal = self.processed_dataframe[signal].values
		energy = self.processed_dataframe['Energy / eV'].values
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
		"""Writes to object.log property to allow the user to observe and changes that have been made to the data since the object was initialised."""
		self._log.append('{time}:\t{message}'.format(time=datetime.datetime.now(), message=message))

	def _SumData(self):
		"""Internal function. Averages data across multiple _MDAFile objects"""
		normalized_dataframe = []
		normalized_dataframe = pd.concat([scan.dataframe for scan in self._MDAlist]).groupby(level=0).mean()
		return normalized_dataframe

	def _ScaleAbs(self, column_id, name, head=5, tail=5, smooth=None):
		"""Internal function. Scales data absolutely between 1 and 0 depending on the wave shape"""
		endloc = self.normalized_dataframe['Rounded Energy / eV'].iloc[-1]
		v_minloc = self.normalized_dataframe[column_id].idxmin()
		new_frame = self.normalized_dataframe[column_id].copy()
		if v_minloc < endloc/2+v_minloc/2:
			new_frame -= new_frame.head(head).mean()
			new_frame /= new_frame.tail(tail).mean()
		else:
			new_frame -= new_frame.tail(tail).mean()
			new_frame /= new_frame.head(head).mean()
		if smooth:
			new_frame = pd.rolling_mean(new_frame, smooth, smooth, center=True) 
		new_frame.columns = ['Rounded Energy / eV', name]
		return new_frame

	def _ScaleRef(self, column_id, name):
		"""Internal function. Scaled data between 1 and 0 irrespective of wave shape"""
		new_frame = self.normalized_dataframe[column_id].copy()
		new_frame -= new_frame.min()
		new_frame /= new_frame.max()
		new_frame.columns = ['Rounded Energy / eV', name]
		return new_frame

class DATFile(object):
	"""So far only one case of this extension appearing. When JF sends reference data as a request."""
	filename = ""
	shortname = ""
	dataframe = ""	

	def __init__(self, filename, shortname=""):
		self.filename = filename
		self.shortname = shortname
		self.dataframe = np.genfromtxt(filename, delimiter='\t', skip_header=1, autostrip=True, unpack=True)

	def plot(self, x=0, y=1, color="red", legend=""):
		if legend:
			plt.plot(self.dataframe[x],self.dataframe[y], color=color, linewidth=1, label=legend)
		elif self.shortname:
			plt.plot(self.dataframe[x],self.dataframe[y], color=color, linewidth=1, label=self.shortname)
		else:
			plt.plot(self.dataframe[x],self.dataframe[y], color=color, linewidth=1, label='no_label')
				
		plt.legend(loc = 1, frameon = False)
		plt.xlim(np.amin(self.dataframe[x]),np.amax(self.dataframe[x]))
		plt.ylim(np.amin(self.dataframe[y]),np.amax(self.dataframe[y])*1.1)
	
	def max_in_range(self, low, high, x=0, y=1, plot=True, do_return=False):
		"""Finds the maximum value of y in a given range of x"""
		y_values = self.dataframe[y][np.logical_and(low < self.dataframe[x], self.dataframe[x] < high)]
		x_values = self.dataframe[x][np.logical_and(low < self.dataframe[x], self.dataframe[x] < high)]
		index_max_y = y_values.argmax()
		max_y = y_values[index_max_y]
		max_x = x_values[index_max_y]
		if plot:
			plt.plot(max_x,max_y,'om', color='red')
		if do_return:
			return max_x, max_y

class ALS6312(object):
	"""Loads .txt data files produced by beamline 6312 at The Advanced Light Source, Lawrence Berkeley National Laboratory. The _MDAFile class is not needed as data was pre-averaged at the beamline before release. Before processing all files were renamed using bash scripts:

	$ for i in *-Avg.txt ; do mv "$i" "${i/-Avg.txt/.txt}" ; done
	$ for i in SigScan* ; do mv "$i" "${i%".txt"}" ; done
	$ for i in SigScan* ; do mv "$i" "${i/#"SigScan"/"SigScan."}" ; done

This renamed the files from "SigScan22222-Avg.txt" to "SigScan.22222". This makes it easier to use the IDC4 object code.

"""
	directory = ""
	basename = ""
	start = ""
	end = ""
	exclude = []
	shortname= ""
	normalized_dataframe = ""
	processed_dataframe = ""
	beamline = "6321, Advanced Light Source"
	_log = []

	#Variables that will change if instrument set up is adjusted
	Energy_id = 'Energy'	
	TFY_id = 'Channeltron'
	TEY_up_id = 'TEY_up'
	TEY_dn_id = 'TEY_dn'
	I0_id = 'Izero'
	REF_id = None
	STD_id = None

	def __init__(self, directory, basename, start, end, exclude=None, shortname="", TFY_smooth=7):
		self._log = []	#reset the log to be empty
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
		self._MDAlist = [_MDAdatafile(directory + basename + '.' + str(ext), flavour='ALS') 
						for ext in range(self.start, self.end+1)
						if ext not in self.exclude]
		self.normalized_dataframe = self._SumData()
		self._AddLog('Object created (__init__)')
		self._AddLog('normalized_dataframe created')		

		#More dataframes can be appended if needed
		self.normalized_dataframe['TFY'] = self._ScaleAbs(self.TFY_id, 'TFY', divisor=self.I0_id)
		self.normalized_dataframe['sTFY'] = self._ScaleAbs(self.TFY_id, 'TFY', divisor=self.I0_id, smooth=TFY_smooth)
		self._AddLog('TFY smoothed: rolling mean window = 7')
		self.normalized_dataframe['TEY'] = self._ScaleRef(self.TEY_up_id, 'TEY', divisor=self.I0_id)	
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
	
	@property
	def log(self):
		"""Prints formatted log output"""
		for line in self._log:

			print (line)

	def plot(self, signal, color="red", legend=''):
		signal = self.processed_dataframe[signal].values
		energy = self.processed_dataframe['Energy / eV'].values
		if legend:
			plt.plot(energy, signal, linewidth = 1, label=legend, color=color)
		elif self.shortname: 
			plt.plot(energy, signal, linewidth = 1, label=self.shortname, color=color)
		else:
			plt.plot(energy, signal, linewidth = 1, label='no_label', color=color)
		plt.axis([np.amin(energy), np.amax(energy), 0, np.amax(signal)*1.1])
		plt.legend(loc = 2, frameon = False).draggable(True)
		#plt.subplots_adjust(hspace=0, wspace=0)	

	def align(self, x_from, x_to):
		"""Applies the appropriate multiplication/division to shift waves in order to account for beam drift and other sources of expermental error. Should be used on a similar peak where x_from=original x location and x_to=desired x location."""		
		energy = self.processed_dataframe['Energy / eV']
		if x_to > x_from:
			energy *= x_to/x_from
		else:
			energy /= x_from/x_to
		self._AddLog(self.shortname + " aligned from " + str(x_from) + "eV to " + str(x_to) + "eV")

	def max_in_range(self, signal, low, high, plot=True, do_return=False):
		"""Finds the maximum value of y in a given range of x"""
		get_signal = self.processed_dataframe[signal].values
		energy = self.processed_dataframe['Energy / eV'].values
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
		"""Writes to object.log property to allow the user to observe and changes that have been made to the data since the object was initialised."""
		self._log.append('{time}:\t{message}'.format(time=datetime.datetime.now(), message=message))

	def _SumData(self):
		"""Internal function. Averages data across multiple _MDAFile objects"""
		normalized_dataframe = []
#		normalized_dataframe = pd.concat([scan.dataframe for scan in self._MDAlist]).groupby(level=0).mean()
		normalized_dataframe = [scan.dataframe for scan in self._MDAlist]
		return normalized_dataframe

	def _ScaleAbs(self, column_id, name, divisor=None, head=5, tail=5, smooth=None):
		"""Internal function. Scales data absolutely between 1 and 0 depending on the wave shape"""
		new_frame = self.normalized_dataframe[column_id].copy()
		if divisor:
			new_frame = new_frame/self.normalized_dataframe[divisor].copy()
		endloc = self.normalized_dataframe['Rounded Energy / eV'].iloc[-1]

		v_minloc = new_frame.idxmin()
		if v_minloc < endloc/2+v_minloc/2:
			new_frame -= new_frame.head(head).mean()
			new_frame /= new_frame.tail(tail).mean()
		else:
			new_frame -= new_frame.tail(tail).mean()
			new_frame /= new_frame.head(head).mean()
		if smooth:
			new_frame = pd.rolling_mean(new_frame, smooth, smooth, center=True) 
		new_frame.columns = ['Rounded Energy / eV', name]
		return new_frame

	def _ScaleRef(self, column_id, name, divisor=None):
		"""Internal function. Scaled data between 1 and 0 irrespective of wave shape"""
		new_frame = self.normalized_dataframe[column_id].copy()
		if divisor:
			new_frame = new_frame/self.normalized_dataframe[divisor].copy()
		new_frame -= new_frame.min()
		new_frame /= new_frame.max()
		new_frame.columns = ['Rounded Energy / eV', name]
		return new_frame
