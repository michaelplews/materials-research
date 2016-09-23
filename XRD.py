# -*- coding: utf-8 -*-
#XRD technique
import csv, os, math, xml.etree.ElementTree as ET, numpy as np
import json
from pprint import pprint

import dm3_lib as tem

from xml.etree import ElementTree
import zipfile

import numpy as np, pandas as pd, matplotlib.pyplot as plt

class BM11CSVfile(object):
	'''Loads CSV files exported from 11-BM .mda Mail-In Service'''
	filename = ""
	shortname = ""
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

			
	def norm_dataframe(self, lower_lim, upper_lim):
		two_theta, intensity = self.dataframe.index, self.dataframe[' intensity']
		if lower_lim and upper_lim:
			max_x, max_y = self.max_in_range(two_theta, intensity, lower_lim, upper_lim)
		else:
			max_y = np.amax(intensity)
		min_y = np.amin(intensity)
		a = max_y-min_y
		intensity[:] = [((x-min_y)/a) * 100 for x in intensity]
		return two_theta, intensity

	def plot(self, normalized=True, color="red", legend="", lower_lim="", upper_lim="", style=""):
		if style:
			style.use=(style)
		if normalized is True:
			two_theta, intensity = self.norm_dataframe(lower_lim, upper_lim)
		else:
			intensity = self.dataframe[' intensity']
			two_theta = self.dataframe.index
		if legend:
			plt.plot(two_theta, intensity, linewidth = 1, label=legend, color=color)
		elif self.shortname:
			plt.plot(two_theta, intensity, linewidth = 1, label=self.shortname, color=color)
		else:
			plt.plot(two_theta, intensity, linewidth = 1, label='no_label', color=color)
		plt.legend(loc = 1).draggable(True)				# fontsize and other parameters removed as we now pickle our pgf styles
		plt.subplots_adjust(hspace=0, wspace=0)

#Taken (with permission) from https://github.com/m3wolf/scimap and edited. Thanks Mark!
class BrukerBrmlFile(object):
	'''Loads data from a Bruker .brml v4 file to an object'''
	filename = ""
	shortname = ""
	
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
	
	@property	
	def norm_dataframe(self):
		two_theta, intensity = self.dataframe
		max = np.amax(intensity)
		min = np.amin(intensity)
		a = max-min
		intensity[:] = [((x-min)/a) * 100 for x in intensity]
		return two_theta, intensity

	def new_source(self, old_wavelength, new_wavelength):
		two_theta, intensity = self.norm_dataframe
		two_theta[:] = [2*math.degrees(np.arcsin((new_wavelength*np.sin(math.radians(x/2)))/old_wavelength)) for x in two_theta]
		return two_theta, intensity

	def plot(self, normalized=True, color="red", legend="", new_source=False, old_wavelength="", new_wavelength="", style=""):
		if style:
			plt.style.use(style)
		if new_source:
			two_theta, intensity = self.new_source(old_wavelength, new_wavelength)		
		elif normalized:
			two_theta, intensity = self.norm_dataframe
		else:
			two_theta, intensity = self.dataframe
		if legend:
			plt.plot(two_theta, intensity, linewidth = 1, label=legend, color=color)
		elif self.shortname:
			plt.plot(two_theta, intensity, linewidth = 1, label=self.shortname, color=color)
		else:
			plt.plot(two_theta, intensity, linewidth = 1, label='no_label', color=color)
		plt.axis([10, 80, 0, 110])
		plt.legend(loc = 1).draggable(True)				# fontsize and other parameters removed as we now pickle our pgf styles
		plt.subplots_adjust(hspace=0, wspace=0)

class XYFile(object):
	'''Class that imports data from ASCII .xy file to an object'''
	filename = ""
	shortname = ""
	
	def __init__(self, filename, shortname=""):
		self.filename = filename
		self.shortname = shortname
		
	@property
	def dataframe(self): #Used by other functions
		dataframe = np.genfromtxt(self.filename, delimiter=None, names=['two_theta', 'intensity'],skip_header=0, autostrip=True)
		return dataframe['two_theta'], dataframe['intensity']
		
	@property
	def norm_dataframe(self):
		two_theta, intensity = self.dataframe
		max = np.amax(intensity)
		min = np.amin(intensity)
		a = max-min
		intensity[:] = [((x-min)/a) * 100 for x in intensity]
		return two_theta, intensity
		
	def plot(self, normalized=True, color="red", legend=""):
		if normalized:
			two_theta, intensity = self.norm_dataframe
		else:
			two_theta, intensity = self.dataframe
		if legend:
			plt.plot(two_theta, intensity, linewidth = 1, label=legend, color=color)
		elif self.shortname:
			plt.plot(two_theta, intensity, linewidth = 1, label=self.shortname, color=color)
		else:
			plt.plot(two_theta, intensity, linewidth = 1, label='no_label', color=color)
		plt.axis([10, 80, 0, 110])
		plt.legend(loc = 1).draggable(True)				# fontsize and other parameters removed as we now pickle our pgf styles
		plt.subplots_adjust(hspace=0, wspace=0)
		
class ICDDXmlFile(object):
	"""Class that imports ICDD xml file containing XRD reference data. 'flavour' ("thousand" or "hundred") optional argument gives the option to not divide all values by 10"""
	filename = ""
	pdf_number = ""
	legend = str("")
	shortname = ""
	xtal_system = ""
	flavour = ""
	
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
	
	@property
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
		plt.legend(loc = 1).draggable(True)				# fontsize and other parameters removed as we now pickle our pgf styles
		if hkl and hkl_values:
			for index in range (0, len(hkl_values)):
				if float(x[index]) < lim:
					plt.annotate(hkl_values[index], xy=(x[index], y[index]), xycoords='data', fontsize=10, horizontalalignment='left', verticalalignment='bottom')
		plt.xlim(float(10),float(lim))
		plt.ylim(float(0),float(110))

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
		plt.legend(loc = 1).draggable(True)				# fontsize and other parameters removed as we now pickle our pgf styles
		if hkl and hkl_values:
			for index in range (0, len(hkl_values)):
				if float(x[index]) < lim:
					plt.annotate(hkl_values[index], xy=(x[index], y[index]), xycoords='data', fontsize=10, horizontalalignment='left', verticalalignment='bottom')
		#plt.xlim(float(10),float(lim))
		#plt.ylim(float(0),float(110))

class MaterProjJSON(object):
	"""Class that imports downloaded XRD JSON files from materialsproject.org for comparison against other XRD patterns"""
	filename = ""
	mp_number = ""
	legend = ""
	shortname = ""
	xtal_system = ""
	flavour = ""
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

	def plot(self, color="red", legend="", hkl=False, lim=80):
		if legend:
			plt.vlines(self.two_theta,-10, self.amplitude, label=legend,colors=color) 
		elif self.legend:
			plt.vlines(self.two_theta,-10, self.amplitude, label=self.legend,colors=color) 
		elif self.shortname:
			plt.vlines(self.two_theta,-10, self.amplitude, label=self.shortname,colors=color) 
		else:
			plt.vlines(self.two_theta,-10, self.amplitude,label='no_label',colors=color) 
		plt.legend(loc = 1).draggable(True)				# fontsize and other parameters removed as we now pickle our pgf styles
		if hkl and self.hkl:
			for index in range (0, len(self.hkl)):
				if float(self.two_theta[index]) < lim and float(self.amplitude[index]) > 10:
					plt.annotate(self.hkl[index], xy=(self.two_theta[index], self.amplitude[index]), xycoords='data', fontsize=10, horizontalalignment='left', verticalalignment='bottom')
		plt.xlim(float(10),float(lim))
		plt.ylim(float(0),float(110))
			  
def retro_plot_xrd(filename, color="red", legend="", number_of_stacked=0):
	#normalize_xrd(filename)
	data = np.genfromtxt(filename + "_normalized.txt", delimiter='\t', names=['x', 'y'],skip_header=1, autostrip=True)
	data['y'] += 10 * number_of_stacked
	if legend:
		plt.plot(data['x'], data['y'], linewidth = 1, label=legend, color=color)
	else:
		plt.plot(data['x'], data['y'], linewidth = 1, label='no_label', color=color)
	plt.axis([10, 80, -10, 110 + number_of_stacked*10])
	plt.legend(loc = 1).draggable(True)				# fontsize and other parameters removed as we now pickle our pgf styles
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
	plt.legend(loc = 1).draggable(True)				# fontsize and other parameters removed as we now pickle our pgf styles
	plt.subplots_adjust(hspace=0)
