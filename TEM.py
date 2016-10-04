# -*- coding: utf-8 -*-
#TEM technique
import csv, os, math, xml.etree.ElementTree as ET, numpy as np
from xml.etree import ElementTree
import zipfile
import numpy as np, pandas as pd, matplotlib.pyplot as plt, dm3_lib as dm3
from matplotlib.colors import Normalize
import matplotlib.patches as patches

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
#		scale bar location dict
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
