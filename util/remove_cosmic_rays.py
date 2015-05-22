"""
Utilities for removing cosmic rays 
"""
import numpy as np

def replace_brightest (img, npixels_to_remove=200 ) :
	"""
	Replace npixels_to_remove brightest pixels from image by 
	less brightest values
	"""
	img = np.ravel(img) 
	N = img.size - npixels_to_remove
	indx = np.argpartition(img, N) 
	img[ indx[N:] ] = img[ indx[N-1] ] 

def replace_bins (img, bins=200) :
	"""
	Perform binning of the image and remove bins that have no or few pixels 
	"""
	hist, bin_edges = np.histogram(img, bins)
	indx = np.argmin(hist).min()
	np.clip(img, 0, bin_edges[max(0,indx-1)], out=img)