import numpy as np
import wx
import h5py
import matplotlib.pyplot as plt
import scipy.stats as stats
from itertools import combinations_with_replacement
from scipy.stats import norm, normaltest

try : filename
except NameError :
	app = wx.App()
	
	openFileDialog = wx.FileDialog (None, "Chose HDF5 file for viewing", "", "", \
				"HDF5 files (*.hdf5)|*.hdf5", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR)

	# Check whether user cancelled
	if openFileDialog.ShowModal() == wx.ID_CANCEL : 
		raise ValueError ("User did not chose the file to view")
	else : filename = openFileDialog.GetPath()

	del app, openFileDialog
	
############################## Load data from file ##############################
with h5py.File (filename, 'r') as F :

	# Load summary data
	optimization_summary = {}
	for key, value in F["optimization_summary"].iteritems() :
		optimization_summary[ key ] = value[...]
		
	# Load reference signal info
	reference_signal = {}
	for key, value in F["reference_signal_summary"].iteritems() :
		reference_signal[ key ] = value[...]
		
############################## Plot QQ plots for every reference signal ##############################

# Find out how to partition the plot into suplots  
points = [ int(np.sqrt(len(reference_signal))) ]
partitions = []

while not len(partitions) :
	partitions = [
		_ for _ in combinations_with_replacement(points, 2 ) 
		if np.prod(_) >= len(reference_signal)
	]
	points.extend( [ min(points) - 1, max(points) + 1] )
	
n1_plots, n2_plots = min(partitions, key=np.prod)

# Plot all QQ plots
n = 1
for key, ref in reference_signal.iteritems() :
	plt.subplot(n1_plots, n2_plots, n) 
	n += 1
	stats.probplot(ref, dist="norm", plot=plt)
	plt.title("QQ plot for %s" % key)
	
plt.show()

############################## Plot correlations ##############################

def Cov (x,y) :
	x = x - x.mean()
	y = y - y.mean()
	cov = np.correlate(x, y, mode='same')
	return cov[len(cov)/2:]

for ref1, ref2 in combinations_with_replacement( reference_signal.itervalues(), 2 ) :
	plt.plot( Cov(ref1, ref2) )
	
plt.legend( [ 
	"%s %s" % ref for ref in combinations_with_replacement( reference_signal.keys(), 2 ) 
] )

plt.show()

