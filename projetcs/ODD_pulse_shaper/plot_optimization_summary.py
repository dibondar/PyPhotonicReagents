
import numpy as np
import wx
import h5py
import matplotlib.pyplot as plt
from itertools import product
import random
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
		optimization_summary[ key ] = value[...]#[:35]
		
	# Load reference signal info
	reference_signal = {}
	for key, value in F["reference_signal_summary"].iteritems() :
		reference_signal[ key ] = value[...]#[:35]
		
############################## Plot data ##############################

plt.subplot(221)
for value in optimization_summary.itervalues() :
	plt.plot( value )

plt.xlabel("iteration number")
plt.ylabel("obj function")
plt.title ("optimization summary")

plt.legend( optimization_summary.keys(), loc=2 )

##############################
plt.subplot(222)
for value in reference_signal.itervalues() :
	plt.plot( value )

plt.legend( ["chanel %s" % x for x in reference_signal.keys()] )

plt.xlabel("iteration number")
plt.ylabel("reference ")
plt.title ("Reference signal fluctuation")

##############################
plt.subplot(223)
obj = optimization_summary["avg"].astype(np.float)

for ref in reference_signal.itervalues() :
	R = obj / ref[:obj.size].astype(np.float)
	R /= R.max()
	plt.plot( R )

plt.legend( [ "avg / ref %s" % x for x in reference_signal.keys() ], loc=4 )

plt.xlabel("iteration number")
plt.title("Obj function normalized to reference signal")

"""
plt.subplot(223)
for ref, obj in product( reference_signal.itervalues(), optimization_summary.itervalues() ) :
	R = ref.astype(np.float)[:obj.size]
	R /= R.max()
	plt.plot( obj/R )

plt.legend( [
	"%s / (chanel %s ref)" % (obj, ref) for ref, obj in 
		product( reference_signal.iterkeys(), optimization_summary.iterkeys() )
], loc=2 )

plt.xlabel("iteration number")
plt.title("Obj function normalized to reference signal")
"""

##############################
plt.subplot(224)

# Iterator over all possible fluorescence matrices form from  
#fluoresence_matrix_ref_signal = product( *[ 
#	combinations(ref, len(reference_signal)) for ref in reference_signal.values() 
#] )

 
def det_fluoresence() :
	while True :
		# Randomly draw a fluorescence matrix out of available reference data
		F = [ random.sample(ref, len(reference_signal)) for ref in reference_signal.values() ]
		yield np.linalg.det(F)
	
det_data = np.fromiter( det_fluoresence(), np.float, 20000 )
S = det_data.size

# Get rid of duplicate data
det_data = np.unique(det_data)
if det_data.size != S : 
	print "There were duplicates in the Monte Carlo simulation"

plt.hist( det_data, bins=100 )

plt.xlabel("value of det")
plt.ylabel("counts")

# If the data passes normality test, then fit the distribution with normal variable
if normaltest(det_data,)[1] < 0.2 : 
	plt.title( "Contribution of noise to the determinant,\n calculations are based on reference signal data")
else :
	# The data seams to be normally distributed
	plt.title( "Contribution of noise to the determinant. Normal fit: mu = %.2e sigma = %.2e" % norm.fit(det_data) )

##############################

plt.show()