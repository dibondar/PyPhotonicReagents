
# Add main directory to enable imports  
if __name__ == '__main__' :
	import os
	os.sys.path.append(os.path.abspath('../..'))
	
	
from libs.gui.iterate_file import IterateFile 
import wx
import numpy as np
from operator import itemgetter 

#########################################################################

class SpectralScanViewer (IterateFile) :
	
	def UpdateFrame (self, event=None) :

		# for each individual (within a given iteration) load fitness and mean spectrum
		loaded_data = [ 
			( 	ind_grp["fitness"][...].sum(), 
				np.mean([ S[...] for S in ind_grp["spectra"].itervalues()], axis=0) )
					for ind_grp in self.GetCurrentFrame()["individuals"].itervalues() 
		]
		
		# Arrange by fitness
		loaded_data.sort(key=itemgetter(0), reverse=True)
		
		# Break into fitness and spectrum
		fitness, mean_spectra = zip(*loaded_data)
		
		# Displaying the spectrum
		self.fig.clear()
		self.fig.patch.set_facecolor('grey')
		axes = self.fig.add_subplot(111, axisbg='grey')
		
			
		for spectrum in mean_spectra :
			axes.plot (spectrum)
			
		axes.legend( ["%.2e" % f for f in fitness] )
		axes.set_title ("Spectrum")
		
		self.canvas.draw()	
		
#########################################################################

if __name__ == '__main__' :

	# Check weather the filename was given as a command line argument	
	import sys
	if len(sys.argv) == 2 : filename = sys.argv[1] 
	else : filename = None

	app = wx.App ()
	SpectralScanViewer ("optimization_iterations", filename=filename)
	app.MainLoop ()
