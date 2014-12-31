
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
	
	def IniLoad(self, data_file) :
	
		# Find out what kind of pulse shapes are saved in the file
		pulse_shaping_option = str( data_file["settings/ODD_GA/pulse_shaping_option"][...] )
		self.plot_ampl_phase =  ( pulse_shaping_option == "amplitude and phase")
		
		if not self.plot_ampl_phase :
				self.plot_title = {	"amplitude only" : "Amplitude shapes", 
									"phase only" : "Phase shapes" }[pulse_shaping_option]
	
	def UpdateFrame (self, event=None) :

		# for each individual (within a given iteration) load fitness and pulse shape 
		loaded_data = [ 
			( 	ind_grp["fitness"][...].sum(), 
				ind_grp["pulse_shape"][...] 	)
					for ind_grp in self.GetCurrentFrame()["individuals"].itervalues() 
		]
		
		# Arrange by fitness
		loaded_data.sort(key=itemgetter(0), reverse=True)
		loaded_data = loaded_data[:5]
		# Break into fitness and pulse_shape
		fitness, pulse_shapes = zip(*loaded_data)
		
		# Displaying the pulse shapes
		self.fig.clear()
		self.fig.patch.set_facecolor('grey')
		
		if self.plot_ampl_phase :
			# Display amplitude and phase separately
			
			ampl_axes = self.fig.add_subplot(211, axisbg='grey')
			ampl_axes.set_title ("Amplitude")
			
			phase_axes = self.fig.add_subplot(212, axisbg='grey')
			phase_axes.set_title ("Phase")
			
			for pulse in pulse_shapes :
				ampl_axes.plot ( pulse[:pulse.size/2] )
				phase_axes.plot ( pulse[pulse.size/2:] )
			
			ampl_axes.legend( ["%.2e" % f for f in fitness] )
			
		else :
			# Phase/Ampl shaping only
			axes = self.fig.add_subplot(111, axisbg='grey')
			
			for pulse in pulse_shapes :
				axes.plot (pulse)
				
			axes.legend( ["%.2e" % f for f in fitness] )
			
			axes.set_title (self.plot_title)
			
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
