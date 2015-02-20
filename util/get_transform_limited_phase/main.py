########################################################################
#
#	Wrapper for ODD discrimination algorithm build on DEAP library (http://deap.readthedocs.org/en/latest/api/algo.html)
#
########################################################################
# Add main directory to enable imports  
if __name__ == '__main__' :
	import os
	os.sys.path.append(os.path.abspath('../..'))
########################################################################

import wx

# Real time plotting
import visvis

# GUI components
from libs.gui.basic_window import BasicWindow

from ga_tab import GA_Tab

# Hardware
from libs.dev.spectrometer_ocean_optics import ManagerOceanOpticsSpectrometer as ManagerSpectrometer
from libs.dev.spectrometer_ocean_optics import OceanOpticsSpectrometerTab as SpectrometerTab
#from libs.dev.camera_istar import ManagerIStarCamera as ManagerSpectrometer
#from libs.dev.camera_istar import IStarCameraTab as SpectrometerTab

from libs.dev.pulse_shaper import ManagerShaper, PulseShaperTab

########################################################################

class SettingsNotebook (wx.Notebook) :
	"""
	GUI for listing all settings
	"""
	def __init__(self, parent, DevSpectrometer, DevPulseShaper ):
		"""
		`DevSpectrometer` is a spectrometer manager
		"""
		wx.Notebook.__init__(self, parent)
		
		self.ODD_GA = GA_Tab(self)
		self.AddPage(self.ODD_GA, "ODD GA")
		
		self.Spectrometer = SpectrometerTab(self, DevSpectrometer)
		self.AddPage (self.Spectrometer, "Spectrometer")
		 
		self.PulseShaper = PulseShaperTab(self, DevPulseShaper)
		self.AddPage (self.PulseShaper, "Pulse shaper")

		# Dictionary to bind names to tabs for saving and loading settings
		self.settings_to_tabs = {"Spectrometer" : self.Spectrometer, 
			"PulseShaper" : self.PulseShaper, "ODD_GA" : self.ODD_GA }		

########################################################################

class ODDExperiment  (BasicWindow) :

	def __init__ (self, parent) :
		# Starting spectrometer
		self.Spectrometer = ManagerSpectrometer()
		self.SpectrometerProc = self.Spectrometer.start()
		
		# Starting pulse shaper
		self.PulseShaper = ManagerShaper()
		self.PulseShaperProc = self.PulseShaper.start()
		
		# Create GUI
		dw, dh = wx.DisplaySize()
		wx.Frame.__init__ (self, parent, title="ODD for multiple fluoresce marker concentration measurements",
								size=(0.9*dw, 0.88*dh) )
		
		self.ConstructGUI ()
		self.Center()
		self.Maximize()
		self.Show ()
		wx.EVT_CLOSE (self, self.on_close)
		
	def __del__ (self) :	
		# Close spectrometer
		self.Spectrometer.exit(); self.SpectrometerProc.join() 
		
		# Close pulse shaper
		self.PulseShaper.exit(); self.PulseShaperProc.join()
	
	def ConstructGUI (self) :
		""" Build GUI """
		self.panel = wx.Panel(self)
		sizer = wx.GridBagSizer ()
		
		############################ Settings Notebook ############################
		self.SettingsNotebook = SettingsNotebook(self.panel, self.Spectrometer, self.PulseShaper)
		sizer.Add(self.SettingsNotebook, pos=(0, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT , border=10)

		############################ Command panel ############################
		boxsizer = wx.BoxSizer (wx.VERTICAL)
		
		# Interactively display spectrum
		boxsizer.Add (self.CreateShowSpectrumButton(), flag=wx.EXPAND, border=5)
		
		# Separator
		boxsizer.Add (wx.StaticText(self.panel), flag=wx.EXPAND, border=5)
	
		# Separator
		boxsizer.Add (wx.StaticText(self.panel), flag=wx.EXPAND, border=5)
		
		# Send random phase to the pulse shaper
		boxsizer.Add (self.CreateRandomPhaseButton(), flag=wx.EXPAND, border=5)
		# Send random amplitude to the pulse shaper
		boxsizer.Add (self.CreateRandomAmplitudeButton(), flag=wx.EXPAND, border=5)
		# Send zero amplitude and zero phase to the pulse shaper
		boxsizer.Add (self.CreateZeroAmplitudeButton(), flag=wx.EXPAND, border=5)
		
		# Open pulse shaper equalizer
		boxsizer.Add (self.CreatePulseShaperEqualizerButton(), flag=wx.EXPAND, border=5)
		
		# Separator
		boxsizer.Add (wx.StaticText(self.panel), flag=wx.EXPAND, border=5)
		
		# Save settings
		boxsizer.Add( self.CreateSaveSettingsButton(), flag=wx.EXPAND, border=5)
		# Load settings
		boxsizer.Add( self.CreateLoadSettingsButton(), flag=wx.EXPAND, border=5)
		
		sizer.Add(boxsizer, pos=(1, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.GROW, border=10)
		########################### End of constructing panel ######################################
		self.panel.SetSizer (sizer)
		
		############################# Setting visvis #######################################
		Figure = app.GetFigureClass()
		self.fig = Figure(self)
		
		boxsizer = wx.BoxSizer (wx.HORIZONTAL)
		boxsizer.Add(self.panel, 0.5, wx.EXPAND)
		boxsizer.Add(self.fig._widget, 2, wx.EXPAND)
		
		#########################################################################################			
		self.SetSizer (boxsizer)
		self.SetAutoLayout(True)
		self.Layout() 
		
#########################################################################

if __name__ == '__main__' :
	app = visvis.use('wx')
	app.Create()
	ODDExperiment (None)
	app.Run()