########################################################################
#
# ODD : Measuring concentration of fluorescent proteins using 
#	motorized filter wheel to shape pulses
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

from odd_tab_filter_wheel import ODD_Tab_FilterWheel

# Hardware
#from libs.dev.spectrometer_ocean_optics import ManagerOceanOpticsSpectrometer as ManagerSpectrometer
#from libs.dev.spectrometer_ocean_optics import OceanOpticsSpectrometerTab as SpectrometerTab
from libs.dev.camera_istar import ManagerIStarCamera as ManagerSpectrometer
from libs.dev.camera_istar import IStarCameraTab as SpectrometerTab

from libs.dev.motorized_filter_wheel import ManagerFilterWheel, FilterWheelTab
from libs.dev.sample_switcher import ManagerSampleSwitcher, SampleSwitcherTab

########################################################################

class SettingsNotebook (wx.Notebook) :
	"""
	GUI for listing all settings
	"""
	def __init__(self, parent, DevSpectrometer, DevSampleSwitcher, FilterWheel ):
		"""
		`DevSpectrometer` is a spectrometer manager
		"""
		wx.Notebook.__init__(self, parent)
		
		self.ODD_FilterWheel = ODD_Tab_FilterWheel(self)
		self.AddPage(self.ODD_FilterWheel, "ODD")
		
		self.Spectrometer = SpectrometerTab(self, DevSpectrometer)
		self.AddPage (self.Spectrometer, "Spectrometer")
		 
		self.SampleSwitcher = SampleSwitcherTab(self, DevSampleSwitcher) 
		self.AddPage (self.SampleSwitcher, "Sample switcher")
	
		self.FilterWheel = FilterWheelTab(self, FilterWheel)
		self.AddPage (self.FilterWheel, "Filter wheel")
	
		# Dictionary to bind names to tabs for saving and loading settings
		self.settings_to_tabs = {
			"Spectrometer" 		: self.Spectrometer,
			"ODD_FilterWheel" 	: self.ODD_FilterWheel,
			"FilterWheel" 		: self.FilterWheel,
			"SampleSwitcher" 	: self.SampleSwitcher  }		

########################################################################

class ODDExperiment_FilterWheel  (BasicWindow) :

	def __init__ (self, parent) :
		# Starting spectrometer
		self.Spectrometer = ManagerSpectrometer()
		self.SpectrometerProc = self.Spectrometer.start()
		
		# Start motorized filter wheel
		self.FilterWheel = ManagerFilterWheel()
		self.FilterWheelProc = self.FilterWheel.start()
		
		# Start sample switcher 
		self.SampleSwitcher = ManagerSampleSwitcher()
		self.ManagerSampleSwitcherProc = self.SampleSwitcher.start()
		
		# Create GUI
		dw, dh = wx.DisplaySize()
		wx.Frame.__init__ (self, parent, title="ODD for multiple fluoresce marker concentration measurements using motorized filter wheel",
								size=(0.9*dw, 0.88*dh) )
		
		self.ConstructGUI ()
		self.Center()
		self.Maximize()
		self.Show ()
		wx.EVT_CLOSE (self, self.on_close)
		
	def __del__ (self) :	
		# Close spectrometer
		self.Spectrometer.exit(); self.SpectrometerProc.join() 
	
		# Close motorized filter wheel
		self.FilterWheel.exit(); self.FilterWheelProc.join()
	
		# Close sample switcher
		self.SampleSwitcher.exit(); self.ManagerSampleSwitcherProc.join()
	
	def ConstructGUI (self) :
		""" Build GUI """
		self.panel = wx.Panel(self)
		sizer = wx.GridBagSizer ()
		
		############################ Settings Notebook ############################
		self.SettingsNotebook = SettingsNotebook(self.panel, self.Spectrometer, self.SampleSwitcher, self.FilterWheel)
		sizer.Add(self.SettingsNotebook, pos=(0, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT , border=10)

		############################ Command panel ############################
		boxsizer = wx.BoxSizer (wx.VERTICAL)
		
		# Interactively display spectrum
		boxsizer.Add (self.CreateShowSpectrumButton(), flag=wx.EXPAND, border=5)
		
		# Separator
		boxsizer.Add (wx.StaticText(self.panel), flag=wx.EXPAND, border=5)
		
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
	ODDExperiment_FilterWheel (None)
	app.Run()