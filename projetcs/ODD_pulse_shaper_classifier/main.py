########################################################################
#
#	Wrapper for ODD discrimination algorithm 
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
from libs.gui.info_tab import InfoTab

from odd_tab import ODD_Tab
from combinatorial_chirp_scan import CombinatorialChirpScan_Tab

# Hardware
#from libs.dev.spectrometer_ocean_optics import ManagerOceanOpticsSpectrometer as ManagerSpectrometer
#from libs.dev.spectrometer_ocean_optics import OceanOpticsSpectrometerTab as SpectrometerTab
from libs.dev.camera_istar import ManagerIStarCamera as ManagerSpectrometer
from libs.dev.camera_istar import IStarCameraTab as SpectrometerTab

from libs.dev.pulse_shaper import ManagerShaper, PulseShaperTab
from libs.dev.sample_switcher import ManagerSampleSwitcher, SampleSwitcherTab
from libs.dev.motorized_filter_wheel import ManagerFilterWheel, FilterWheelTab

########################################################################

class SettingsNotebook (wx.Notebook) :
	"""
	GUI for listing all settings
	"""
	def __init__(self, parent, DevSpectrometer, DevSampleSwitcher, DevPulseShaper, DevFilterWheel ):
		"""
		`DevSpectrometer` is a spectrometer manager
		"""
		wx.Notebook.__init__(self, parent)
		
		self.InfoTab = InfoTab(self) 
		self.AddPage(self.InfoTab, "Info")
		
		self.ODD_Tab = ODD_Tab(self)
		self.AddPage(self.ODD_Tab, "ODD")
		
		self.CombinatorialChirpScan_Tab = CombinatorialChirpScan_Tab(self)
		self.AddPage(self.CombinatorialChirpScan_Tab, "Combinatorial scan")
		
		self.FilterWheel =  FilterWheelTab(self, DevFilterWheel)
		self.AddPage (self.FilterWheel, "Filter wheel")
		
		self.Spectrometer = SpectrometerTab(self, DevSpectrometer)
		self.AddPage (self.Spectrometer, "Spectrometer")
		 
		self.SampleSwitcher = SampleSwitcherTab(self, DevSampleSwitcher) 
		self.AddPage (self.SampleSwitcher, "Sample switcher")
		 
		self.PulseShaper = PulseShaperTab(self, DevPulseShaper)
		self.AddPage (self.PulseShaper, "Pulse shaper")

		# Dictionary to bind names to tabs for saving and loading settings
		self.settings_to_tabs = { "Spectrometer" : self.Spectrometer, 
			"PulseShaper" : self.PulseShaper, "ODD_Tab" : self.ODD_Tab, 
			"CombinatorialScan" : self.CombinatorialChirpScan_Tab,
			"SampleSwitcher" : self.SampleSwitcher, "Info" : self.InfoTab,
			"FilterWheel" : self.FilterWheel }		

########################################################################

class ODDExperiment  (BasicWindow) :

	def __init__ (self, parent) :
		# Starting spectrometer
		self.Spectrometer = ManagerSpectrometer()
		self.SpectrometerProc = self.Spectrometer.start()
		
		# Starting pulse shaper
		self.PulseShaper = ManagerShaper()
		self.PulseShaperProc = self.PulseShaper.start()
		
		# Start sample switcher 
		self.SampleSwitcher = ManagerSampleSwitcher()
		self.ManagerSampleSwitcherProc = self.SampleSwitcher.start()
		
		# Start filter wheel
		self.FilterWheel = ManagerFilterWheel()
		self.FilterWheelProc = self.FilterWheel.start()
		
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
	
		# Close sample switcher
		self.SampleSwitcher.exit(); self.ManagerSampleSwitcherProc.join()
		
		# Close filter wheel
		self.FilterWheel.exit(); self.FilterWheelProc.join()
	
	def ConstructGUI (self) :
		""" Build GUI """
		self.panel = wx.Panel(self)
		sizer = wx.GridBagSizer ()
		
		############################ Settings Notebook ############################
		self.SettingsNotebook = SettingsNotebook(self.panel, self.Spectrometer, self.SampleSwitcher, self.PulseShaper, self.FilterWheel)
		sizer.Add(self.SettingsNotebook, pos=(0, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT , border=10)

		############################ Command panel ############################
		boxsizer = wx.BoxSizer (wx.VERTICAL)
		
		# Interactively display spectrum
		boxsizer.Add (self.CreateShowSpectrumButton(), flag=wx.EXPAND, border=5)
		
		# Separator
		boxsizer.Add (wx.StaticText(self.panel), flag=wx.EXPAND, border=5)
	
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