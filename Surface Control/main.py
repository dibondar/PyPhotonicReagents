########################################################################
#
#	Main program for surface control
#
########################################################################

import wx
import h5py
import multiprocessing
from consts import * 
import numpy as np
# Hardware
from spectrometer_shamrock import ManagerShamrockSpectrometer, ShamrockSpectrometerTab
from pulse_shaper import ManagerShaper, PulseShaperTab
from thorlabs_apt_moving_stage import ManagerThorlabsAPTMovingStage
# Procedure
from rectangular_scan import ManagerRectangularScan, RectangularScanTab

# Real time plotting
import visvis
########################################################################

# starting <ManagerShamrockSpectrometer> in a separate process
#class ManagerShamrockSpectrometer(BaseManager): pass
#ManagerShamrockSpectrometer.register('ShamrockSpectrometer', ShamrockSpectrometer)

########################################################################

class SettingsNotebook (wx.Notebook) :
	"""
	GUI for listing all settings
	"""
	def __init__(self, parent):
		wx.Notebook.__init__(self, parent)
		 
		self.Spectrometer = ShamrockSpectrometerTab(self)
		self.AddPage (self.Spectrometer, "Spectra settings")
		 
		self.PulseShaper = PulseShaperTab(self)
		self.AddPage (self.PulseShaper, "Pulse shaper settings")

		self.RectangularScan = RectangularScanTab(self)
		self.AddPage (self.RectangularScan, "Rectangular scan")
		 
		# Dictionary to bind names to tabs for saving and loading settings
		self.settings_to_tabs = {"Spectrometer" : self.Spectrometer, "PulseShaper" : self.PulseShaper,
			"RectangularScan" : self.RectangularScan }
		 
########################################################################

class SurfaceControlExperiment (wx.Frame) :
	"""
	Application for running experiments
	"""
	def __init__ (self, parent) :
		# Starting spectrometer
		self.Spectrometer = ManagerShamrockSpectrometer()
		self.SpectrometerProc = self.Spectrometer.start()
		
		# Starting pulse shaper
		self.PulseShaper = ManagerShaper()
		self.PulseShaperProc = self.PulseShaper.start()
		
		# Starting moving stages (by specifying their serial numbers)
		self.MovingStageX = ManagerThorlabsAPTMovingStage(83843642)
		self.MovingStageXProc = self.MovingStageX.start()
				
		self.MovingStageY = ManagerThorlabsAPTMovingStage(83843641)
		self.MovingStageYProc = self.MovingStageY.start()

		# Create GUI
		dw, dh = wx.DisplaySize()
		wx.Frame.__init__ (self, parent, title="Surface control experiment", size=(0.9*dw, 0.88*dh) )
		
		self.ConstructGUI ()
		self.Center()
		self.Show ()
		wx.EVT_CLOSE (self, self.on_close)
	
	def on_close (self, event):
		"""
		Windows is about to be closed. Stop all timers.
		"""
		self.StopAllJobs ()
		self.Destroy ()	
	
	def ConstructGUI (self) :
		""" Build GUI """
		self.panel = wx.Panel(self)
		sizer = wx.GridBagSizer ()
		
		############################ Settings Notebook ############################
		self.SettingsNotebook = SettingsNotebook(self.panel)
		sizer.Add(self.SettingsNotebook, pos=(0, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT , border=10)

		############################ Command panel ############################
		boxsizer = wx.BoxSizer (wx.VERTICAL)
		
		# Test button
		test_button = wx.Button (self.panel, label="Test") 
		
		def OnTestButton (event) :
			self.PulseShaper.Initialize( self.SettingsNotebook.PulseShaper.GetSettings() )
			self.PulseShaper.Test()
			
		self.Bind (wx.EVT_BUTTON, OnTestButton, test_button)
		boxsizer.Add (test_button, flag=wx.EXPAND|wx.TOP, border=5)
		
		# Interactively display spectrum
		self.show_spectrum_button = wx.Button (self.panel)
		self.show_spectrum_button.__start_label__ = "Show spectrum"
		self.show_spectrum_button.__stop_label__ = "STOP measuring spectrum"
		self.show_spectrum_button.SetLabel (self.show_spectrum_button.__start_label__)
		self.Bind (wx.EVT_BUTTON, self.MeasureSingleSpectrum, self.show_spectrum_button)
		#self.show_spectrum_button.Bind(wx.EVT_BUTTON, self.MeasureSingleSpectrum)
		boxsizer.Add (self.show_spectrum_button, flag=wx.EXPAND, border=5)
		
		################## Rectangular scan button ##################
		self.rectangular_scan_button = wx.Button (self.panel)
		self.rectangular_scan_button.Bind (wx.EVT_LEFT_DOWN, self.PerformMeasurments)
		self.rectangular_scan_button.Bind (wx.EVT_LEFT_DCLICK, self.PerformMeasurments)
		boxsizer.Add(self.rectangular_scan_button, flag=wx.EXPAND, border=5)
		# Define labels
		self.rectangular_scan_button.__start_label__ 	= "Rectangular scan"
		self.rectangular_scan_button.__pause_label__ 	= "PAUSE scan"
		self.rectangular_scan_button.__resume_label__	= "RESUME scan"
		self.rectangular_scan_button.__stop_label__ 	= "STOP scan"
		self.rectangular_scan_button.SetLabel (self.rectangular_scan_button.__start_label__)
		# Specify the measurements settings
		self.rectangular_scan_button.__measurmenet_manager__ 		= ManagerRectangularScan
		self.rectangular_scan_button.__measurmenet_manager_args__	= (self.Spectrometer, self.MovingStageX, self.MovingStageY)
		self.rectangular_scan_button.__tab_settings__				= "RectangularScan"
		#self.rectangular_scan_button.__post_process__ = self.RectangularScanPostProcess 
		
		# Save settings
		self.save_settings_button = wx.Button (self.panel, label="Save settings...")
		self.Bind (wx.EVT_BUTTON, self.SaveSettings, self.save_settings_button)
		boxsizer.Add(self.save_settings_button, flag=wx.EXPAND|wx.TOP, border=5)
		
		# Load settings
		self.load_settings_button = wx.Button (self.panel, label="Load settings...")
		self.Bind (wx.EVT_BUTTON, self.LoadSettings, self.load_settings_button)
		boxsizer.Add(self.load_settings_button, flag=wx.EXPAND|wx.TOP, border=5)
		
		sizer.Add(boxsizer, pos=(1, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT , border=10)
		########################### End of constructing panel ######################################
		self.panel.SetSizer (sizer)
		
		############################# Setting visvis #######################################
		Figure = app.GetFigureClass()
		self.fig = Figure(self)
		
		boxsizer = wx.BoxSizer (wx.HORIZONTAL)
		boxsizer.Add(self.panel, 1, wx.EXPAND)
		boxsizer.Add(self.fig._widget, 2, wx.EXPAND)

		self.SetSizer (boxsizer)
		self.SetAutoLayout(True)
		self.Layout() 	
		
	def __del__ (self) :	
		# Close moving stages
		self.MovingStageX.exit(); self.MovingStageXProc.join()
		self.MovingStageY.exit(); self.MovingStageYProc.join()
		
		# Close spectrometer
		self.Spectrometer.exit(); self.SpectrometerProc.join() 
		
		# Close pulse shaper
		self.PulseShaper.exit(); self.PulseShaperProc.join()
		
	def PerformMeasurments (self, event) :
		"""
		A universal wrapper for performing measurements
		"""
		# Extracting which button was clicked
		try :
			button = event.GetEventObject()
			# Mouse double clicking stops scanning
			if event.GetEventType() == wx.wxEVT_LEFT_DCLICK  : button.SetLabel (button.__stop_label__)
		except AttributeError : button = event

		if button.GetLabel() == button.__start_label__ :
			self.StopAllJobs ()
			
			# get spectrometer's settings
			settings = self.SettingsNotebook.Spectrometer.GetSettings()
			# Initiate spectrometer
			if self.Spectrometer.SetSettings(settings) == RETURN_FAIL : return
			
			# Extract the measurements settings 
			tab = self.SettingsNotebook.settings_to_tabs[button.__tab_settings__]
			settings = tab.GetSettings()
			if "filename" not in settings : settings["filename"] = "results.hdf5"
			
			# Save the global settings
			result = self.SaveSettings( default_filename=settings["filename"], title=button.__start_label__)
			if not isinstance(result, basestring) : 
				# User did not chose the file name
				return
			
			if settings["filename"] != result :
				# Update the file name if user chosen different
				settings["filename"] = result
				tab.SetSettings (settings)
			
			# Saving the filename for postprocesing 
			button.__results_filename__ = result
			
			# Start scanning via the corresponding manager
			button.__running_manager__ = button.__measurmenet_manager__(settings, *button.__measurmenet_manager_args__)
			
			# Start timer to monitor weather measurement is over
			TIMER_ID = wx.NewId()
			button.__scanning_timer__ = wx.Timer (self, TIMER_ID)
			button.__scanning_timer__.Start (2000) # check every 2 seconds
			
			def check_weather_scanning_finished (event) : 
				if not button.__running_manager__.is_running () : 
					button.SetLabel (button.__stop_label__); self.PerformMeasurments (button)
			
			wx.EVT_TIMER (self, TIMER_ID, check_weather_scanning_finished)
			
			# Changing the button's label 
			button.SetLabel (button.__pause_label__)
			
		elif button.GetLabel() == button.__pause_label__ :
			button.__running_manager__.pause(); button.SetLabel (button.__resume_label__)
			
		elif button.GetLabel() == button.__resume_label__ :
			button.__running_manager__.resume(); button.SetLabel (button.__pause_label__)
			
		elif button.GetLabel() == button.__stop_label__ :
			# Stop timer
			button.__scanning_timer__.Stop()
			del button.__scanning_timer__
			# Stop measurements
			button.__running_manager__.stop(); 
			del button.__running_manager__
			button.SetLabel (button.__start_label__)
			
			try : # Start post processing, if present
				button.__post_process__(button.__results_filename__)
			except AttributeError : pass
			
		else : raise ValueError ("Unrecognised button label")
		
	def StopAllJobs (self) :
		"""
		Stop  tasks 
		"""
		for control in self.panel.GetChildren() :
			try :
				if isinstance(control, wx.Button) and control.GetLabel() != control.__start_label__ :
					control.SetLabel (control.__stop_label__)
					control.GetEventHandler().ProcessEvent(wx.PyCommandEvent(wx.EVT_BUTTON.typeId, control.GetId()))
			except AttributeError : pass
		
	def MeasureSingleSpectrum (self, event=None) :
		"""
		Button <self.show_spectrum_button> was clicked
		"""
		button = self.show_spectrum_button
		
		if button.GetLabel() == button.__start_label__ :
			self.StopAllJobs()
			# get spectrometer's settings
			spect_settings = self.SettingsNotebook.Spectrometer.GetSettings()
			
			# Initiate spectrometer
			if self.Spectrometer.SetSettings(spect_settings) == RETURN_FAIL : return
			self.wavelengths = self.Spectrometer.GetWavelengths()
			
			# Clearing the figure
			visvis.clf()
		
			def draw_spectrum (event) :
				"""Timer function """
				spectrum = self.Spectrometer.AcquiredData() 
				if spectrum == RETURN_FAIL : return
				# Display the spectrum
				
				############### Take the log of spectrum ##########
				#spectrum = spectrum / float(spectrum.max())
				#np.log10(spectrum, out=spectrum)
				##############################
				
				ax = visvis.gca()
				ax.Clear()	
				visvis.plot (self.wavelengths, spectrum)
				visvis.xlabel("wavelength (nm)")
				visvis.ylabel("counts")
				
				# Display the current temperature
				visvis.title ("Temperature %d (C)" % self.Spectrometer.GetTemperature() )
				
			# Set up timer to draw spectrum
			TIMER_ID = wx.NewId()
			self.spectrum_timer =  wx.Timer (self, TIMER_ID)
			self.spectrum_timer.Start (spect_settings["exposure_time"])
			
			# Change button's label
			button.SetLabel (button.__stop_label__)
			wx.EVT_TIMER (self, TIMER_ID, draw_spectrum)
			
		elif button.GetLabel() == button.__stop_label__ :
			# Stopping timer
			self.spectrum_timer.Stop()
			del self.spectrum_timer
			# Change button's label
			button.SetLabel (button.__start_label__) 
			
		else : raise ValueError("Label is not recognized") 
		
	def SaveSettings (self, event=None, default_filename = "settings.hdf5", title="Open HDF5 file to save settings" ) :
		"""
		Button <self.save_settings_button> was clicked
		"""
		openFileDialog = wx.FileDialog(self, title, "", default_filename, "HDF5 files (*.hdf5)|*.hdf5", 
							wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT | wx.FD_CHANGE_DIR)
		# Check whether user cancelled
		if openFileDialog.ShowModal() == wx.ID_CANCEL: return None	
		
		with h5py.File (openFileDialog.GetPath(), 'w') as file_settings :
			# create general settings 
			parameters_grp = file_settings.create_group("settings")
			# Loop over all settings tab
			for SettingsTabName, SettingsTab in self.SettingsNotebook.settings_to_tabs.items() :
				# Save all settings on a given tab
				grp = parameters_grp.create_group(SettingsTabName)
				for key, value in SettingsTab.GetSettings().items() : grp[key] = value
		
		# return valid filename
		return openFileDialog.GetPath()
		
	def LoadSettings (self, event) :
		"""
		Button <self.load_settings_button> was clicked. This method is closely related to <self.SaveSettings>
		"""
		openFileDialog = wx.FileDialog(self, "Open HDF5 file to load settings", "", "",
                                       "HDF5 files (*.hdf5)|*.hdf5", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR)
		# Check whether user canceled
		if openFileDialog.ShowModal() == wx.ID_CANCEL: return	
		
		self.StopAllJobs()
		
		with h5py.File (openFileDialog.GetPath(), 'r') as file_settings :
			for SettingsTabName, SettingsTab in file_settings["settings"].items() :
				self.SettingsNotebook.settings_to_tabs[SettingsTabName].SetSettings(SettingsTab)
		
		
#########################################################################
if __name__ == '__main__' :
	multiprocessing.freeze_support()
	app = visvis.use('wx')
	app.Create()
	SurfaceControlExperiment (None)
	app.Run()
