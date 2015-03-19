#######################################################################################
#
# This file contains implementation of the abstract class `BasicWindow` that
# implements basic functionality, such as saving and loading settings, stopping all jobs
#
#######################################################################################

from libs.dev.consts import * 
import wx, h5py, visvis
import numpy as np
import functools

from libs.gui.shaper_equalizer import PulseShaperEqualizer

def GetSmartAutoScaleRange (data, previous_bound=None) : 
	"""
	Return the tuple specifying plotting range for `data`.
	"""
	if not isinstance(previous_bound, tuple) or len(previous_bound) == 0:
		return ( min(data), max(data) )

	min_r = min(previous_bound)
	min_d = min(data)
	
	max_r = max(previous_bound) 
	max_d = max(data)
	
	if min_r > min_d or max_r < max_d : 
		# Bounds needs to be updated
	
		if min_d > 0 :
			min_r = 0.9*min_d
		else : 
			min_r = 1.1*min_d
			
		if max_d > 0 :
			max_r = 1.1*(max_d - min_r) + min_r
		else :
			max_r = 0.9*(max_d - min_r) + min_r
	
	return (min_r, max_r)

def SaveSettings (event=None, SettingsNotebook=None, filename = "settings.hdf5", title="Open HDF5 file to save settings", OpenDialog=True ) :
	"""
	Method for saving setting 
	"""
	if OpenDialog :
		# Ask user to select the file
		openFileDialog = wx.FileDialog(SettingsNotebook, title, "", filename, "HDF5 files (*.hdf5)|*.hdf5", 
						wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT | wx.FD_CHANGE_DIR)
		# Check whether user cancelled
		if openFileDialog.ShowModal() == wx.ID_CANCEL: return None	
		filename = openFileDialog.GetPath()
		
	with h5py.File (filename, 'a') as file_settings :

		# Crete the grope if it does not exist
		try : parameters_grp = file_settings["settings"] 
		except KeyError : parameters_grp = file_settings.create_group("settings")
			
		# Loop over all settings tab
		for SettingsTabName, SettingsTab in SettingsNotebook.settings_to_tabs.items() : 
			# Save all settings on a given tab
			try : del parameters_grp[SettingsTabName]
			except KeyError : pass
			grp = parameters_grp.create_group(SettingsTabName)
				
			for key, value in SettingsTab.GetSettings().items() : grp[key] = value
		
	# return valid filename
	return filename
	
class BasicWindow (wx.Frame) :
	"""
	Abstract class `BasicWindow` implements basic functionality, such as 
	saving and loading settings, stopping all jobs.
	
	The child of this abstract class must implement:
	
		`self.panel` (type: wx.Panel) panel with wx.Buttons for operating some tasks
		`self.SettingsNotebook` (type: wx.Notebook) containing setting tabs 
		`self.SettingsNotebook.settings_to_tabs`(type: dictionary) binding of setting names to setting tabs
	"""
	def on_close (self, event):
		"""
		Windows is about to be closed. Stop all jobs.
		"""
		self.StopAllJobs ()
		self.Destroy ()	
	
	def StopAllJobs (self) :
		"""
		Stop all tasks 
		"""
		for control in self.panel.GetChildren() :
			try :
				if isinstance(control, wx.Button) and control.GetLabel() != control.__start_label__ :
					control.SetLabel (control.__stop_label__)
					control.GetEventHandler().ProcessEvent(wx.PyCommandEvent(wx.EVT_BUTTON.typeId, control.GetId()))
					control.GetEventHandler().ProcessEvent(wx.PyCommandEvent(wx.EVT_LEFT_DOWN.typeId, control.GetId()))
			except AttributeError : pass
	
	#######################################################################################
	
	def CreateSaveSettingsButton (self) :
		"""
		Add a wx.Button to `self.panel` that will be binded to method `self.SaveSettings`
		"""
		self.save_settings_button = wx.Button (self.panel, label="Save settings...")
		self.SaveSettings = functools.partial(SaveSettings, SettingsNotebook=self.SettingsNotebook)
		self.Bind (wx.EVT_BUTTON, self.SaveSettings, self.save_settings_button)
		return self.save_settings_button
	
	#######################################################################################
	
	def CreateLoadSettingsButton (self) :
		"""
		Add a wx.Button to `self.panel` that will be binded to method `self.LoadSettings`
		"""
		self.load_settings_button = wx.Button (self.panel, label="Load settings...")
		self.Bind (wx.EVT_BUTTON, self.LoadSettings, self.load_settings_button)
		return self.load_settings_button
		
	def LoadSettings (self, event=None, title="Open HDF5 file to load settings") :
		"""
		Load settings. This method is closely related to <self.SaveSettings>
		"""
		openFileDialog = wx.FileDialog(self, title, "", "",
                                       "HDF5 files (*.hdf5)|*.hdf5", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR)
		# Check whether user cancelled
		if openFileDialog.ShowModal() == wx.ID_CANCEL: return None
		
		self.StopAllJobs()
		
		with h5py.File (openFileDialog.GetPath(), 'r') as file_settings :
			for SettingsTabName, SettingsTab in file_settings["settings"].items() :
				try :
					self.SettingsNotebook.settings_to_tabs[SettingsTabName].SetSettings(SettingsTab)
				except KeyError :
					print "Load Settings Error: Settings %s are ignored" %  SettingsTabName
					
		# return valid filename
		return openFileDialog.GetPath()
	
	#######################################################################################
	#
	#	Advanced features
	#
	#######################################################################################
	
	def CreateShowSpectrumButton (self) :
		"""
		Add a wx.Button to `self.panel` that will be binded to method `self.InteractivelyMeasureSpectrum`
		"""
		self.show_spectrum_button = wx.Button (self.panel)
		self.show_spectrum_button.__start_label__ = "Show spectrum"
		self.show_spectrum_button.__stop_label__ = "STOP measuring spectrum"
		self.show_spectrum_button.SetLabel (self.show_spectrum_button.__start_label__)
		self.Bind (wx.EVT_BUTTON, self.StartInteractivelyMeasureSpectrum, self.show_spectrum_button)
		return self.show_spectrum_button
		
	def StartInteractivelyMeasureSpectrum (self, event) :
		"""
		This method display spectrum
		"""
		button = self.show_spectrum_button
		
		if button.GetLabel() == button.__start_label__ :
			self.StopAllJobs()
			# get spectrometer's settings
			spect_settings = self.SettingsNotebook.Spectrometer.GetSettings()
			
			# Initiate spectrometer
			if self.Spectrometer.SetSettings(spect_settings) == RETURN_FAIL : return
			
			try : self.wavelengths = self.Spectrometer.GetWavelengths()
			except AttributeError : self.wavelengths = None
			
			# Clearing the figure
			visvis.cla(); visvis.clf();		
			
			# Set up timer to draw spectrum
			TIMER_ID = wx.NewId()
			self.spectrum_timer =  wx.Timer (self, TIMER_ID)
			self.spectrum_timer.Start (spect_settings["exposure_time"])
			wx.EVT_TIMER (self, TIMER_ID, self.DrawSpectrum)
			
			# Chose plotting options 
			try :
				self.is_autoscaled_spectrum = not(spect_settings["fix_vertical_axis"])
			except KeyError :
				self.is_autoscaled_spectrum = True
				
			if not self.is_autoscaled_spectrum :
				self.spectrum_plot_limits = ( spect_settings["vertical_axis_min_val"],
											spect_settings["vertical_axis_max_val"] )
			
			# Change button's label
			button.SetLabel (button.__stop_label__)
		
		elif button.GetLabel() == button.__stop_label__ :		
			# Stopping timer
			self.spectrum_timer.Stop()
			del self.spectrum_timer
			# Delete the parameter for auto-scaling
			del self.spectrum_plot_limits, self.is_autoscaled_spectrum
			
			# Delate visvis objects
			try : del self.__interact_2d_spectrum__
			except AttributeError : pass
			try : del self.__interact_1d_spectrum__
			except AttributeError : pass
			
			# Change button's label
			button.SetLabel (button.__start_label__) 
			
		else : raise ValueError("Label is not recognized") 
	
	def DrawSpectrum (self, event) :
		"""
		Draw spectrum interactively
		"""		
		spectrum = self.Spectrometer.AcquiredData() 
		if spectrum == RETURN_FAIL : return
		# Display the spectrum
		if len(spectrum.shape) > 1:
			try :
				self.__interact_2d_spectrum__.SetData(spectrum)
			except AttributeError :
				visvis.cla(); visvis.clf(); 	
				# Spectrum is a 2D image
				visvis.subplot(211)
				self.__interact_2d_spectrum__ = visvis.imshow(spectrum, cm=visvis.CM_JET)
				visvis.subplot(212)
				
			# Plot a vertical binning
			spectrum = spectrum.sum(axis=0)
			
		# Linear spectrum
		try :
			self.__interact_1d_spectrum__.SetYdata(spectrum)	
		except AttributeError :
			if self.wavelengths is None : 
				self.__interact_1d_spectrum__ = visvis.plot (spectrum, lw=3)
				visvis.xlabel ("pixels")
			else : 
				self.__interact_1d_spectrum__ = visvis.plot (self.wavelengths, spectrum, lw=3)
				visvis.xlabel("wavelength (nm)")
			visvis.ylabel("counts")
			
		if self.is_autoscaled_spectrum :
			# Smart auto-scale linear plot
			try :
				self.spectrum_plot_limits = GetSmartAutoScaleRange(spectrum, self.spectrum_plot_limits)
			except AttributeError :
				self.spectrum_plot_limits = GetSmartAutoScaleRange(spectrum)
		
		visvis.gca().SetLimits ( rangeY=self.spectrum_plot_limits )
		
		# Display the current temperature
		try : visvis.title ("Temperature %d (C)" % self.Spectrometer.GetTemperature() )
		except AttributeError : pass

	#######################################################################################
	
	def CreatePulseShaperEqualizerButton (self) :
		"""
		Add a wx.Button to `self.panel` that will be binded to method `self.AmplitudePhaseEqualizer`
		"""
		self.pulse_shaper_ampl_phase_equalizer_button = wx.Button (self.panel, label="Pulse shaper equalizer...")			
		self.Bind (wx.EVT_BUTTON, self.OpenPulseShaperEqualizer, self.pulse_shaper_ampl_phase_equalizer_button)
		return self.pulse_shaper_ampl_phase_equalizer_button
		
	def OpenPulseShaperEqualizer (self, event) :
		"""
		Start amplitude phase shaper equalizer
		"""
		self.PulseShaper.StopDevice()
		# Initiate pulse shaper
		if self.PulseShaper.Initialize( self.SettingsNotebook.PulseShaper.GetSettings() ) == RETURN_FAIL : return
		# Open a separate window
		PulseShaperEqualizer(self.PulseShaper).start()
		# Start spectrometer, if the corresponding button is added
		try : self.InteractivelyMeasureSpectrum(self.show_spectrum_button)
		except AttributeError : pass
	
	#######################################################################################
	
	def CreateRandomPhaseButton (self) :
		"""
		Add a wx.Button to `self.panel` that will be binded to method `self.SetRandomPhase`
		"""
		self.set_random_phase_button = wx.Button (self.panel, label="Set random phase")
		self.Bind (wx.EVT_BUTTON, self.SetRandomPhase, self.set_random_phase_button)
		return self.set_random_phase_button
		
	def SetRandomPhase (self, event) :
		"""
		Send random phase to the pulse shaper
		"""
		self.PulseShaper.StopDevice()
		# Initiate pulse shaper
		if self.PulseShaper.Initialize( self.SettingsNotebook.PulseShaper.GetSettings() ) == RETURN_FAIL : return
		# Get number of logical pixels
		param_num = self.PulseShaper.GetParamNumber()
		# Set random phase
		self.PulseShaper.SetAmplPhase( np.ones(param_num), np.random.rand(param_num) )

	#######################################################################################
	
	def CreateRandomAmplitudeButton (self) :
		"""
		Add a wx.Button to `self.panel` that will be binded to method `self.SetRandomAmplitude`
		"""
		self.set_random_ampl_button = wx.Button (self.panel, label="Set random amplitude")
		self.Bind (wx.EVT_BUTTON, self.SetRandomAmplitude, self.set_random_ampl_button)
		return self.set_random_ampl_button
		
	def SetRandomAmplitude (self, event) :
		"""
		Send random amplitude to the pulse shaper
		"""
		self.PulseShaper.StopDevice()
		# Initiate pulse shaper
		if self.PulseShaper.Initialize( self.SettingsNotebook.PulseShaper.GetSettings() ) == RETURN_FAIL : return
		# Get number of logical pixels
		param_num = self.PulseShaper.GetParamNumber()
		# Set random amplitude
		self.PulseShaper.SetAmplPhase( np.random.rand(param_num), np.zeros(param_num) )
		
	#######################################################################################
	
	def CreateZeroAmplitudeButton (self) :
		"""
		Add a wx.Button to `self.panel` that will be binded to method `self.SetZeroAmplitude`
		"""
		self.set_zero_ampl_button = wx.Button (self.panel, label="Set zero amplitudes")
		self.Bind (wx.EVT_BUTTON, self.SetZeroAmplitude, self.set_zero_ampl_button)
		return self.set_zero_ampl_button
		
	def SetZeroAmplitude (self, event) :
		"""
		Set zeros values of amplitude
		"""
		self.PulseShaper.StopDevice()
		# Initiate pulse shaper
		if self.PulseShaper.Initialize( self.SettingsNotebook.PulseShaper.GetSettings() ) == RETURN_FAIL : return
		# Get number of logical pixels
		param_num = self.PulseShaper.GetParamNumber()
		# Set random amplitude
		self.PulseShaper.SetAmplPhase( np.zeros(param_num), np.zeros(param_num) )
		
