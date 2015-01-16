########################################################################
#
# 	optimal dynamic discrimination experiment using motorized filter wheel
#	for pulse shaping
# 
########################################################################

import wx
import h5py
import itertools
import numpy as np
from scipy.optimize import nnls
from scipy.ndimage.filters import gaussian_filter

import visvis 

from libs.gui.hardware_control import HardwareGUIControl
from libs.gui.basic_window import SaveSettings 

from libs.dev.consts import * 

########################################################################
		
class _SpectrumPostProcess :
	"""
	Container for prepossessing spectrum
	"""
	################ Methods of post processing spectra ################
	@classmethod
	def VertBinned (cls, spectrum) :
		if len(spectrum.shape) == 2 : return spectrum.sum(axis=0)
		else : return spectrum
	
	@classmethod
	def Smooth (cls, spectrum) :
		return gaussian_filter(spectrum, sigma=3)
	
	@classmethod
	def Sum (cls, spectrum) :
		return np.array( [spectrum.sum()] )
	
	@classmethod
	def VertBinnedSmooth (cls, spectrum) :
		return cls.Smooth( cls.VertBinned(spectrum) )
	
	@classmethod
	def SmoothSum (cls, spectrum) :
		return cls.Sum( cls.VertBinnedSmooth(spectrum) )
	
	################ Service implementations ################
	@classmethod
	def _GetDict (cls) :
		return { 	"vertical binned spectrum" 			: cls.VertBinned,
					"smooth vertical binned spectrum" 	: cls.VertBinnedSmooth,
					"total intensity" 					: cls.Sum,
					"smooth total intensity"			: cls.SmoothSum }
	
	@classmethod
	def Select (cls, key) :
		return cls._GetDict()[key]
	
	@classmethod
	def GetChoices (cls) :
		return cls._GetDict().keys()

########################################################################

class ODD_Tab_FilterWheel (HardwareGUIControl) :
	"""
	GUI for ODD algorithm 
	"""
	def __init__ (self, parent) :
	
		HardwareGUIControl.__init__(self, parent)
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		######################################################################
		# List of positions of channels
		sizer.Add (wx.StaticText(self, label="Channel number with pure samples for learning"), flag=wx.LEFT, border=5)
		self.chanel_odd_experiment_ctrl = wx.TextCtrl (self, value="", style=wx.TE_MULTILINE|wx.EXPAND)
		self.chanel_odd_experiment_ctrl.__label__ = "channels"
		sizer.Add (self.chanel_odd_experiment_ctrl, flag=wx.EXPAND, border=5)
		
		# List of filters to be used in the motorized filter wheel
		sizer.Add (wx.StaticText(self, label="\nOptical filters in motorized wheel"), flag=wx.LEFT, border=5)
		self.filters_odd_experiment_ctrl = wx.TextCtrl (self, value="", style=wx.TE_MULTILINE|wx.EXPAND)
		self.filters_odd_experiment_ctrl.__label__ = "filters"
		sizer.Add (self.filters_odd_experiment_ctrl, flag=wx.EXPAND, border=5)
		
		# List of positions of channels_mixtures
		sizer.Add (wx.StaticText(self, label="\nChannel number to measure concentrations"), flag=wx.LEFT, border=5)
		self.channel_mixtues_ctrl = wx.TextCtrl (self, value="", style=wx.TE_MULTILINE|wx.EXPAND)
		self.channel_mixtues_ctrl.__label__ = "channels_mixtures"
		sizer.Add (self.channel_mixtues_ctrl, flag=wx.EXPAND, border=5)
		
		# Spectrum post-processing options
		sizer.Add ( wx.StaticText(self, label="\nSpectrum post-processing options") )
		choices = _SpectrumPostProcess.GetChoices()
		spectrum_options_ctrl = wx.ComboBox (self, choices=choices, value=choices[0], style=wx.CB_READONLY )
		spectrum_options_ctrl.__label__ = "spectrum_postprocess"
		sizer.Add (spectrum_options_ctrl,  flag=wx.EXPAND, border=5)	
		
		# Separator
		sizer.Add (wx.StaticText(self), border=5)
		
		# Record background signal
		background_signal_button = wx.Button (self, label="Record background")
		background_signal_button.Bind ( wx.EVT_BUTTON, self.RecordBackground )
		sizer.Add (background_signal_button, flag=wx.EXPAND, border=5)
		
		# Number of measurements taken sequentially
		sizer.Add (wx.StaticText(self, label="\nNumber of sequential measurements"), flag=wx.LEFT, border=5)
		num_seq_measurments_ctrl = wx.SpinCtrl (self, value="1", min=1)
		num_seq_measurments_ctrl.__label__ = "num_seq_measurments"
		sizer.Add (num_seq_measurments_ctrl, flag=wx.EXPAND, border=5)
		
		# Record spectra manually
		manualy_record_spectra_button = wx.Button (self, label="Manually record spectra")
		manualy_record_spectra_button.Bind( wx.EVT_BUTTON, self.ManualyRecordSpectra) 
		sizer.Add (manualy_record_spectra_button, flag=wx.EXPAND, border=5)
		
		# Measure concentration button
		measure_concentr_button = wx.Button (self)  
		measure_concentr_button._start_label 	= "Measure concentration"
		measure_concentr_button._start_method 	= self.MeasureConcentration
		measure_concentr_button._stop_label 	= "STOP measuring concentration"
		measure_concentr_button._stop_method 	= self.Stop_MeasureConcentration
		measure_concentr_button.SetLabel (measure_concentr_button._start_label)
		measure_concentr_button.Bind (wx.EVT_BUTTON, measure_concentr_button._start_method)
		sizer.Add(measure_concentr_button, flag=wx.EXPAND, border=5) 
		
		######################################################################
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		self.CreateSettingsDict()
	
	def RecordBackground (self, event=None) :
		"""
		Record background spectrum
		"""
		# Create pseudonyms 
		self.DevSpectrometer 	= self.parent.Spectrometer.dev
		self.DevSampleSwitcher	= self.parent.SampleSwitcher.dev
		
		# Initiate spectrometer
		settings = self.parent.Spectrometer.GetSettings()
		if self.DevSpectrometer.SetSettings(settings) == RETURN_FAIL : return
		
		# Initiate sample switcher
		settings = self.parent.SampleSwitcher.GetSettings()
		if self.DevSampleSwitcher.Initialize(settings) == RETURN_FAIL : return
		
		# Saving the name of channels 
		ga_settings = self.GetSettings()
		self.channels = sorted(eval( "(%s,)" % ga_settings["channels"] ))
		if self.DevSampleSwitcher.GetChannelNum()-1 < max(self.channels) :
			raise ValueError ("Error: Some channels specified are not accessible by sample switcher.")
			
		# Record background for each available channels 
		self.background_signal = {}
		for channel in self.channels :
			self.DevSampleSwitcher.MoveToChannel(channel)
			self.background_signal[ channel ] = self.DevSpectrometer.AcquiredData().astype(np.int)
	
	def ResetLogs (self) : 
		"""
		Reset all logs when new experiment begins
		"""
		self.optimization_log 		= { }
		self.reference_signal 		= { }
		self.log_reference_signal	= { }
		self.emission_spectra 		= { }
		self.N_emission_spectra		= { }
	
	def GetSampleSpectrum (self, channel) :
		"""
		Measure sample fluorescence spectra 
		"""
		# Get spectra
		spectrum = self.DevSpectrometer.AcquiredData().astype(np.int)
		
		# The following block is to obtain a super long time averaged 
		# emission spectra of molecules in each channel
		
		# The mean is calculated iteratively 
		# see, e.g., http://www.heikohoffmann.de/htmlthesis/node134.html
		try :
			self.N_emission_spectra[channel] += 1
			self.emission_spectra[channel] += ( spectrum - self.emission_spectra[channel] )/ self.N_emission_spectra[channel]
		except KeyError :
			self.emission_spectra[channel] = spectrum.astype(np.float)
			self.N_emission_spectra[channel] = 1
			
		# Subtract the background
		spectrum  -= self.background_signal[channel]
		
		return self.SpectrumPostProcess(spectrum)
		
	def CheckBackground (self) :
		"""
		Check whether the background signal is recorded and ready to be used. 
		"""
		try : 
			# The background signal must be consistent with self.channels
			for channel in self.channels :
				if channel not in self.background_signal : raise AttributeError
		except AttributeError :
			def SetBackgroundZero () :
				self.background_signal = dict( (channel, 0) for channel in self.channels ) 
				
			options = { "record background now" : self.RecordBackground, 
						"continue optimization without recording background" : SetBackgroundZero }
						
			dlg = wx.SingleChoiceDialog (self, 'Background sygnal has not been recoreded. Select one of the following option', 
				'Background signal not found', options.keys(), wx.CHOICEDLG_STYLE ) 
			dlg.Center()
			if dlg.ShowModal() == wx.ID_OK :
				options[ dlg.GetStringSelection() ]()
			else :
				# user cancel
				return
	
	def MeasureSpectra (self) :
		"""
		Measure spectra
		"""
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return
		
		print "Measuring spectra\n"
			
		# Initializing the dictionary where all measurements will be saved
		self.measurements = dict( itertools.product(self.channels, [{}] ) )
			
		# Loop over samples
		for channel in self.channels :
			# Go to channel
			self.DevSampleSwitcher.MoveToChannel(channel)
			
			# Loop over pulse shapes that are given by the filter
			for filter in self.filters :
			
				# Change filter
				self.FilterWheel.SetFilter( filter )
				
				wx.Yield()
				# abort, if requested 
				if self.need_abort : return
		
				# Acquiring data
				spectrum = self.GetSampleSpectrum(channel)
				self.measurements[channel][filter] = spectrum
		
				print "(channel = %d, filter = %d) Total fluorescence :  %1.3e" %(channel, filter, spectrum.sum())
			print "\n"
			
	def MeasureConcentration (self, event) :
		"""
		Perform concentration measurement 
		"""
		# Create pseudonyms of necessary devices 
		self.DevSpectrometer 	= self.parent.Spectrometer.dev
		self.FilterWheel		= self.parent.FilterWheel.dev
		self.DevSampleSwitcher	= self.parent.SampleSwitcher.dev
		
		# Save global settings and get the name of log file
		measurement_filename = SaveSettings(SettingsNotebook=self.parent, 
								title="Select file to save concentration measurements", filename="concentrations.hdf5")
		if measurement_filename is None : return
		
		####################### Initiate devices #############################
		
		# Initiate spectrometer
		settings = self.parent.Spectrometer.GetSettings()
		if self.DevSpectrometer.SetSettings(settings) == RETURN_FAIL : return
		
		# Initiate filter wheel
		settings = self.parent.FilterWheel.GetSettings()
		if self.FilterWheel.Initialize(settings) == RETURN_FAIL : return
		
		# Initiate sample switcher
		settings = self.parent.SampleSwitcher.GetSettings()
		if self.DevSampleSwitcher.Initialize(settings) == RETURN_FAIL : return
		
		settings = self.GetSettings()
		
		# Load channels with mixtures 
		channels_mixtures 	= sorted( eval( "(%s,)" % settings["channels_mixtures"] ) )
		channels_pure		= sorted( eval( "(%s,)" % settings["channels"] ) )
		
		# Saving the name of all channels containing pure and mixed samples
		self.channels = sorted(set( channels_pure + channels_mixtures  ))
		if self.DevSampleSwitcher.GetChannelNum()-1 < max(self.channels) :
			raise ValueError ("Error: Some channels specified are not accessible by sample switcher.")
		
		# Load filters
		self.filters	= sorted( eval( "(%s,)" % settings["filters"] ) ) 
		if min(self.filters) < 1 or max(self.filters) > self.FilterWheel.GetNumFilters() :
			raise ValueError ("Error: Some filters are not accessible by motorized wheel")
		
		# Set post-processing spectrum function 
		self.SpectrumPostProcess = _SpectrumPostProcess.Select( settings["spectrum_postprocess"] )
		
		# Check whether the background signal array is present
		self.CheckBackground()
		
		# Adjusting button's settings
		self.need_abort = False
		
		button = event.GetEventObject()
		button.SetLabel (button._stop_label)
		button.Bind( wx.EVT_BUTTON, button._stop_method)
		button.SetBackgroundColour('red')

		################ Concentration determination ################
		
		# Measure spectra from all channels 
		self.ResetLogs()
		self.MeasureSpectra()
		
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return 
		
		# Form fluorescent matrix of pure samples 
		fluoresence_matrix = np.vstack( [
			np.hstack( [self.measurements[channel][filter].flat for filter in self.filters] ) for channel in channels_pure  
		] ).T
		
		# This is for debugging purposes
		self.MeasureSpectra()
		
		# Calculate concentrations and save the results 
		with h5py.File (measurement_filename, 'a') as measurement_file :
			
			#################### Saving log information ####################
			# Delete existing summary group
			try : del measurement_file["optimization_summary"] 
			except KeyError : pass
			
			# Create new group
			summary_group = measurement_file.create_group("optimization_summary")
			for key, val in self.optimization_log.iteritems() :
				summary_group[key] = val
			"""
			# Save the reference signal summary 
			try : del measurement_file["reference_signal_summary"] 
			except KeyError : pass
			
			summary_group = measurement_file.create_group("reference_signal_summary")
			for key, val in self.log_reference_signal.iteritems() :
				summary_group[ "channel_%d" % key ] = val
			"""	
			# Save super long time averages of emission spectra  
			try : del measurement_file["emission_spectra"] 
			except KeyError : pass
			
			emission_spectra_grp = measurement_file.create_group("emission_spectra")
			for key, val in self.emission_spectra.iteritems() :
				emission_spectra_grp[ "channel_%d" % key ] = val
			
			##################### Saving raw measurements ######################
			try : del measurement_file["raw_measurements"] 
			except KeyError : pass
			
			raw_measurements_grp = measurement_file.create_group("raw_measurements")
			for channel, value in self.measurements.items() :
				channel_grp = raw_measurements_grp.create_group("channel_%d" % channel)
				for filter, spectrum in value.items() :
					channel_grp["filter_%d" % filter ] = spectrum
		
			##########################################################################
			
			# Create grope where results of concentration measurements will be saved 
			try : del measurement_file["concentration_measurements"]
			except KeyError : pass
			concentration_measurements_grp = measurement_file.create_group("concentration_measurements")
			
			# Perform calculations
			for channel in channels_mixtures :
				# Create group for individual channel containing mixtures 
				channel_grp = concentration_measurements_grp.create_group("channel_%d" % channel)
			
				# Grope fluorescence spectral data for the studied mixture 
				mixture_fluorescence = np.hstack( [ self.measurements[channel][filter].flat for filter in self.filters ] )
			
				# least square solution using all measured fluorescence data
				channel_grp["least square concentrations"] = \
					np.linalg.lstsq( fluoresence_matrix, mixture_fluorescence )[0]
				
				# non-negative least square solution
				channel_grp["non-negative least square concentrations"] = \
					nnls( fluoresence_matrix, mixture_fluorescence )[0]
					
				# print the results
				print "\nResults of concentration determination for channel %d" % channel
				for method, results in channel_grp.iteritems() :
					print "\t %s: \n\t\t %s" % (method, "; ".join("%.2e" % r for r in results[...]) )
		
		# Adjusting button's settings
		button = event.GetEventObject()
		button.SetLabel (button._start_label)
		button.Bind( wx.EVT_BUTTON, button._start_method)
		button.SetBackgroundColour('')
	
	def Stop_MeasureConcentration (self, event) :
		"""
		Stop measuring concentration 
		"""
		self.need_abort = True
		
		# Adjusting button's settings
		button = event.GetEventObject()
		button.SetLabel (button._start_label)
		button.Bind( wx.EVT_BUTTON, button._start_method)
		button.SetBackgroundColour('')
		
	def ManualyRecordSpectra (self, event) :
		"""
		Manually record spectra
		"""
		# Create pseudonyms of necessary devices 
		self.DevSpectrometer 	= self.parent.Spectrometer.dev
		self.DevSampleSwitcher	= self.parent.SampleSwitcher.dev
		
		# Save global settings and get the file name
		measurement_filename = SaveSettings(SettingsNotebook=self.parent, 
				title="Select file to save manual measurements of spectra", filename="manual_measruments.hdf5")
		if measurement_filename is None : return
		
		####################### Initiate devices #############################
		
		# Initiate spectrometer
		settings = self.parent.Spectrometer.GetSettings()
		if self.DevSpectrometer.SetSettings(settings) == RETURN_FAIL : return
		
		# Initiate sample switcher
		settings = self.parent.SampleSwitcher.GetSettings()
		if self.DevSampleSwitcher.Initialize(settings) == RETURN_FAIL : return
		
		settings = self.GetSettings()
		num_seq_measurments = settings["num_seq_measurments"]
		
		# Load channels 
		channels_mixtures 	= sorted( eval( "(%s,)" % settings["channels_mixtures"] ) )
		channels_pure		= sorted( eval( "(%s,)" % settings["channels"] ) )
		
		# Saving the name of all channels containing pure and mixed samples
		self.channels = sorted(set( channels_pure + channels_mixtures  ))
		if self.DevSampleSwitcher.GetChannelNum()-1 < max(self.channels) :
			raise ValueError ("Error: Some channels specified are not accessible by sample switcher.")
		
		# Set post-processing spectrum function 
		self.SpectrumPostProcess = _SpectrumPostProcess.Select( settings["spectrum_postprocess"] )
		
		# Check whether the background signal array is present
		self.CheckBackground()
		self.ResetLogs()
		
		# Save all measurements in file 
		with h5py.File (measurement_filename, 'a') as measurement_file :
			
			try : del measurement_file["manually_measured_spectra"]
			except KeyError : pass
			spectra_grp = measurement_file.create_group("manually_measured_spectra")

			# Create group for each C
			enumerated_channel_grp = [ (C, spectra_grp.create_group("channel_%d" % C)) for C in self.channels ]
				
			while True :
				# Ask user how to label this measurement
				dlg = wx.TextEntryDialog (self, "Enter the measurement label", caption="Select label for spectrum")
				if dlg.ShowModal() != wx.ID_OK : break
				label = dlg.GetValue()
				
				for channel, channel_grp in enumerated_channel_grp :
					if label in channel_grp :
						wx.MessageDialog(self, "This measurements is ignored because the label has already been used", 
												caption="Signal not recorder").ShowModal() 
						break
					else :
						# Go to channel
						self.DevSampleSwitcher.MoveToChannel(channel)
						
						measurments_seq_grp = channel_grp.create_group( label )
						for num in xrange(num_seq_measurments) :
							# Save measurements
							spectrum = self.GetSampleSpectrum (channel)
							measurments_seq_grp[ str(num) ] = spectrum
							print "Channel %d / %s: total fluorescence %1.3e" % (channel, label, spectrum.sum()) 	
					print "\n"
		
			
			"""
			for channel in self.channels :
				# Go to channel
				self.DevSampleSwitcher.MoveToChannel(channel)
				channel_grp = spectra_grp.create_group("channel_%d" % channel)
				
				while True :
					# Ask user how to label this measurement
					dlg = wx.TextEntryDialog (self, "Enter label for the current measurement of channel %d" % channel, 
												caption="Select label for spectrum")
					dlg.Center()
					if dlg.ShowModal() != wx.ID_OK: break
					else :
						label = dlg.GetValue()
						if label in channel_grp :
							wx.MessageDialog(self, "This measurements is ignored because the label has already been used", 
												caption="Signal not recorder").ShowModal() 
						else :
							# Save measurements
							spectrum = self.GetSampleSpectrum (channel)
							channel_grp[ label ] = spectrum
							print "Channel %d / %s: total fluorescence %1.3e" % (channel, label, spectrum.sum()) 
			
				print "\n"
			"""
			#################### Saving log information ####################
			# Save super long time averages of emission spectra  
			try : del measurement_file["emission_spectra"] 
			except KeyError : pass
			
			emission_spectra_grp = measurement_file.create_group("emission_spectra")
			for key, val in self.emission_spectra.iteritems() :
				emission_spectra_grp[ "channel_%d" % key ] = val