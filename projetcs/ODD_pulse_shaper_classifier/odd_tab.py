########################################################################
#
#	Scanning chirps 
#
########################################################################

import visvis 
import wx
import h5py
from wx.lib.agw.floatspin import FloatSpin as wxFloatSpin
import numpy as np
from scipy.ndimage.filters import gaussian_filter

from libs.gui.hardware_control import HardwareGUIControl
from libs.gui.basic_window import SaveSettings 

from libs.dev.consts import * 

########################################################################

class ODD_Tab (HardwareGUIControl) :
	"""
	GUI to scan chirp
	"""
	def __init__ (self, parent) :
		HardwareGUIControl.__init__(self, parent)
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		# List of positions of channels
		sizer.Add (wx.StaticText(self, label="Channel number with pure samples for learning"), flag=wx.LEFT, border=5)
		self.chanel_odd_experiment_ctrl = wx.TextCtrl (self, value="", style=wx.TE_MULTILINE|wx.EXPAND)
		self.chanel_odd_experiment_ctrl.__label__ = "channels"
		sizer.Add (self.chanel_odd_experiment_ctrl, flag=wx.EXPAND, border=5)
	
		sizer.Add (wx.StaticText(self, label="\nMax amplitude (ND filter)"), flag=wx.LEFT, border=5)
		max_ampl_ctrl = wxFloatSpin(self, min_val=0, max_val=1, increment=0.01, value=1., digits=3)
		max_ampl_ctrl.__label__ = "max_ampl"
		sizer.Add (max_ampl_ctrl , flag=wx.EXPAND, border=5)
	
		################ Parameters of reference mask ####################
		sb_sizer = wx.StaticBoxSizer( wx.StaticBox(self, label="Reference mask parameters: "),  wx.VERTICAL )
		
		# Min value of coefficient
		sb_sizer.Add (wx.StaticText(self, label="min coefficient"), flag=wx.LEFT, border=5)
		coeff_min_ctrl = wxFloatSpin (self, min_val=-10, max_val=10, increment=0.01, value=-0.9, digits=3)
		coeff_min_ctrl.__label__ = "coeff_min"
		sb_sizer.Add (coeff_min_ctrl , flag=wx.EXPAND, border=5)
		
		# Max value of coefficient
		sb_sizer.Add (wx.StaticText(self, label="max coefficient"), flag=wx.LEFT, border=5)
		coeff_max_ctrl = wxFloatSpin (self, min_val=0, max_val=10, increment=0.01, value=0.9, digits=3)
		coeff_max_ctrl.__label__ = "coeff_max"
		sb_sizer.Add (coeff_max_ctrl , flag=wx.EXPAND, border=5)
		
		# Number of coeff scans
		sb_sizer.Add (wx.StaticText(self, label="number of scans"), flag=wx.LEFT, border=5)
		coeff_num_ctrl = wx.SpinCtrl (self, value="20", min=10, max=100000)
		coeff_num_ctrl.__label__ = "coeff_num"
		sb_sizer.Add (coeff_num_ctrl , flag=wx.EXPAND, border=5)
		
		# Polynomial basis type
		sb_sizer.Add (wx.StaticText(self, label="\nPolynomial basis type"), flag=wx.LEFT, border=5)
		self.polynomial_bases = {
			"Chebyshev"	:	np.polynomial.chebyshev.Chebyshev,
			"Legendre"	:	np.polynomial.legendre.Legendre,
			"Laguerre"	:	np.polynomial.laguerre.Laguerre,
			"Hermite"	:	np.polynomial.hermite.Hermite,
			"Monomials" :	np.polynomial.polynomial.Polynomial
		}  
		choices = self.polynomial_bases.keys()
		polynomial_bais_ctrl = wx.ComboBox (self, choices=choices, value=choices[0], style=wx.CB_READONLY )
		polynomial_bais_ctrl.__label__ = "polynomial_basis"
		sb_sizer.Add (polynomial_bais_ctrl,  flag=wx.EXPAND, border=5)
		
		# Polynomial order
		sb_sizer.Add (wx.StaticText(self, label="\npolynomial order"), flag=wx.LEFT, border=5)
		poly_order_ctrl = wx.SpinCtrl (self, value="2", min=0, max=100000)
		poly_order_ctrl.__label__ = "polynomial_order"
		sb_sizer.Add (poly_order_ctrl , flag=wx.EXPAND, border=5)
		
		sizer.Add (sb_sizer, flag=wx.EXPAND, border=10)
		################################################
		
		# Scan button
		self.get_coeff_scan_btn = wx.Button (self)
		self.get_coeff_scan_btn._start_label  = "Scan polynomial coefficient"
		self.get_coeff_scan_btn._start_method = self.DoScannning
		self.get_coeff_scan_btn._stop_label  = "STOP scanning"
		self.get_coeff_scan_btn._stop_method = self.StopScannning
		self.get_coeff_scan_btn.SetLabel( self.get_coeff_scan_btn._start_label )
		self.get_coeff_scan_btn.Bind( wx.EVT_BUTTON, self.get_coeff_scan_btn._start_method )
		sizer.Add (self.get_coeff_scan_btn, flag=wx.EXPAND, border=5)
		
		# Separator
		sizer.Add (wx.StaticText(self), border=5)
		
		# Record background signal
		background_signal_button = wx.Button (self, label="Record background")
		background_signal_button.Bind ( wx.EVT_BUTTON, self.RecordBackground )
		sizer.Add (background_signal_button, flag=wx.EXPAND, border=5)
		
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
		odd_settings = self.GetSettings()
		self.channels = sorted(eval( "(%s,)" % odd_settings["channels"] ))
		if self.DevSampleSwitcher.GetChannelNum()-1 < max(self.channels) :
			raise ValueError ("Error: Some channels specified are not accessible by sample switcher.")
			
		# Record background for each available channels 
		self.background_signal = {}
		for channel in self.channels :
			self.DevSampleSwitcher.MoveToChannel(channel)
			self.background_signal[ channel ] = self.DevSpectrometer.AcquiredData().astype(np.int)
	
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
						"continue without recording background" : SetBackgroundZero }
						
			dlg = wx.SingleChoiceDialog (self, 'Background sygnal has not been recoreded. Select one of the following option', 
				'Background signal not found', options.keys(), wx.CHOICEDLG_STYLE ) 
			
			if dlg.ShowModal() == wx.ID_OK :
				options[ dlg.GetStringSelection() ]()
			else :
				# user cancel
				return
	
	def GetSampleSpectrum (self, channel) :
		"""
		Measure sample fluorescence spectra 
		"""
		# Get spectra
		spectrum = self.DevSpectrometer.AcquiredData().astype(np.int)
		
		"""
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
		"""		
		# Subtract the background
		spectrum  -= self.background_signal[channel]
		
		#return self.SpectrumPostProcess(spectrum)
		return spectrum
		
	def DoScannning (self, event) :
		"""
		Perform scanning of different phase mask 
		"""
		# Create pseudonyms of necessary devices 
		self.DevSpectrometer 	= self.parent.Spectrometer.dev
		self.DevPulseShaper		= self.parent.PulseShaper.dev
		self.DevSampleSwitcher	= self.parent.SampleSwitcher.dev
			
		# Save global settings and get the name of log file
		self.log_filename = SaveSettings(SettingsNotebook=self.parent, 
								title="Select file to save phase mask scanning", filename="scanning_phase_mask.hdf5")
		if self.log_filename is None : return
		
		####################### Initiate devices #############################
		
		# Initiate spectrometer
		settings = self.parent.Spectrometer.GetSettings()
		if self.DevSpectrometer.SetSettings(settings) == RETURN_FAIL : return
		
		# Initiate pulse shaper
		settings = self.parent.PulseShaper.GetSettings()
		if self.DevPulseShaper.Initialize(settings) == RETURN_FAIL : return
		
		# Get number of optimization variables 
		self.num_pixels = self.DevPulseShaper.GetParamNumber()
		if self.num_pixels == RETURN_FAIL :
			raise RuntimeError ("Optimization cannot be started since calibration file was not loaded")
		
		# Initiate sample switcher
		settings = self.parent.SampleSwitcher.GetSettings()
		if self.DevSampleSwitcher.Initialize(settings) == RETURN_FAIL : return
		
		# Saving the name of channels 
		odd_settings = self.GetSettings()
		self.channels = sorted(eval( "(%s,)" % odd_settings["channels"] ))
		if self.DevSampleSwitcher.GetChannelNum()-1 < max(self.channels) :
			raise ValueError ("Error: Some channels specified are not accessible by sample switcher.")
		
		# Check whether the background signal array is present
		self.CheckBackground()
		
		#####################################################################
		
		# Get range of coefficient 
		coeff_range = np.linspace( odd_settings["coeff_min"], odd_settings["coeff_max"], odd_settings["coeff_num"] )
		
		# List all polynomial coefficients
		N = odd_settings["polynomial_order"]
		poly_coeffs = np.zeros( (coeff_range.size*N, N+1) )
		for n in range(1,N+1) :
			poly_coeffs[(n-1)*coeff_range.size:n*coeff_range.size, n ] = coeff_range
		
		# Chose max amplitude
		max_ampl = odd_settings["max_ampl"]*np.ones(self.num_pixels)
		
		# Arguments of the basis 
		X = np.linspace(-1., 1., self.num_pixels)
		
		# Retrieve the basis type
		polynomial_basis = self.polynomial_bases[ odd_settings["polynomial_basis"] ]
		
		# Adjusting button's settings
		button = event.GetEventObject()
		button.SetLabel (button._stop_label)
		button.SetBackgroundColour('red')
		button.Bind( wx.EVT_BUTTON, button._stop_method)
		self.need_abort = False
		
		#####################################################################
		
		# Start scanning
		with  h5py.File (self.log_filename, 'a') as log_file :
			for channel in self.channels :
				# Move to a selected channel 
				self.DevSampleSwitcher.MoveToChannel(channel)
		
				# abort, if requested 
				wx.Yield()
				if self.need_abort : break
				
				# Looping over pulse shapes
				for scan_num, coeff in enumerate(poly_coeffs) :
				
					# Calculate new phase
					phase = polynomial_basis(coeff)(X)
				
					# Set the pulse shape
					self.DevPulseShaper.SetAmplPhase(max_ampl, phase) 
		
					# Save phase in radians 
					ampl, phase = self.DevPulseShaper.GetUnwrappedAmplPhase(max_ampl, phase)
					if scan_num == 0 :
						# Initialize the array
						phases_rad = np.zeros( (len(poly_coeffs), phase.size), dtype=phase.dtype )
						amplitudes = np.zeros_like(phases_rad)
						amplitudes[:] = ampl.min()
					
					# Save phase
					phases_rad[scan_num] = phase
					amplitudes[scan_num] = ampl
					
					# abort, if requested 
					wx.Yield()
					if self.need_abort : break
					
					# Get spectrum
					spectrum = self.GetSampleSpectrum(channel)
					# Vertical binning
					spectrum = ( spectrum.sum(axis=0) if len(spectrum.shape) == 2 else spectrum )
					
					if scan_num == 0 :
						# Initialize the array
						spectra = np.zeros( (len(poly_coeffs), spectrum.size), dtype=spectrum.dtype )
						spectra[:] = spectrum.min()
						
					# Save the spectrum
					spectra[scan_num] = spectrum
					
					# Display the currently acquired data
					try : 
						spectra_2d_img.SetData(spectra)
						phases_rad_2d_img.SetData( phases_rad%(2*np.pi) )
						#amplitudes_2d_img.SetData( amplitudes )
					except NameError :
						visvis.cla(); visvis.clf()
						
						visvis.subplot(121)
						spectra_2d_img = visvis.imshow(spectra, cm=visvis.CM_JET)
						visvis.ylabel ('scans'); visvis.xlabel ('wavelegth')
						visvis.title("spectral scan")
						
						visvis.subplot(122)
						phases_rad_2d_img = visvis.imshow( phases_rad%(2*np.pi), cm=visvis.CM_JET)
						visvis.title("phase shapes")
						
						#visvis.subplot(133)
						#amplitudes_2d_img = visvis.imshow(amplitudes, cm=visvis.CM_JET)
						#visvis.title("amplitudes")
						
				# Save the data for the given channel
				try : del log_file[ str(channel) ]
				except KeyError : pass
				
				channel_grp = log_file.create_group( str(channel) )
				channel_grp["spectra"] 		= spectra
				channel_grp["phases_rad"]	= phases_rad
				channel_grp["amplitudes"]	= amplitudes
				channel_grp["poly_coeffs"]	= poly_coeffs
				
		# Readjust buttons settings
		self.StopScannning(event)
		
	def StopScannning (self, event) :
		"""
		Stop scanning
		"""
		self.need_abort = True 
		
		# Adjusting button's settings
		button = event.GetEventObject()
		button.SetLabel (button._start_label)
		button.SetBackgroundColour('')
		button.Bind( wx.EVT_BUTTON, button._start_method)