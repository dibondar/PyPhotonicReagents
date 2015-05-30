########################################################################
#
#	Scanning chirps 
#
########################################################################

import visvis 
import wx
import wx.grid
import h5py
from itertools import repeat, product, chain, islice, cycle
from wx.lib.agw.floatspin import FloatSpin as wxFloatSpin
import numpy as np
from scipy.ndimage.filters import gaussian_filter

from libs.gui.hardware_control import HardwareGUIControl
from libs.gui.basic_window import SaveSettings 

from libs.dev.consts import * 

########################################################################

class CombinatorialChirpScan_Tab (HardwareGUIControl) :
	"""
	GUI to scan chirp combinatorially
	"""
	def __init__ (self, parent) :
		HardwareGUIControl.__init__(self, parent)
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		# List of positions of channels
		sizer.Add (wx.StaticText(self, label="Channel number with pure samples for learning"), flag=wx.LEFT, border=5)
		self.chanel_odd_experiment_ctrl = wx.TextCtrl (self, value="", style=wx.TE_MULTILINE|wx.EXPAND)
		self.chanel_odd_experiment_ctrl.__label__ = "channels"
		sizer.Add (self.chanel_odd_experiment_ctrl, flag=wx.EXPAND, border=5)
		
		# Filters (on filter wheel) to be used 
		sizer.Add (wx.StaticText(self, label="Filter numbers"), flag=wx.LEFT, border=5)
		self.filters_ctrl = wx.TextCtrl (self, value="", style=wx.TE_MULTILINE|wx.EXPAND)
		self.filters_ctrl.__label__ = "filters"
		sizer.Add (self.filters_ctrl, flag=wx.EXPAND, border=5)
	
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
		coeff_num_ctrl = wx.SpinCtrl (self, value="20", min=2, max=100000)
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
		
		# Initial polynomial order
		sb_sizer.Add (wx.StaticText(self, label="\nInitial polynomial order"), flag=wx.LEFT, border=5)
		min_poly_order_ctrl = wx.SpinCtrl (self, value="2", min=0, max=100000)
		min_poly_order_ctrl.__label__ = "min_polynomial_order"
		sb_sizer.Add (min_poly_order_ctrl , flag=wx.EXPAND, border=5)
		
		# Final polynomial order
		sb_sizer.Add (wx.StaticText(self, label="\nFinal polynomial order"), flag=wx.LEFT, border=5)
		max_poly_order_ctrl = wx.SpinCtrl (self, value="4", min=0, max=100000)
		max_poly_order_ctrl.__label__ = "max_polynomial_order"
		sb_sizer.Add (max_poly_order_ctrl , flag=wx.EXPAND, border=5)
		
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
		settings = self.GetSettings()
		self.channels = sorted(eval( "(%s,)" % settings["channels"] ))
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
		self.DevFilterWheel		= self.parent.FilterWheel.dev
		
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
		
		# Initialize the filter wheel
		settings = self.parent.FilterWheel.GetSettings()
		if self.DevFilterWheel.Initialize(settings) == RETURN_FAIL : return
		
		# Saving the name of channels 
		settings = self.GetSettings()
		self.channels = sorted( eval( "(%s,)" % settings["channels"] ) )
		if self.DevSampleSwitcher.GetChannelNum()-1 < max(self.channels) :
			raise ValueError ("Error: Some channels specified are not accessible by sample switcher.")
			
		# Saving the name of filter 
		self.filters = sorted( eval( "(%s,)" % settings["filters"] ) )
		if self.DevFilterWheel.GetNumFilters() < max(self.filters) :
			raise ValueError ("Error: Some filters specified are not accessible by filter wheel.")
		
		# Check whether the background signal array is present
		self.CheckBackground()
		
		#####################################################################
		
		# Get range of coefficient 
		coeff_range = np.linspace( settings["coeff_min"], settings["coeff_max"], settings["coeff_num"] )
		
		min_N = settings["min_polynomial_order"]
		max_N = settings["max_polynomial_order"]
		
		# Create all polynomial coefficients for scanning
		coeff_iter = product( *chain(repeat([0], min_N), repeat(coeff_range, max_N+1-min_N)) )
		poly_coeffs = np.array( list(coeff_iter) )
		
		# Chose max amplitude
		max_ampl = settings["max_ampl"]*np.ones(self.num_pixels)
		 
		# Arguments of the basis 
		X = np.linspace(-1., 1., self.num_pixels)
		
		# Retrieve the basis type
		polynomial_basis = self.polynomial_bases[ settings["polynomial_basis"] ]
		
		# Adjusting button's settings
		button = event.GetEventObject()
		button.SetLabel (button._stop_label)
		button.SetBackgroundColour('red')
		button.Bind( wx.EVT_BUTTON, button._stop_method)
		self.need_abort = False
		
		#####################################################################
		
		# Start scanning
		with  h5py.File (self.log_filename, 'a') as log_file :
		
			# Allocate memory for fluorescence signal of all proteins
			fluorescences =  dict( 
				(key, np.zeros(len(poly_coeffs), dtype=np.int)) for key in product(self.channels, self.filters) 
			)
			
			# Allocate memory for reference fluoresence
			ref_fluorescences = dict( (C, []) for C in self.channels ) 
			
			# Iterator for polynomial coefficients
			poly_coeffs_iter = enumerate(poly_coeffs)
			
			while True :
				# Chunk up the data
				chunk_poly_coeffs = list( islice(poly_coeffs_iter, len(coeff_range)) )
				
				# about, if there is nothing to iterate over
				if len(chunk_poly_coeffs) == 0 : break
				
				# abort, if requested
				if self.need_abort : break
		
				for channel in self.channels :
					# Move to a selected channel 
					self.DevSampleSwitcher.MoveToChannel(channel)
			
					# abort, if requested 
					if self.need_abort : break
					
					# Looping over pulse shapes
					for scan_num, coeff in chunk_poly_coeffs :
					
						# Calculate new phase
						phase = polynomial_basis(coeff)(X)
					
						# Set the pulse shape
						self.DevPulseShaper.SetAmplPhase(max_ampl, phase) 
						
						# abort, if requested 
						if self.need_abort : break
						
						# Looping over filters in filter wheels
						for filter in self.filters :
							self.DevFilterWheel.SetFilter(filter)
							
							# abort, if requested
							wx.Yield()
							if self.need_abort : break
							
							# Get spectrum sum, i.e., fluorescence
							spectrum_sum = self.GetSampleSpectrum(channel).sum()
							
							# Save the spectrum
							fluorescences[ (channel, filter) ][scan_num] = spectrum_sum
					
					# Record reference fluorescence 
					self.DevPulseShaper.SetAmplPhase(max_ampl, np.zeros_like(max_ampl)) 
					
					# Go to filters with max transmission (ASSUMING that it has lowest number) 
					self.DevFilterWheel.SetFilter( min(self.filters) )
					
					spectrum_sum = self.GetSampleSpectrum(channel).sum()
					ref_fluorescences[channel].append( spectrum_sum )
					
				# Display the currently acquired data
				try : 
					for key, img in fluorescence_imgs.items() :
						img.SetYdata( fluorescences[key] )
				except NameError :
					visvis.cla(); visvis.clf()
						
					# Iterator of colours
					colour_iter = cycle(['r', 'g', 'b', 'k', 'y'])
						
					fluorescence_imgs = dict( 
						(K, visvis.plot( F, lw=1, lc=colour_iter.next() ) ) for K, F  in fluorescences.items() 
					)
					visvis.xlabel("phase masks")
					visvis.ylabel("Integrated fluorescence")
					visvis.title("Measured fluorescence")
							
			################ Scanning is over, save the data ########################
			
			# Delete and create groups corresponding to channel
			for channel in self.channels :
				try : del log_file[ "channel_%d" % channel ]
				except KeyError : pass
				gchannel_grp = log_file.create_group( "channel_%d" % channel )
				gchannel_grp["ref_fluorescence"] = ref_fluorescences[channel]
				gchannel_grp["poly_coeffs"]	= poly_coeffs
				
			for channel_filter, fluorescence in  fluorescences.items() :	
				log_file["channel_%d/filter_%d_fluorescence" % channel_filter ]	= fluorescence
				
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