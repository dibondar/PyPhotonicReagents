########################################################################
#
#	Implementation of Multiphoton Intrapulse Interference Phase Scan (MIIPS)
#	Vadim V. Lozovoy, Igor Pastirk, and Marcos Dantus, Optics Letters 29, 775 (2004) 
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


class MIIPS_Tab (HardwareGUIControl) :
	"""
	GUI for Multiphoton Intrapulse Interference Phase Scan (MIIPS)
	"""
	def __init__ (self, parent) :
		
		HardwareGUIControl.__init__(self, parent)
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		################ Parameters of reference mask ####################
		
		# Hint
		sizer.Add (wx.StaticText(self, label="\nReference mask parameters: "), flag=wx.LEFT, border=5)
		
		# Min value of coefficient
		sizer.Add (wx.StaticText(self, label="\nmin coefficient"), flag=wx.LEFT, border=5)
		coeff_min_ctrl = wxFloatSpin (self, min_val=-10, max_val=10, increment=0.01, value=-0.9, digits=3)
		coeff_min_ctrl.__label__ = "coeff_min"
		sizer.Add (coeff_min_ctrl , flag=wx.EXPAND, border=5)
		
		# Max value of coefficient
		sizer.Add (wx.StaticText(self, label="max coefficient"), flag=wx.LEFT, border=5)
		coeff_max_ctrl = wxFloatSpin (self, min_val=0, max_val=10, increment=0.01, value=0.9, digits=3)
		coeff_max_ctrl.__label__ = "coeff_max"
		sizer.Add (coeff_max_ctrl , flag=wx.EXPAND, border=5)
		
		# Number of coeff scans
		sizer.Add (wx.StaticText(self, label="number of scans"), flag=wx.LEFT, border=5)
		coeff_num_ctrl = wx.SpinCtrl (self, value="10", min=1, max=100000)
		coeff_num_ctrl.__label__ = "coeff_num"
		sizer.Add (coeff_num_ctrl , flag=wx.EXPAND, border=5)
		
		# Polynomial basis type
		sizer.Add (wx.StaticText(self, label="\nPolynomial basis type"), flag=wx.LEFT, border=5)
		self.polynomial_bases = {
			"Chebyshev"	:	np.polynomial.chebyshev.chebval,
			"Legendre"	:	np.polynomial.legendre.legval,
			"Monomials" :	np.polyval
		}  
		choices = self.polynomial_bases.keys()
		polynomial_bais_ctrl = wx.ComboBox (self, choices=choices, value=choices[0], style=wx.CB_READONLY )
		polynomial_bais_ctrl.__label__ = "polynomial_basis"
		sizer.Add (polynomial_bais_ctrl,  flag=wx.EXPAND, border=5)
		
		# Polynomial order
		sizer.Add (wx.StaticText(self, label="\npolynomial order"), flag=wx.LEFT, border=5)
		poly_order_ctrl = wx.SpinCtrl (self, value="2", min=0, max=100000)
		poly_order_ctrl.__label__ = "polynomial_order"
		sizer.Add (poly_order_ctrl , flag=wx.EXPAND, border=5)
		
		# Scan button
		get_coeff_scan = wx.Button (self, label="Scan polynomial coefficient")
		get_coeff_scan.Bind( wx.EVT_BUTTON, self.GetDeltaScan )
		sizer.Add (get_coeff_scan , flag=wx.EXPAND, border=5)
		
		#######################################################
		
		# Separator
		sizer.Add (wx.StaticText(self), border=5)
		
		# Select current value of coefficient 
		sizer.Add (wx.StaticText(self, label="Value of current coefficient"), flag=wx.LEFT, border=5)
		current_coeff_val = wxFloatSpin (self, min_val=-10, max_val=10, increment=0.01, value=0., digits=3)
		current_coeff_val.__label__ = "current_coeff_val"
		sizer.Add (current_coeff_val , flag=wx.EXPAND, border=5)
		
		# Button to add the correction
		correct_phase = wx.Button (self, label="Add correction")
		correct_phase.Bind( wx.EVT_BUTTON, self.AddCurrentCorrection )
		sizer.Add (correct_phase, flag=wx.EXPAND, border=5)
		
		# Separator
		sizer.Add (wx.StaticText(self), border=5)
		
		# Buttons to save correcting procedure
		save_correction_log = wx.Button (self, label="Save correction log")
		save_correction_log.Bind( wx.EVT_BUTTON, self.SaveCorrectionLog )
		sizer.Add (save_correction_log, flag=wx.EXPAND, border=5)
		
		# Buttons to plot reference phase
		plot_reference_phase = wx.Button (self, label="Plot reference phase")
		plot_reference_phase.Bind( wx.EVT_BUTTON, self.PlotReferencePhase  )
		sizer.Add (plot_reference_phase, flag=wx.EXPAND, border=5)
		
		# Buttons to copy reference phase to clipboard
		reference_phase2clipboard = wx.Button (self, label="Copy reference phase to clipboard")
		reference_phase2clipboard.Bind( wx.EVT_BUTTON, self.ReferencePhase2Clipboard  )
		sizer.Add (reference_phase2clipboard, flag=wx.EXPAND, border=5)
		
		# Separator
		sizer.Add (wx.StaticText(self), border=5)
		
		# Clear all button
		clear_all = wx.Button (self, label="Clear all")
		clear_all.Bind( wx.EVT_BUTTON, self.ClearAllData )
		sizer.Add (clear_all, flag=wx.EXPAND, border=5)
		
		# Clear all info
		self.ClearAllData (None)
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		self.CreateSettingsDict()
	
	def ReferencePhase2Clipboard (self, event) :
		"""
		Copy the values of the reference phase into the Clipboard as text 
		"""
		if wx.TheClipboard.Open() :
			txt_data = str( self.GetReferencePhase()[0] )
			wx.TheClipboard.SetData(wx.TextDataObject(txt_data)) 
			wx.TheClipboard.Close()
		else :
			wx.MessageBox("Unable to open the clipboard", "Error")
		
	def ClearAllData (self, event) :
		"""
		Reset all data
		"""
		# List of all scans 
		self.scan_log = []
		
		# List of corrections
		self.corrections = np.zeros(1) 
		
	def SaveCorrectionLog (self, event) :
		"""
		Save the log of corrections
		"""
		# Save global settings and get the name of log file
		self.miips_filename = SaveSettings(SettingsNotebook=self.parent, 
								title="Select file to save MIIPS scan", filename="miips_scan.hdf5")
		if self.miips_filename is None : return
		
		# Save scan into the file
		with h5py.File(self.miips_filename, 'a') as F :
			try : del F["scan_log"]
			except KeyError : pass
			
			# Save scanning log
			scan_log_grp = F.create_group("scan_log")
			for num, scan in enumerate(self.scan_log) :
				scan_num_grp = scan_log_grp.create_group( str(num) )
				for key, val in scan.iteritems() :
					scan_num_grp[key] = val
					
			# Save corrections
			try : del F["reference_phase"]
			except KeyError : pass
			reference_phase_grp = F.create_group("reference_phase")
			reference_phase_grp["corrections"] = self.corrections
			reference_phase_grp["reference_phase"] = self.GetReferencePhase()[0]
			
	def AddCurrentCorrection (self, event) :
		"""
		Add current correction
		"""
		miips_settings = self.GetSettings()
		
		N = miips_settings["polynomial_order"]
		
		# Pad array if necessary
		if self.corrections.size < N+1 :
			self.corrections = np.lib.pad(self.corrections, (0,N+1-self.corrections.size), 'constant', constant_values=0)
		
		coeff_val = 0.5*miips_settings["current_coeff_val"]
		self.corrections[N] += coeff_val
		if N :
			self.corrections[0] += coeff_val
	
		# plot the data
		self.PlotReferencePhase(None)
	
	def PlotReferencePhase (self, event) :
		"""
		Plot reference phase
		"""
		visvis.cla(); visvis.clf(); 
		
		visvis.subplot(211)
		visvis.plot( self.GetReferencePhase()[0] )
		visvis.ylabel ('Reference phase')
		visvis.xlabel ('puls shaper pixel number')
		visvis.title ('Current reference phase')
		
		visvis.subplot(212)
		visvis.plot( self.corrections, lc='r',ms='*', mc='r' )
		visvis.ylabel ('Value of coefficients')
		visvis.xlabel ('coefficient number')
		visvis.title ('Value of corrections')
		
	def GetReferencePhase (self, miips_settings=None) :
		"""
		Retrieve the phase correction 
		"""
		if miips_settings is None :
			miips_settings = self.GetSettings()
		
		# Retrieve the basis type
		polyval = self.polynomial_bases[ miips_settings["polynomial_basis"] ]
		
		# Arguments of the basis 
		X = np.linspace(-1., 1., self.num_pixels)
		
		# Get the phase correction
		reference_phase = polyval(X, self.corrections) 
	
		# Check for consistency
		np.clip(reference_phase, 0, 1, out=reference_phase)
		
		return reference_phase, X, polyval
		
	def GetDeltaScan (self, event) :
		"""
		Measure spectra by varying parameter delta in reference phase mask
		"""
		# Create pseudonyms of necessary devices 
		self.DevSpectrometer 	= self.parent.Spectrometer.dev
		self.DevPulseShaper		= self.parent.PulseShaper.dev
		
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
		
		miips_settings = self.GetSettings()
		#####################################################################
		
		# Get range of coefficient 
		coeff_range = np.linspace( miips_settings["coeff_min"], miips_settings["coeff_max"], miips_settings["coeff_num"] )
	
		# Get reference phase based on the corrections already obtained
		reference_phase, X, polyval = self.GetReferencePhase(miips_settings)
		
		# Get new polynomial
		new_coeffs = np.zeros(miips_settings["polynomial_order"] + 1)
		new_coeffs[-1] = 1
		current_polynomial = polyval(X, new_coeffs)
	
		# Chose max amplitude
		ones_ampl = np.ones(self.num_pixels)
		
		for scan_num, coeff in enumerate(coeff_range) :
			# Get current phase mask
			current_phase  = 0.5*coeff*( current_polynomial + 1 )
			current_phase += reference_phase
			
			# Check for consistency
			np.clip(current_phase, 0, 1, out=current_phase)
			
			# Set the pulse shape
			self.DevPulseShaper.SetAmplPhase(ones_ampl, current_phase) 
			
			# Make sure software is not frozen
			wx.Yield()
			
			# Get spectrum
			spectrum = self.DevSpectrometer.AcquiredData()
			
			# Vertical binning
			spectrum = ( spectrum.sum(axis=0) if len(spectrum.shape) == 2 else spectrum )
		
			try : spectral_scans
			except NameError :
				# Allocate space for spectral scans
				spectral_scans = np.zeros( (coeff_range.size, spectrum.size), dtype=spectrum.dtype )
			
			# Save spectrum
			spectral_scans[ scan_num ] = spectrum
			
			# Display the currently acquired data
			try : 
				spectral_scan_2d_img.SetData(spectral_scans)
				# spectral_scan_2d_img.SetClim( spectral_scans.min(), spectral_scans.max() )
			except NameError :
				visvis.cla(); visvis.clf(); 	
				spectral_scan_2d_img = visvis.imshow(spectral_scans, cm=visvis.CM_JET)
				visvis.ylabel ('coefficeints')
				visvis.xlabel ('wavelegth')
				
		#######################################################
		
		# Get current wavelength
		wavelength = self.DevSpectrometer.GetWavelengths()
		
		# Adding the scan to log
		self.scan_log.append( {
			"spectral_scans" 	:	spectral_scans,
			"reference_phase" 	:	reference_phase,
			"current_polynomial":	current_polynomial,
			"coeff_range"		: 	coeff_range,
			"wavelength"		:	wavelength
		} )
		
		# Determine an optimal value of coeff
		max_sum_val = np.argmax( gaussian_filter(spectral_scans, 2).sum(axis=1) )
		self.current_coeff_val.SetValue( coeff_range[max_sum_val] )
		spectral_scans[max_sum_val] = spectral_scans.min()
		
		################ Plotting results ####################
		spectral_scan_2d_img.SetData(spectral_scans)
		