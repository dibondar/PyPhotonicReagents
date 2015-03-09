########################################################################
#
#	Implementation of Multiphoton Intrapulse Interference Phase Scan (MIIPS)
#	Vadim V. Lozovoy, Igor Pastirk, and Marcos Dantus, Optics Letters 29, 775 (2004) 
# 
########################################################################

import visvis 
import wx
from wx.lib.agw.floatspin import FloatSpin as wxFloatSpin
import numpy as np

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
		sizer.Add (wx.StaticText(self, label="\nReference mask parameters: \nf(w) = alpha*cos( gamma * w + delta )"), flag=wx.LEFT, border=5)
		
		# alpha parameter ctr
		sizer.Add (wx.StaticText(self, label="alpha = "), flag=wx.LEFT, border=5)
		alpha_ctrl = wxFloatSpin (self, min_val=0, max_val=1, increment=0.01, value=0.8, digits=3)
		alpha_ctrl.__label__ = "alpha"
		sizer.Add (alpha_ctrl , flag=wx.EXPAND, border=5)
		
		# gamma parameter ctr
		sizer.Add (wx.StaticText(self, label="gamma = "), flag=wx.LEFT, border=5)
		gamma_ctrl = wxFloatSpin (self, min_val=0, max_val=1, increment=0.01, value=0.1, digits=3)
		gamma_ctrl.__label__ = "gamma"
		sizer.Add (gamma_ctrl , flag=wx.EXPAND, border=5)
		
		# Min value of delta
		sizer.Add (wx.StaticText(self, label="\nmin delta (2*Pi)"), flag=wx.LEFT, border=5)
		delta_min_ctrl = wxFloatSpin (self, min_val=0, max_val=10, increment=0.01, value=0., digits=3)
		delta_min_ctrl.__label__ = "delta_min"
		sizer.Add (delta_min_ctrl , flag=wx.EXPAND, border=5)
		
		# Max value of delta
		sizer.Add (wx.StaticText(self, label="max delta (2*Pi)"), flag=wx.LEFT, border=5)
		delta_max_ctrl = wxFloatSpin (self, min_val=0, max_val=10, increment=0.01, value=2., digits=3)
		delta_max_ctrl.__label__ = "delta_max"
		sizer.Add (delta_max_ctrl , flag=wx.EXPAND, border=5)
		
		# Max value of delta
		sizer.Add (wx.StaticText(self, label="number of delta scans"), flag=wx.LEFT, border=5)
		delta_num_ctrl = wx.SpinCtrl (self, value="10", min=1, max=100000)
		delta_num_ctrl.__label__ = "delta_num"
		sizer.Add (delta_num_ctrl , flag=wx.EXPAND, border=5)
		
		# Separator
		sizer.Add (wx.StaticText(self), border=5)
		
		#######################################################
		self.get_delta_scan = wx.Button (self, label="Get delta scans")
		self.get_delta_scan.Bind( wx.EVT_BUTTON, self.GetDeltaScan )
		sizer.Add (self.get_delta_scan , flag=wx.EXPAND, border=5)
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		self.CreateSettingsDict()
		
		
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
		
		# Check whether the background signal array is present
		#self.CheckBackground()
		
		miips_settings = self.GetSettings()
		
		#####################################################################
		
		# Get range of delta
		delta_range = np.linspace( miips_settings["delta_min"], miips_settings["delta_max"], miips_settings["delta_num"] )
		delta_range *= (2*np.pi)
		
		gamma_w = miips_settings["gamma"] * np.arange(self.num_pixels)
		ones_ampl = np.ones(self.num_pixels)
		
		# Get spectral scans
		delta_scans = []
		
		for delta in delta_range :
			
			# Make sure software not frosen
			wx.Yield()
		
			# reference phase: f(w) = alpha*cos( gamma * w + delta )
			reference_phase = np.cos( gamma_w + delta )
			# make it in the rage [0,1]
			reference_phase -= reference_phase.min()
			reference_phase /= reference_phase.max()
		
			# scale by the value of the magnitude
			reference_phase *= miips_settings["alpha"]
			
			# Set the pulse shape
			self.DevPulseShaper.SetAmplPhase(ones_ampl, reference_phase) 
			
			# Get spectrum
			spectrum = self.DevSpectrometer.AcquiredData()
			
			# Vertical binning
			spectrum = ( spectrum.sum(axis=0) if len(spectrum.shape) == 2 else spectrum )
		
			# Save spectrum
			delta_scans.append( spectrum )
		
		delta_scans = np.array(delta_scans)
		
		################ Plotting results ####################
		visvis.cla(); visvis.clf(); 
		
		visvis.imshow( delta_scans, cm=visvis.CM_JET )
		visvis.ylabel ('delta'); visvis.xlabel ('wavelegth')