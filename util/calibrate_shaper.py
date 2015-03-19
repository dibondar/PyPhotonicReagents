"""
Calibrate shaper
"""
# Add main directory to enable imports  
if __name__ == '__main__' :
	import os
	os.sys.path.append(os.path.abspath('..'))

import wx, h5py, time
import numpy as np
from scipy.interpolate import pchip_interpolate, PchipInterpolator
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import argrelmin

# Real time plotting
import visvis

# GUI components
from libs.gui.basic_window import BasicWindow, SaveSettings
from libs.gui.hardware_control import HardwareGUIControl

# Hardware
from libs.dev.spectrometer_ocean_optics import ManagerOceanOpticsSpectrometer as ManagerSpectrometer
from libs.dev.spectrometer_ocean_optics import OceanOpticsSpectrometerTab as SpectrometerTab
#from libs.dev.camera_istar import ManagerIStarCamera as ManagerSpectrometer
#from libs.dev.camera_istar import IStarCameraTab as SpectrometerTab

from libs.dev.pulse_shaper import *

########################################################################
#
#	The following are calibration curves provided by the manufacturer.
#	They will be used as initial guesses.
#
########################################################################


voltages_trial = np.array([600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950,
975, 1000, 1050, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
2250, 2500, 2750, 3000, 3250, 3500, 3750, 4095])

# ~633nm (provided by CRi)
vis_modulation = np.array([1882, 1842.8, 1789.9, 1730.6, 1660.9, 1595.2, 1543.6, 1484.2, 1414.1, 1357.2, 1300.1,
1236.5, 1184, 1138.2, 1090, 1042.7, 991.37, 928.2, 850.1, 728.15, 633, 563.15, 507.28,
461.38, 422.72, 389.57, 360.27, 337.86, 297, 265.25, 236.99, 214.65, 196.24, 180.61, 167.26, 151.62])

# Converting wavelength modulation into phase modulation
vis_modulation *= (2*np.pi / 633.)

# ~800nm (provided by CRi)
nir_modulation = np.array([2907.5, 2858.1, 2742.7, 2651.4, 2545.5, 2456, 2357.1, 2261.1, 2168.6, 2074.5,
1989.1, 1899, 1818.2, 1738.4, 1665.2, 1582.5, 1530.8, 1408.2, 1301.5, 1116.1, 971.62, 868.72, 781.03, 711.87,
648.03, 602.2, 564.07, 526.33, 454.12, 401.6, 358.47, 322.02, 299.93, 274.91, 254.31, 230.19])

# Converting wavelength modulation into phase modulation
nir_modulation *= (2*np.pi / 800. )

########################################################################

class CalibrateShaperTab (HardwareGUIControl) :
	"""
	Settings for shaper calibration
	"""
	def __init__(self, parent) :
		HardwareGUIControl.__init__(self, parent)

		sizer = wx.BoxSizer(wx.VERTICAL)
		
		# Starting pixel to calibrate 
		sizer.Add (wx.StaticText(self, label="Initial pixel"), flag=wx.LEFT, border=5)
		initial_pixel_ctr = wx.SpinCtrl (self, value="0", min=0, max=640)
		initial_pixel_ctr.SetLabel("initial pixel")
		sizer.Add (initial_pixel_ctr, flag=wx.EXPAND, border=5)
		
		# Final pixel to calibrate
		sizer.Add (wx.StaticText(self, label="Final pixel"), flag=wx.LEFT, border=5)
		final_pixel_ctr = wx.SpinCtrl (self, value="640", min=0, max=640)
		final_pixel_ctr.SetLabel("final pixel")
		sizer.Add (final_pixel_ctr, flag=wx.EXPAND, border=5)
		
		# Size of pixel bundle
		sizer.Add (wx.StaticText(self, label="\nPixels to bundle"), flag=wx.EXPAND, border=5)
		self.pixel_bundle_width = wx.SpinCtrl(self, value="1", min=1, max=640)
		self.pixel_bundle_width.SetLabel ("pixel bundle width")
		sizer.Add (self.pixel_bundle_width, flag=wx.EXPAND, border=5)
		
		# Correspondence between wavelength and pulse shaper pixel number 
		# entered as dictionary without the "{}"
		sizer.Add (wx.StaticText(self, 
			label="\nPixel to wavelength entered as \n<pixel1> : <lambda1>,\n<pixel2> : <lambda2>, etc.")
		, flag=wx.EXPAND, border=5)
		pixel2lambda_ctrl = wx.TextCtrl (self, value="", style=wx.TE_MULTILINE|wx.EXPAND)
		pixel2lambda_ctrl.__label__ = "pixel_to_lamba"
		sizer.Add (pixel2lambda_ctrl, flag=wx.EXPAND, border=5)
		
		# Initial voltage
		sizer.Add (wx.StaticText(self, label="\nInitial voltage"), flag=wx.LEFT, border=5)
		initial_voltage_ctr = wx.SpinCtrl (self, value="500", min=0, max=PULSESHAPER_MAX_VAL)
		initial_voltage_ctr.SetLabel("initial voltage")
		sizer.Add (initial_voltage_ctr, flag=wx.EXPAND, border=5)
		
		# Final voltage
		sizer.Add (wx.StaticText(self, label="Final voltage"), flag=wx.LEFT, border=5)
		final_voltage_ctr = wx.SpinCtrl (self, value="800", min=0, max=PULSESHAPER_MAX_VAL)
		final_voltage_ctr.SetLabel("final voltage")
		sizer.Add (final_voltage_ctr, flag=wx.EXPAND, border=5)
		
		# Voltage step size
		sizer.Add (wx.StaticText(self, label="Voltage scanning step"), flag=wx.LEFT, border=5)
		voltage_step_ctr = wx.SpinCtrl (self, value="10", min=0, max=PULSESHAPER_MAX_VAL)
		voltage_step_ctr.SetLabel("voltage step")
		sizer.Add (voltage_step_ctr, flag=wx.EXPAND, border=5)
	
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		self.CreateSettingsDict()
		
########################################################################

class SettingsNotebook (wx.Notebook) :
	"""
	GUI for listing all settings
	"""
	def __init__(self, parent, DevSpectrometer, DevPulseShaper):
		wx.Notebook.__init__(self, parent)
		
		self.CalibrateShaper = CalibrateShaperTab(self)
		self.AddPage(self.CalibrateShaper, "Calibrate shaper")
		
		self.Spectrometer = SpectrometerTab(self, DevSpectrometer)
		self.AddPage (self.Spectrometer, "OO Spectrometer settings")
		 
		self.PulseShaper = PulseShaperTab(self, DevPulseShaper)
		self.AddPage (self.PulseShaper, "Pulse shaper settings")

		# Dictionary to bind names to tabs for saving and loading settings
		self.settings_to_tabs = {"Spectrometer" : self.Spectrometer, 
			"PulseShaper" : self.PulseShaper, "CalibrateShaper" : self.CalibrateShaper }
			
########################################################################

class CalibrateShaper (BasicWindow) :

	def __init__ (self, parent) :
		# Starting spectrometer
		self.Spectrometer = ManagerSpectrometer()
		self.SpectrometerProc = self.Spectrometer.start()
		
		# Starting pulse shaper
		self.PulseShaper = ManagerShaper()
		self.PulseShaperProc = self.PulseShaper.start()
		
		# Create GUI
		dw, dh = wx.DisplaySize()
		wx.Frame.__init__ (self, parent, title="Pulse shaper calibration with Ocean Optics Spectrometer",
								size=(0.9*dw, 0.88*dh) )
		
		self.ConstructGUI ()
		self.Center()
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
		
		# Vary one pixel button
		boxsizer.Add (wx.StaticText(self.panel, label="\nPixel to vary"), flag=wx.EXPAND, border=5)
		self.pixel_to_vary = wx.SpinCtrl(self.panel, value="320", min=0, max=640)
		boxsizer.Add (self.pixel_to_vary, flag=wx.EXPAND, border=5)
		
		self.vary_pixel_bundle_button = wx.Button (self.panel)
		self.vary_pixel_bundle_button.__start_label__ = "Vary pixel bundle"
		self.vary_pixel_bundle_button.__stop_label__ = "STOP varying"
		self.vary_pixel_bundle_button.SetLabel (self.vary_pixel_bundle_button.__start_label__)
		self.Bind (wx.EVT_BUTTON, self.VaryPixelBundle, self.vary_pixel_bundle_button)
		boxsizer.Add (self.vary_pixel_bundle_button, flag=wx.EXPAND, border=5)
		
		# Separator
		boxsizer.Add (wx.StaticText(self.panel), flag=wx.EXPAND, border=5)
		
		################## Calibrate button ##################
		self.calibrate_button = wx.Button (self.panel)
		self.calibrate_button.Bind (wx.EVT_LEFT_DOWN, self.PerformCalibration)
		self.calibrate_button.Bind (wx.EVT_LEFT_DCLICK, self.PerformCalibration)
		boxsizer.Add(self.calibrate_button, flag=wx.EXPAND, border=5)
		# Define labels
		self.calibrate_button.__start_label__ 	= "Calibrate pulse shaper"
		self.calibrate_button.__pause_label__ 	= "PAUSE calibration"
		self.calibrate_button.__resume_label__	= "RESUME calibration"
		self.calibrate_button.__stop_label__ 	= "STOP calibration"
		self.calibrate_button.SetLabel (self.calibrate_button.__start_label__)
		
		# Extract phases functions
		self.extract_phase_function = wx.Button (self.panel, label="Extract calibration phase functions")
		self.Bind (wx.EVT_BUTTON, self.ExtractPhaseFunc, self.extract_phase_function)
		boxsizer.Add (self.extract_phase_function, flag=wx.EXPAND, border=5)
		
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
		
	def ScanVoltage (self) :
		"""
		Using the iterator <self.scan_pixel_voltage_pair> record the spectral response 
		by applying the voltages
		"""
		# Pause calibration, if user requested
		try : 
			if self.pause_calibration : return
		except AttributeError : return
				
		try :
			param = self.scan_pixel_voltage_pair.next()
			self.PulseShaper.SetUniformMasks(*param)
			
			# Getting spectrum
			spectrum = self.Spectrometer.AcquiredData() 
			# Save the spectrum
			try : self.SpectraGroup["voltages_%d_%d" % param] = spectrum
			except RuntimeError : print "There was RuntimeError while saving scan voltages_%d_%d" % param
			
			# Plot the spectra
			visvis.gca().Clear()
			
			visvis.plot (self.wavelengths, spectrum)
			visvis.xlabel("wavelength (nm)")
			visvis.ylabel("counts")
			
			# Scanning progress info
			self.scanned += 1.
			percentage_completed = 100.*self.scanned/self.scan_length
			seconds_left = ( time.clock() - self.initial_time )*(100./percentage_completed - 1.)
			# convert to hours:min:sec
			m, s = divmod(seconds_left, 60)
			h, m = divmod(m, 60)
			
			title_info = param + (percentage_completed, h, m, s)
			visvis.title ("Scanning spectrum by applying voltages %d/%d. Progress: %.1f%% completed. Time left: %d:%02d:%02d." %  title_info)
			
			self.fig.DrawNow()
			
			# Measure the next pair
			wx.CallAfter(self.ScanVoltage)
		except StopIteration :
			# Perform processing of the measured data
			wx.CallAfter(self.ExtractPhaseFunc, filename=self.calibration_file.filename)
			# All voltages are scanned
			self.StopAllJobs()
			# Sop using the shaper
			self.PulseShaper.StopDevice()
	
	def FitSpectralScans (self, scans, voltages, pixels_edges) :
		"""
		Perform fitting to the pulse shaper's mask transmission coefficient 
		"""
		
		def FitIndividualPixel (voltages, pixel, modulation, p0=None) :
			"""
			Find fit for individual pixel with initial guess for phase function given by `modulation`.
			`voltages` voltage values for which `pixel` was measured
			`pixel` measured transmission (to be fitted) 
			`p0` is the initial guess for fitting parametres
			"""	
			def GetVM (V0, V1, m0, m1) :
				"""
				Return voltage and phase modulation from parameters
				"""
				M = m0 + m1*modulation
				V = np.linspace(V0, V1, M.size)
				return V, M
			def FittedTransmission  (voltages, offset, amplitude, *params) :
				"""
				Return the transmission function for a shaper
				"""
				V, M = GetVM(*params)
				return amplitude*np.cos( pchip_interpolate(V,M,voltages) )**2 + offset
			
			# Set fitting parameters to their default values
			if p0 is None : p0 = [0., 1., voltages.min(), voltages.max(), 0., 1.]
			
			# Fitting the transmission
			try : popt, _ = curve_fit(FittedTransmission, voltages, pixel, p0=p0)
			except RuntimeError : popt = p0
				
			# Get fitting error
			fitting_error = np.sum( ( FittedTransmission(voltages, *popt) - pixel )**2 )
						
			return fitting_error, GetVM(*popt[2:]), popt
		
		############################################################################
		# Selecting the voltage range
		V_min = max( voltages_trial.min(), voltages.min() )
		V_max = min( voltages_trial.max(), voltages.max() )
		indx = np.nonzero( (V_min <= voltages)&(voltages <= V_max) )
		
		# Number of calibration points lying within the voltage region
		num_vol_trial = np.sum( (V_min <= voltages_trial)&(voltages_trial <= V_max) ) 
		if num_vol_trial < 2 : num_vol_trial = 2
		
		# Re-sample modulation provided by CRi so that the voltage is equidistantly spaced
		resampled_vis_modulation = pchip_interpolate(voltages_trial, vis_modulation,
			np.linspace(V_min, V_max, min(len(voltages), num_vol_trial) )
		)
		resampled_nir_modulation = pchip_interpolate(voltages_trial, nir_modulation,
			np.linspace(V_min, V_max, min(len(voltages), num_vol_trial) )
		)
		
		# Normalizing scans
		scans -= scans.min(axis=0); scans /= scans.max(axis=0)
		
		# Bin the spectrum into pixels
		spectral_slices = map( lambda begin,end : scans[:,begin:end].mean(axis=1), pixels_edges[:-1], pixels_edges[1:] ) 
		
		# List containing calibration data for each pixel
		calibration = []
		
		# Initial guesses for fitting
		vis_p0 = None; nir_p0 = None;
		
		# Fit individual pulse shaper pixels them
		for pixel_num, pixel in enumerate(spectral_slices) :
		
			# Smoothing and normalizing each pixel 
			#pixel = gaussian_filter(pixel,sigma=1)
			pixel -= pixel.min(); pixel /= pixel.max()
			
			# Fit the pixel by using the vis calibration curve as the initial guess
			vis_err, vis_calibration, vis_p0 = FitIndividualPixel(voltages[indx], pixel[indx], 
																	resampled_vis_modulation, vis_p0)
			
			# Fit the pixel by using the nir calibration curve as the initial guess
			nir_err, nir_calibration, nir_p0 = FitIndividualPixel(voltages[indx], pixel[indx],
																	resampled_nir_modulation, nir_p0)
			
			# Choose the best fit
			if nir_err > vis_err : 
				calibation_voltage, calibration_phase = vis_calibration
				fit_err = vis_err
			else : 
				calibation_voltage, calibration_phase = nir_calibration
				fit_err = nir_err
				
			###################### Plot ######################## 
			visvis.clf()
			
			# Plot measured data
			visvis.plot( voltages, pixel, lc='r',ms='*', mc='r')
		 
			# Plot fitted data
			plot_voltages = np.linspace( calibation_voltage.min(), calibation_voltage.max(), 500) 
			transmission_fit = np.cos( pchip_interpolate(calibation_voltage, calibration_phase, plot_voltages) )**2
			visvis.plot( plot_voltages, transmission_fit, lc='b')
			
			visvis.title ('Calibrating pixel %d / %d' % (pixel_num, len(spectral_slices)-1) )
			visvis.legend(['measured', 'fitted'])
			visvis.xlabel ('voltages')
			visvis.ylabel ('Transmission coefficient')
			 
			self.fig.DrawNow()
			############ Save the calibration data ##################
			calibration.append( ( calibation_voltage, calibration_phase, fit_err ) )
		
		return calibration
	
	def GetDispersionPhaseCurves (self, wavelengths, scans, voltages, pixels_edges, 
		calibration, pixel2lambda_func, pulse_shaper_pixel_num ) :
		"""
		Surface calibration:
		Find dispersion and phase curves based on a result of the method `self.FitSpectralScans`
			`calibration` is the return of `self.FitSpectralScans`
		"""
		import operator		
		# Find the calibration function that was best fitted
		best_calibation_voltage, best_calib_phase, _ = min(calibration, key=operator.itemgetter(2))
		
		# Form a function representing calibration curve (best fitted)
		# Normalizing calibration curve
		best_calib_phase -= best_calib_phase.min()
		best_calib_phase /= best_calib_phase.max()
		best_calibration_curve = PchipInterpolator (best_calibation_voltage, best_calib_phase)
		
		# Selecting the voltage range
		V_min = max( best_calibation_voltage.min(), voltages.min() )
		V_max = min( best_calibation_voltage.max(), voltages.max() )
		V_indx = np.nonzero( (V_min <= voltages)&(voltages <= V_max) )
		
		# Select scanning range that corresponds to a valid region of voltages and wavelength
		scans = scans[V_indx[0], pixels_edges[0]:pixels_edges[-1]]
		
		# Wavelength position of pixels in a new sliced scan
		wavelengths_cut = wavelengths[ pixels_edges[0]:pixels_edges[-1] ]
		
		# Central wavelength of each pixels
		logical_pixel_lambda = 0.5*( wavelengths[pixels_edges[1:]] + wavelengths[pixels_edges[:-1]] )
		
		# Construct the initial guess for the dispersion curve
		# use the cubic polynomial interpolation
		# normalize it to the best fitted calibration curve
		best_min 	= best_calibration_curve(V_min)
		best_max	= best_calibration_curve(V_max)
		
		# values of the calibration curves
		calibration_values = best_calibration_curve( voltages[V_indx] )[:,np.newaxis]
		
		def TransmissionC (calibration_values, params) :
			"""
			Transmission coefficient for fitting
			"""
			offset 		= params[:len(params)/2]
			multiplier 	= params[len(params)/2:]
			phase = calibration_values * np.polyval(multiplier, wavelengths_cut)
			phase += np.polyval(offset, wavelengths_cut) 
			return np.cos( phase )**2
			
		def Fit_TransmissionC (calibration_values, *params) :
			return np.ravel( TransmissionC(calibration_values, params) )
			
		# Initial guess for fitting parameters
		c_min, c_max = zip(*[ pchip_interpolate(c[0],c[1],[V_min, V_max]) for c in calibration ])
		c_min = np.array(c_min); c_max = np.array(c_max)
		
		# Use different power fits 
		power_fits = []
		for power in [1, 3, 5, 7] :
			offset 		= np.polyfit(logical_pixel_lambda, (best_max*c_min-best_min*c_max)/(best_max-best_min), power)
			multiplier	= np.polyfit(logical_pixel_lambda, (c_max-c_min)/(best_max-best_min), power)
			p0=np.append(offset, multiplier)
			
			try :
				popt, _ = curve_fit(Fit_TransmissionC, calibration_values, np.ravel(scans), p0=p0)
			except RuntimeError : popt = p0
			
			# Calculate the Transmission coefficients for plotting 
			TC_fitted = TransmissionC(calibration_values, popt)
		
			# Calculate fitting error
			error = np.sum( (TC_fitted - scans)**2 )
			power_fits.append( (error, popt, TC_fitted) )
		
		# Select the best power fit
		_, popt, TC_fitted = min(power_fits)
		
		# Extracted the fitted parameters
		offset		= popt[:len(popt)/2]
		multiplier	= popt[len(popt)/2:] 
		
		# Get wavelength for each physical pixel in the pulse shaper
		physical_pixel_lambda = pixel2lambda_func(np.arange(pulse_shaper_pixel_num))
		
		# Calculate offset and multiplier for each physical pixel
		offset 		= np.polyval(offset, physical_pixel_lambda)
		multiplier	= np.polyval(multiplier, physical_pixel_lambda)
		
		# Return results 
		return TC_fitted, scans, {	"offset" 					: offset,
									"multiplier"				: multiplier,
									"calibration_curve_voltage" : best_calibation_voltage, 
									"calibration_curve_phase" 	: best_calib_phase }
			
	def ExtractPhaseFunc (self, event=None, filename=None) :
		"""
		This function is the last stage of calibration, 
		when measured data is mathematically processed to obtain the calibration curves.
		"""
		
		# If filename is not specified, open the file dialogue
		if filename is None : 
			filename = self.LoadSettings (title="Load calibration file...")
			# Update the name of pulse shaper calibration file
			import os
			self.SettingsNotebook.PulseShaper.SetSettings({"calibration_file_name" : os.path.abspath(filename)})
			
		visvis.clf()
		
		# Loading the file calibration file
		with h5py.File(filename, 'a') as calibration_file :
			############### Loading data ####################
			wavelengths 		= calibration_file["calibration_settings/wavelengths"][...]
			fixed_voltage 		= calibration_file["calibration_settings/fixed_voltage"][...]
			pulse_shaper_pixel_num = calibration_file["calibration_settings/pulse_shaper_pixel_num"][...]
			initial_pixel	 	= calibration_file["settings/CalibrateShaper/initial_pixel"][...]
			final_pixel 		= calibration_file["settings/CalibrateShaper/final_pixel"][...]
			pixel_bundle_width 	= calibration_file["settings/CalibrateShaper/pixel_bundle_width"][...]
			pixel_to_lamba		= str(calibration_file["settings/CalibrateShaper/pixel_to_lamba"][...])
			# Convert to dict
			pixel_to_lamba 		= eval( "{%s}" % pixel_to_lamba )
			
			# Loading scans of slave and masker masks
			master_mask_scans = []; slave_mask_scans = []
			for key, spectrum in calibration_file["spectra_from_uniform_masks"].items() :
				master_volt, slave_volt = map(int, key.split('_')[-2:])
				if master_volt == fixed_voltage : slave_mask_scans.append( (slave_volt, spectrum[...]) )
				if slave_volt == fixed_voltage	: master_mask_scans.append( (master_volt, spectrum[...]) )
			
			# Sort by voltage
			master_mask_scans.sort(); slave_mask_scans.sort()
			
			# Extract spectral scans and voltages for each mask
			master_mask_voltage, master_mask_scans = zip(*master_mask_scans) 
			master_mask_voltage = np.array(master_mask_voltage); master_mask_scans = np.array(master_mask_scans)
			
			slave_mask_voltage, slave_mask_scans = zip(*slave_mask_scans)
			slave_mask_voltage = np.array(slave_mask_voltage); slave_mask_scans = np.array(slave_mask_scans)
			
			################### Find the edges of pixels #################
			
			# function that converts pixel number to pulse shaper
			# `pixel_to_lamba` is a dictionary  with key = pixel, value = lambda
			deg = min(2,len(pixel_to_lamba)-1)
			pixel2lambda_func = np.poly1d( np.polyfit(pixel_to_lamba.keys(), pixel_to_lamba.values(), deg=deg ) )
			#lambda2pixel_func = np.poly1d( np.polyfit(pixel_to_lamba.values(), pixel_to_lamba.keys(), deg=deg ) )
			
			# Get pixels_edges in shaper pixels 
			pixels_edges_num = np.arange(initial_pixel, final_pixel, pixel_bundle_width)
			if pixels_edges_num[-1] < final_pixel :
				pixels_edges_num = np.append( pixels_edges_num, [final_pixel])
			
			# Get pixels_edges in logical_pixel_lambda
			pixels_edges = pixel2lambda_func(pixels_edges_num)
			
			# Get pixels_edges in positions of the spectrometer spectra
			pixels_edges = np.abs(wavelengths - pixels_edges[:,np.newaxis]).argmin(axis=1)
			
			# Sorting
			indx = np.argsort(pixels_edges)
			pixels_edges = pixels_edges[indx]
			pixels_edges_num = pixels_edges_num[indx]
			
			# Plot
			visvis.cla(); visvis.clf()
			visvis.plot( pixel_to_lamba.values(), pixel_to_lamba.keys(), ls=None, ms='*', mc='g', mw=15)
			visvis.plot ( wavelengths[pixels_edges], pixels_edges_num, ls='-',lc='r')
			visvis.xlabel('wavelength (nm)')
			visvis.ylabel ('pulse shaper pixel')
			visvis.legend( ['measured', 'interpolated'])
			self.fig.DrawNow()
			 
			################ Perform fitting of the phase masks #################	
			master_mask_fits = self.FitSpectralScans(  master_mask_scans, master_mask_voltage, pixels_edges )
			slave_mask_fits = self.FitSpectralScans(  slave_mask_scans, slave_mask_voltage, pixels_edges )
			
			################# Perform fitting of dispersion and phase functions #################
			master_mask_TC_fitted, master_mask_scans, master_mask_disp_ph_curves = \
				self.GetDispersionPhaseCurves (wavelengths, master_mask_scans, master_mask_voltage, 
					pixels_edges, master_mask_fits, pixel2lambda_func, pulse_shaper_pixel_num ) 
			
			slave_mask_TC_fitted, slave_mask_scans, slave_mask_disp_ph_curves = \
				self.GetDispersionPhaseCurves (wavelengths, slave_mask_scans, slave_mask_voltage, 
					pixels_edges, slave_mask_fits, pixel2lambda_func, pulse_shaper_pixel_num ) 
			
			################ Saving fitting parameters ####################
			
			#################################################################
			# Save surface calibration
			try : del calibration_file["calibrated_surface"]
			except KeyError : pass
			
			CalibratedSurfaceGroupe = calibration_file.create_group("calibrated_surface")
			
			master_mask_calibrated_surface = CalibratedSurfaceGroupe.create_group("master_mask")
			for key, value in master_mask_disp_ph_curves.items() :
				master_mask_calibrated_surface[key] = value
			
			slave_mask_calibrated_surface = CalibratedSurfaceGroupe.create_group("slave_mask")
			for key, value in slave_mask_disp_ph_curves.items() :
				slave_mask_calibrated_surface[key] = value
			
			#################################################################
			# Clear the group, if it exits 
			try : del calibration_file["calibrated_pixels"]
			except KeyError : pass
			
			# This group is self consistent, thus the redundancy in saved data (with respect to other group in the file)
			CalibratedPixelsGroupe = calibration_file.create_group ("calibrated_pixels")
			CalibratedPixelsGroupe["pixels_edges"] 			= pixels_edges
			
			# Finding spectral bounds for each pixel
			pixels_spectral_bounds = zip( wavelengths[pixels_edges[:-1]], wavelengths[pixels_edges[1:]] )
			
			# Pulse shaper pixel bounds
			pixels_bounds = zip(pixels_edges_num[:-1], pixels_edges_num[1:])
			
			PixelsGroup = CalibratedPixelsGroupe.create_group("pixels")
			for pixel_num, master_mask_calibration, slave_mask_calibration, \
				spectral_bound, pixel_bound in zip( range(len(master_mask_fits)), \
						master_mask_fits, slave_mask_fits, pixels_spectral_bounds, pixels_bounds  ) :
				pixel =	PixelsGroup.create_group("pixel_%d" % pixel_num)
				pixel["spectral_bound"] 		= spectral_bound
				pixel["pixel_bound"]			= pixel_bound
				# Saving fitted calibration data in the tabular form
				pixel["voltage_master_mask"] 	= master_mask_calibration[0]
				pixel["phase_master_mask"] 		= master_mask_calibration[1]
				pixel["voltage_slave_mask"] 	= slave_mask_calibration[0]
				pixel["phase_slave_mask"]		= slave_mask_calibration[1]
						
						
		################ Plotting results of calibration ####################
		visvis.cla(); visvis.clf(); 
		
		visvis.subplot(2,2,1)
		visvis.imshow( master_mask_scans, cm=visvis.CM_JET )
		visvis.title ("Master mask (measured data)")
		visvis.ylabel ('voltage'); visvis.xlabel ('wavelength (nm)')
		
		visvis.subplot(2,2,3)
		visvis.imshow( master_mask_TC_fitted, cm=visvis.CM_JET )
		visvis.title ("Master mask (fitted data)")
		visvis.ylabel ('voltage'); visvis.xlabel ('wavelength (nm)')
		
		visvis.subplot(2,2,2)
		visvis.imshow( slave_mask_scans, cm=visvis.CM_JET )
		visvis.title ("Slave mask (measured data)")
		visvis.ylabel ('voltage'); visvis.xlabel ('wavelength (nm)')
		
		visvis.subplot(2,2,4)
		visvis.imshow( slave_mask_TC_fitted, cm=visvis.CM_JET )
		visvis.title ("Slave mask (fitted data)")
		visvis.ylabel ('voltage'); visvis.xlabel ('wavelength (nm)')
		
		
	def PerformCalibration (self, event=None) :
		"""
		<self.calibrate_button> was clicked.
		The calibration of pulse shaper for 
		"""
		# Creating the reference
		button = self.calibrate_button 
	
		try :
			# Mouse double clicking stops scanning
			if event.GetEventType() == wx.wxEVT_LEFT_DCLICK  : button.SetLabel (button.__stop_label__)
		except AttributeError : pass
			
		if button.GetLabel() == button.__start_label__ :
			self.StopAllJobs ()
			
			filename = self.SaveSettings(title= button.__start_label__, filename="pulse_shaper_calibration.hdf5")
			if filename is None : return
			
			# Update the name of pulse shaper calibration file
			import os
			self.SettingsNotebook.PulseShaper.SetSettings({"calibration_file_name" : os.path.abspath(filename)})
			
			# get spectrometer's settings
			settings = self.SettingsNotebook.Spectrometer.GetSettings()
			# Initiate spectrometer
			if self.Spectrometer.SetSettings(settings) == RETURN_FAIL : return
			# Initiate pulse shaper
			settings = self.SettingsNotebook.PulseShaper.GetSettings()
			if self.PulseShaper.Initialize(settings) == RETURN_FAIL : return
			
			# Open the file for saving calibration data
			self.calibration_file = h5py.File (filename, 'a')

			# The HDF5 group where spectral scans will be saved
			try : del self.calibration_file["spectra_from_uniform_masks"]
			except KeyError : pass
			self.SpectraGroup = self.calibration_file.create_group ("spectra_from_uniform_masks")
			
			# The HDF5 group where the calibration parameters are stored
			try : del self.calibration_file["calibration_settings"]
			except KeyError : pass
			self.CalibrationGroup = self.calibration_file.create_group ("calibration_settings")
			
			# Get wavelengths
			self.wavelengths = self.Spectrometer.GetWavelengths()
			self.CalibrationGroup["wavelengths"] = self.wavelengths
			
			# Save zero mask
			self.fixed_mask = self.PulseShaper.GetZeroMask()
			
			# Save the pulse shaper resolution (physical pixels)
			self.CalibrationGroup["pulse_shaper_pixel_num"] = self.PulseShaper.GetPixelNumber()
			
			# Define the iterator for scanning the voltage difference
			
			# The general idea is that by measuring spectra we measure the transmission
			# coefficients T(V1, V2) = cos^2 ( phi_m(V1) - phi_s(V2) )
			# However, we do not need to measure T for all V1 and V2 if we use the following
			# strategy:
			# [ phi_m(V1) - phi_s(V2) ] = [ phi_m(V1) - phi_s(V0) ] + [ phi_m(V0) - phi_s(V2) ]
			#								- [ phi_m(V0) - phi_s(V0) ] 
			# where Vo is an an arbitrary fixed voltage
			
			# Set the iterator
			settings = self.SettingsNotebook.CalibrateShaper.GetSettings()
			initial_voltage = settings["initial_voltage"]
			final_voltage = settings["final_voltage"]
			voltage_step 	= settings["voltage_step"]
			
			self.fixed_voltage = int( 0.5*(final_voltage + initial_voltage) )  # This is V0
			self.CalibrationGroup["fixed_voltage"] = self.fixed_voltage 
			
			if self.fixed_voltage > final_voltage or self.fixed_voltage < initial_voltage : 
				voltage = np.arange(initial_voltage, final_voltage+1, voltage_step) 
			else :
				# Skip the fixed voltage
				voltage = np.append( np.arange(initial_voltage, self.fixed_voltage, voltage_step), 
										np.arange(self.fixed_voltage+1, final_voltage+1, voltage_step) )
										
			from itertools import product, chain
			self.scan_pixel_voltage_pair = chain( 	product(voltage, [self.fixed_voltage] ),
													product([self.fixed_voltage], voltage ) )
													#[(self.fixed_voltage, self.fixed_voltage)]   ) 
													
			# Variables for progress info
			self.scan_length = 2*len(voltage) + 1
			self.scanned = 0.
			self.initial_time = time.clock()
			
			# Changing the button's label 
			button.SetLabel (button.__pause_label__)
			
			visvis.clf()
			
			# Start calibrating
			self.pause_calibration = False
			wx.CallAfter(self.ScanVoltage)
	
		elif button.GetLabel() == button.__pause_label__ :
			self.pause_calibration = True; button.SetLabel (button.__resume_label__)
		
		elif button.GetLabel() == button.__resume_label__ :
			self.pause_calibration = False
			wx.CallAfter(self.ScanVoltage)
			button.SetLabel (button.__pause_label__)

		elif button.GetLabel() == button.__stop_label__ :
			# Closing the file
			self.calibration_file.close()
			# Delete all the attributes associated with scanning 
			del self.pause_calibration, self.calibration_file, self.SpectraGroup, self.CalibrationGroup
			
			button.SetLabel (button.__start_label__)
			
		else : raise ValueError ("Unrecognised button's label")
			
	def ShowSpectra_by_VaryingPixelBundle (self) :
		"""
		This method is affiliated to the method <self.VaryPixelBundle>
		"""
		# Exit if the iterator is not defined
		try : self.pixel_bundel_value_iter
		except AttributeError : return
		
		try :
			voltage = self.pixel_bundel_value_iter.next()
			
			# Set the mask for pixel bundle
			width = self.SettingsNotebook.CalibrateShaper.pixel_bundle_width.GetValue() / 2
			start_pixel_bundle = self.pixel_to_vary.GetValue()
			mask = np.copy(self.fixed_mask)
			mask[ max(start_pixel_bundle-width, 0):min(mask.size, start_pixel_bundle+width) ] = voltage
			self.PulseShaper.SetMasks( mask, self.fixed_mask)
		
			# Getting spectrum
			spectrum = self.Spectrometer.AcquiredData()
		
			# Plot the spectra
			visvis.gca().Clear()
			visvis.plot (self.wavelengths, spectrum)
			visvis.xlabel("wavelength (nm)")
			visvis.ylabel("counts")
			visvis.title ("Voltage %d / %d " % (voltage, self.fixed_mask[0]) )
			self.fig.DrawNow()
			
			# Going to the next iteration
			wx.CallAfter (self.ShowSpectra_by_VaryingPixelBundle)
		except StopIteration :
			# Finish the job
			self.StopAllJobs()
	
	def VaryPixelBundle (self, event) :
		"""
		<self.vary_pixel_bundle_button> was clicked
		"""
		button = self.vary_pixel_bundle_button
		if button.GetLabel() == button.__start_label__ :
			self.StopAllJobs ()
			
			# get spectrometer's settings
			settings = self.SettingsNotebook.Spectrometer.GetSettings()
			# Initiate spectrometer
			if self.Spectrometer.SetSettings(settings) == RETURN_FAIL : return
			# Initiate pulse shaper
			settings = self.SettingsNotebook.PulseShaper.GetSettings()
			if self.PulseShaper.Initialize(settings) == RETURN_FAIL : return
		
			# Get wavelengths
			self.wavelengths = self.Spectrometer.GetWavelengths()

		
			# Set the iterator 
			settings = self.SettingsNotebook.CalibrateShaper.GetSettings()
			voltage_step 	= settings["voltage_step"]
			initial_voltage = settings["initial_voltage"]
			final_voltage	= settings["final_voltage"]
			self.pixel_bundel_value_iter = iter(xrange(initial_voltage, final_voltage+voltage_step, voltage_step))

			# Save fixed mask
			self.fixed_mask = self.PulseShaper.GetZeroMask()
			self.fixed_mask += PULSESHAPER_MAX_VAL
			
			# Start variation
			wx.CallAfter (self.ShowSpectra_by_VaryingPixelBundle)
			
			# Change button's label
			button.SetLabel (button.__stop_label__)
			
		elif button.GetLabel() == button.__stop_label__ :
			del self.pixel_bundel_value_iter
			button.SetLabel (button.__start_label__) 
			
		else : raise ValueError("Label is not recognized") 
		

		
#########################################################################
if __name__ == '__main__' :
	app = visvis.use('wx')
	app.Create()
	CalibrateShaper (None)
	app.Run()