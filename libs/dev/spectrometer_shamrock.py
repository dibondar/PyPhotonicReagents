########################################################################
#
#	This module contains classes for controlling and GUI representation of
#	spectrometer Shamrok SR303i connected to camera iDUS DU401A-BV
#
########################################################################

from libs.gui.hardware_control import HardwareGUIControl
from libs.dev.basic_device import BasicDevice

import ctypes, wx, os, multiprocessing
from multiprocessing.sharedctypes import RawArray
import numpy as np
from libs.dev.consts import * 


########################################################################
#
#	Manager class that communicates with the process where
#	the device resides
#
########################################################################

class ManagerShamrockSpectrometer :
	"""
	Class that manges the Shamrock spectrometer
	"""
	def __init__ (self) :
		# Create the lock for device 
		self.lock = multiprocessing.Lock()
		# Create a pipe for communication
		self.parent_connection, self.child_connection = multiprocessing.Pipe()
		# Create the buffer 
		self.spectrometer_buffer = ShamrockSpectrometer.AllocateBuffer()
		
	def __del__ (self) :
		self.parent_connection.close()
		self.child_connection.close()

	def start(self) :
		"""
		Start the process controlling the spectrometer
		"""
		p = ShamrockSpectrometer(self.child_connection, self.spectrometer_buffer)
		p.start()
		return p
	
	def run(self, command, arguments=None) :
		"""
		Send the command to the spectrometer through the pipe
		"""
		self.lock.acquire()
		self.parent_connection.send( (command, arguments) )
		result = self.parent_connection.recv()
		self.lock.release()
		return result
		
	def exit(self) :
		"""
		Close the process
		"""
		return self.run("Exit")

	def SetSettings (self, settings) :
		"""
		Set settings for the spectrometer and camera
		"""
		return self.run("SetSettings", settings)
	
	def AcquiredData (self) :
		"""
		Get the spectral data
		"""
		if self.run("AcquiredData") == RETURN_FAIL : print "Spectrometer Acquisition failed"
		return np.frombuffer(self.spectrometer_buffer, dtype=np.dtype(ctypes.c_long))

	def GetWavelengths (self) : 
		"""
		Get the calibration data.
		"""
		return self.run("GetWavelengths")
	
	def GetTemperature (self) : return self.run("GetTemperature")
	
########################################################################
#
#	Process where the device resides
#
########################################################################

########################################################################

# Constants from ShamrockCIF.H
SHAMROCK_SUCCESS = 20202 

# Constants from ATMCD32D.H
DRV_SUCCESS 		= 20002
DRV_VXDNOTINSTALLED = 20003
DRV_ERROR_FILELOAD 	= 20006
DRV_ERROR_PAGELOCK 	= 20010
DRV_ERROR_ACK 		= 20013
DRV_ACQ_BUFFER 		= 20018
DRV_KINETIC_TIME_NOT_MET = 20022
DRV_ACCUM_TIME_NOT_MET 	= 20023
DRV_NO_NEW_DATA 	= 20024
DRV_SPOOLERROR 		= 20026
DRV_TEMP_OFF 		= 20034
DRV_TEMP_NOT_STABILIZED = 20035
DRV_TEMP_STABILIZED = 20036
DRV_TEMP_NOT_REACHED = 20037
DRV_TEMP_DRIFT 		= 20040
DRV_FLEXERROR 		= 20053
DRV_P1INVALID 		= 20066
DRV_P2INVALID 		= 20067
DRV_INIERROR 		= 20070
DRV_COFERROR 		= 20071
DRV_ACQUIRING		= 20072
DRV_IDLE			= 20073
DRV_TEMPCYCLE 		= 20074
DRV_NOT_INITIALIZED = 20075
DRV_USBERROR 		= 20089
DRV_ERROR_NOCAMERA 	= 20990
DRV_NOT_SUPPORTED 	= 20991

########################################################################

AndorDriverFolder = "C:/Program Files/Andor SOLIS/Drivers/"
ShamrockDriverFolder = "C:/Program Files/Andor SOLIS/Drivers/Shamrock/"

########################################################################

# Transcript of Error messages

CameraErrorMsg = { DRV_VXDNOTINSTALLED : "VxD not loaded",
	DRV_INIERROR : "Unable to load 'DETECTOR.INI'", DRV_COFERROR : "Unable to load '*.COF'",
	DRV_FLEXERROR : "Unable to load '*.RBF'", DRV_ERROR_ACK : "Unable to communicate with card",
	DRV_ERROR_FILELOAD : "Unable to load '*.COF' or '*.RBF' files", DRV_ERROR_PAGELOCK : "Unable to acquire lock on requested memory",
	DRV_USBERROR : "Unable to detect USB device or not USB2.0", DRV_ERROR_NOCAMERA : "No camera found",
	DRV_NOT_INITIALIZED : "System not initialized", DRV_ACQUIRING : "Acquisition in progress",
	DRV_P1INVALID : "Invalid readout mode passed", DRV_TEMPCYCLE : "Executing temperature cycle",
	DRV_ACCUM_TIME_NOT_MET : "Unable to meet Accumulate cycle time", DRV_KINETIC_TIME_NOT_MET : "Unable to meet Kinetic cycle time",
	DRV_ACQ_BUFFER : "Computer unable to read the data via the ISA slot", DRV_SPOOLERROR : "Overflow of the spool buffer",
	DRV_P2INVALID : "Array size is incorrect", DRV_NO_NEW_DATA : "No acquisition has taken place",
	DRV_TEMP_OFF : "Temperature is OFF", DRV_TEMP_STABILIZED : "Temperature has stabilized at set point",
	DRV_TEMP_NOT_REACHED : "Temperature has not reached set point", DRV_TEMP_DRIFT : "Temperature had stabilized but has since drifted",
	DRV_TEMP_NOT_STABILIZED : "Temperature reached but not stabilized", DRV_NOT_SUPPORTED : "Capability not supported"
}

class ShamrockSpectrometer (BasicDevice) :
	"""
	Control spectrometer and camera hardware from a separate process
	"""	
	# Resolution of the camera image (these values may need to be adjusted when porting the code)
	xpixels 	= 1600
	ypixels 	= 200
	
	def __init__ (self, pipe, spectrometer_buffer) :
		"""
		Initialize the spectrometer and camera
		"""
		BasicDevice.__init__(self, pipe)
		# saving the buffer where the spectrum will be saved
		self.buffer = spectrometer_buffer
		
	@classmethod
	def AllocateBuffer (cls) :
		"""
		This static method allocates buffer that corresponds to Full Vertical Binning readout mode
		"""
		return RawArray (ctypes.c_long, cls.xpixels )
		
	def InitializeShamrock (self) :
		"""
		Initialize Shamrock spectrometer
		"""
		# Expanding PATH environmental variable
		os.environ["PATH"] += os.pathsep + ShamrockDriverFolder
		
		# Loading the spectrometer driver
		self.ShamrockLib = ctypes.WinDLL ("ShamrockCIF.dll")
		
		# Initializing
		if self.ShamrockLib.ShamrockInitialize(ShamrockDriverFolder) != SHAMROCK_SUCCESS : raise RuntimeError ("Error in ShamrockInitialize!")
		
		"""
		# Verifying that there is a single Shamrock spectrometer
		totalSpectrometer = ctypes.c_int()
		if self.ShamrockLib.ShamrockGetNumberDevices( ctypes.byref(totalSpectrometer) ) != SHAMROCK_SUCCESS : raise RuntimeError ("Error in ShamrockGetNumberDevices!")
		if totalSpectrometer.value > 1: raise RuntimeError ("More than one Shamrock spectrometer!")
		"""
		
	def InitializeCamera (self) :
		"""
		Initialize iDus camera 
		"""
		# Expanding PATH environmental variable
		os.environ["PATH"] += os.pathsep + AndorDriverFolder
		
		# Loading the spectrometer driver
		self.CameraLib = ctypes.WinDLL ("atmcd32d.dll")
		
		# Verifying that there is a single Andor camera
		totalCameras = ctypes.c_int()
		if self.CameraLib.GetAvailableCameras( ctypes.byref(totalCameras) ) != DRV_SUCCESS : 
			raise RuntimeError ("Error in GetAvailableCameras")

		if totalCameras.value > 1 : 
			raise RuntimeError ("More than one Andor camera is present")
	
		# Initialize the camera	
		result = self.CameraLib.Initialize(AndorDriverFolder) 
		if result != DRV_SUCCESS : raise RuntimeError ("Error in Initialize: %s " % CameraErrorMsg[result])
		
		# Find out the number of pixels for figures
		__xpixels__ = ctypes.c_int(); __ypixels__ = ctypes.c_int()
		if self.CameraLib.GetDetector( ctypes.byref(__xpixels__), ctypes.byref(__ypixels__) ) != DRV_SUCCESS : raise RuntimeError ("Error in GetDetector")
		self.max_x_pixels = __xpixels__.value; self.max_y_pixels = __ypixels__.value
	
		# Check whether the static properties coincide with actual resolution
		if self.xpixels != self.max_x_pixels or self.ypixels != self.max_y_pixels :
			raise ValueError ("Static properties <xpixels> and <ypixels> of class <Spectrometer> have wrong values. Correct values are %d and %d. Source code must be modified." % ( self.max_x_pixels, self.max_y_pixels))
	
		return RETURN_SUCCESS
		
	def SetSettings (self, settings) :
		"""
		Assign settings 
		"""	
		self.SetCameraSettings(settings); self.SetShamrockSettings(settings)
		return RETURN_SUCCESS
		
	def SetCameraSettings (self, settings) :
		################ Camera settings ################
		# check the temperature range
		mintemp = ctypes.c_int(); maxtemp = ctypes.c_int()
		result = self.CameraLib.GetTemperatureRange( ctypes.byref(mintemp), ctypes.byref(maxtemp) )
		if result != DRV_SUCCESS : raise RuntimeError ("Error in GetTemperatureRange: %s " % CameraErrorMsg[result])
		temperature = settings["temperature"]
		if temperature > maxtemp.value or temperature < mintemp.value : raise RuntimeError("Requested temperature is out of range")
	
		# Set the temperature
		if self.CameraLib.CoolerON() != DRV_SUCCESS : raise RuntimeError ("Error in CoolerON")
		if self.CameraLib.SetTemperature (temperature) != DRV_SUCCESS : raise RuntimeError ("Error in SetTemperature")
		
		# Set single scan acquisition mode
		result = self.CameraLib.SetAcquisitionMode(1)
		if result != DRV_SUCCESS : raise RuntimeError ("Error in SetAcquisitionMode: %s " % CameraErrorMsg[result])
		
		# Set Full Vertical Binning readout mode
		result = self.CameraLib.SetReadMode(0)
		if result != DRV_SUCCESS : raise RuntimeError ("Error in SetReadMode: %s " % CameraErrorMsg[result])
		
		# Set exposure time (this must be set at the end)
		exposure_time = float(settings["exposure_time"])/1000.
		result = self.CameraLib.SetExposureTime( ctypes.c_float(exposure_time) )
		if result != DRV_SUCCESS : raise RuntimeError ("Error in SetExposureTime: %s " % CameraErrorMsg[result] )
		
		# Retrieve Acquisition timing and compare with requested values
		exposure = ctypes.c_float(); accumulate = ctypes.c_float(); kinetic = ctypes.c_float()
		result = self.CameraLib.GetAcquisitionTimings( ctypes.byref(exposure), ctypes.byref(accumulate), ctypes.byref(kinetic) )
		if result != DRV_SUCCESS : raise RuntimeError ("Error in GetAcquisitionTimings: %s " % CameraErrorMsg[result] )
		exposure = exposure.value; accumulate = accumulate.value; kinetic = kinetic.value
		
		if not np.isclose(exposure_time,exposure,rtol=1e-3) : 
			raise RuntimeError ("Requested exposure time cannot be set. Nearest available value is %f (s)"% exposure)
			
		return RETURN_SUCCESS
		
	def SetShamrockSettings (self, settings) :
		""" Shamrock settings """
		
		################ Grating ################
		if self.ShamrockLib.ShamrockSetGrating(0, settings["grating"]) != SHAMROCK_SUCCESS : 
			raise RuntimeError ("Error in ShamrockSetGrating")
		grating = ctypes.c_int()
		if self.ShamrockLib.ShamrockGetGrating(0,ctypes.byref(grating)) != SHAMROCK_SUCCESS :
			raise RuntimeError ("Error in ShamrockGetGrating")
		if grating.value != settings["grating"] : raise ValueError ("Grating was not properly set")
		
		################ Print Wavelength Limit ################
		min_wavelength = ctypes.c_float(); max_wavelength = ctypes.c_float() 
		if self.ShamrockLib.ShamrockGetWavelengthLimits(0, grating, ctypes.byref(min_wavelength), ctypes.byref(max_wavelength)) != SHAMROCK_SUCCESS :
			raise ValueError ("Error in ShamrockGetWavelengthLimits")
		print "Grating %d resolves from %f to %f (nm)" % (grating.value, min_wavelength.value, max_wavelength.value)
		
		################ Wavelength ################
		if self.ShamrockLib.ShamrockSetWavelength(0, ctypes.c_float(settings["wavelength"])) != SHAMROCK_SUCCESS :
			raise ValueError ("Error in ShamrockSetWavelength")
		wavelength = ctypes.c_float()
		if self.ShamrockLib.ShamrockGetWavelength(0, ctypes.byref(wavelength)) != SHAMROCK_SUCCESS : raise ValueError ("Error in ShamrockGetWavelength")
		
		if not np.isclose(settings["wavelength"],wavelength.value,rtol=1e-3) : 
			print "Warning: Wavelength %f (nm) requested, but %f (nm) set\n" % (settings["wavelength"], wavelength.value)
		
		################ Slit width ################
		if self.ShamrockLib.ShamrockSetSlit(0, ctypes.c_float(settings["slit_width"])) != SHAMROCK_SUCCESS : 
			raise ValueError ("Error in ShamrockSetSlit")
		slit_width = ctypes.c_float()
		if self.ShamrockLib.ShamrockGetSlit(0, ctypes.byref(slit_width)) != SHAMROCK_SUCCESS : raise ValueError ("Error in ShamrockGetSlit")
		
		if not np.isclose(settings["slit_width"],slit_width.value,rtol=1e-3) :
			raise ValueError ("Slit width was not properly set")
	
		return RETURN_SUCCESS
		
	def GetTemperature (self, arguments=None) :
		"""
		Get current temperature
		"""
		temperature = ctypes.c_int()
		result = self.CameraLib.GetTemperature( ctypes.byref(temperature) ) 
		#if result != DRV_TEMP_STABILIZED : print (CameraErrorMsg[result])
		return temperature.value
	
	def GetWavelengths (self, arguments=None) :
		"""
		Return Shamrock wavelengths calibration
		"""
		npixels = len(self.buffer)
		
		if self.ShamrockLib.ShamrockSetNumberPixels(0,npixels) != SHAMROCK_SUCCESS : raise ValueError("Error in ShamrockSetNumberPixels")
		
		# Get the pixel size for the camera
		pixel_width = ctypes.c_float(); pixel_height = ctypes.c_float()
		if self.CameraLib.GetPixelSize(ctypes.byref(pixel_width), ctypes.byref(pixel_height)) !=  DRV_SUCCESS :
			raise ValueError ("Error in GetPixelSize")
		
		# Specify the pixel width
		if self.ShamrockLib.ShamrockSetPixelWidth(0,pixel_width) != SHAMROCK_SUCCESS : raise ValueError("Error in ShamrockSetPixelWidth")
		
		# Get wave length per pixel
		wavelegths = np.zeros(npixels, dtype=np.dtype(ctypes.c_float))
		wavelegths_ptr = wavelegths.ctypes.data_as(ctypes.POINTER(ctypes.c_float))
		if self.ShamrockLib.ShamrockGetCalibration(0,wavelegths_ptr,npixels) != SHAMROCK_SUCCESS : raise ValueError("Error in ShamrockGetCalibration")
		
		return wavelegths
	
	def StopDevice(self) :
		"""
		Closing the camera
		"""
		if self.CameraLib.CoolerOFF() != DRV_SUCCESS : print ("Error in CoolerOFF")
		if self.CameraLib.ShutDown() != DRV_SUCCESS : print ("Error in ShutDown")
		if self.ShamrockLib.ShamrockClose() != SHAMROCK_SUCCESS : print ("Error in ShamrockClose")
		
		return RETURN_SUCCESS
		
	def AcquiredData (self, arguments=None) :
		"""
		Acquire data from spectrometer
		"""
		# Check status 
		status = ctypes.c_int()
		if self.CameraLib.GetStatus( ctypes.byref(status) ) != DRV_SUCCESS : raise RuntimeError("Error in GetStatus")
		status = status.value	
		
		if status != DRV_IDLE : raise RuntimeError ("Status Error: %s" % CameraErrorMsg[status])
		
		# Camera is ready to accept commands, then begin acquisition
		if self.CameraLib.StartAcquisition() != DRV_SUCCESS : raise RuntimeError("Error in StartAcquisition")
	
		# Watling till acquisition finishes	
		result = self.CameraLib.WaitForAcquisition()
		if result != DRV_SUCCESS : raise RuntimeError ("Error in WaitForAcquisition: %s" % CameraErrorMsg[result])
		
		# Moving the Data into the buffer
		result = self.CameraLib.GetAcquiredData(self.buffer, len(self.buffer)) 
		if result != DRV_SUCCESS : raise RuntimeError ("Error in GetAcquiredData: %s" % CameraErrorMsg[result])
		
		return RETURN_SUCCESS

	def run (self) :
		"""
		Overloaded function provided by BasicDevice
		"""	
		# Initialize the devices
		self.InitializeShamrock(); self.InitializeCamera()

		BasicDevice.run(self)
		
		# Closing the device
		self.StopDevice()
		
########################################################################

class ShamrockSpectrometerTab (HardwareGUIControl) :
	"""
	This class represents a GUI controlling properties of the spectrometer.
	"""
	def __init__(self, parent) :
		HardwareGUIControl.__init__(self, parent)
		
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		################################################
		
		sizer.Add (wx.StaticText(self, label="Exposure time (ms)"), flag=wx.LEFT, border=5)
		exposure_time_ctr = wx.SpinCtrl (self, value="20", min=1, max=1e6)
		exposure_time_ctr.SetLabel("Exposure time")
		sizer.Add (exposure_time_ctr, flag=wx.EXPAND, border=5)

		################ Temperature ################
		sizer.Add (wx.StaticText(self, label="\nTemperature"), flag=wx.LEFT, border=5)
		temperature_ctr = wx.SpinCtrl (self, value="-10", min=-50, max=50)
		temperature_ctr.SetLabel ("Temperature")
		sizer.Add (temperature_ctr, flag=wx.EXPAND, border=10)
		
		################ Shamrock settings ################
		
		################ Grating ################
		sizer.Add (wx.StaticText(self, label="\nGrating #"), flag=wx.LEFT, border=5)
		grating_ctr = wx.SpinCtrl (self, value="1", min=1, max=3)
		grating_ctr.SetLabel("grating")
		sizer.Add (grating_ctr, flag=wx.EXPAND, border=5)
		
		################ Wavelength ################
		sizer.Add (wx.StaticText(self, label="\nWavelength (nm)"), flag=wx.LEFT, border=5)
		wavelength_ctr = wx.SpinCtrl (self, value="610", min=10, max=1000)
		wavelength_ctr.SetLabel("wavelength")
		sizer.Add (wavelength_ctr, flag=wx.EXPAND, border=5)
		
		################ Slit width ################
		sizer.Add (wx.StaticText(self, label="\nSlit width (um)"), flag=wx.LEFT, border=5)
		slit_width_ctr = wx.SpinCtrl (self, value="500", min=0, max=2000)
		slit_width_ctr.SetLabel("slit width")
		sizer.Add (slit_width_ctr, flag=wx.EXPAND, border=5)
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		
		self.CreateSettingsDict()
	
	