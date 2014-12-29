########################################################################
#
#	This module contains classes for controlling and GUI representation of
#	OceanOptics spectrometer 
#
########################################################################

from libs.gui.hardware_control import HardwareGUIControl
from libs.dev.basic_device import BasicDevice

import ctypes, wx, os, multiprocessing
import numpy as np
from libs.dev.consts import * 

########################################################################
#
#	Manager class that communicates with the process where
#	the device resides
#
########################################################################

class ManagerOceanOpticsSpectrometer :
	"""
	Class that manges the OceanOptics spectrometer
	"""
	def __init__ (self) :
		# Create the lock for device 
		self.lock = multiprocessing.Lock()
		# Create a pipe for communication
		self.parent_connection, self.child_connection = multiprocessing.Pipe()
		
	def __del__ (self) :
		self.parent_connection.close()
		self.child_connection.close()
		
	def start(self) :
		"""
		Start the process controlling the spectrometer
		"""
		p = OceanOpticsSpectrometer(self.child_connection)
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
		return self.run("AcquiredData")

	def GetWavelengths (self) : 
		"""
		Get the calibration data.
		"""
		return self.run("GetWavelengths")

########################################################################
#
#	Process where the device resides
#
########################################################################

class OceanOpticsSpectrometer (BasicDevice) :
	"""
	Control spectrometer from a separate process
	"""	

	def Initialize (self, arguments=None) :
		"""
		Initialize spectrometer
		"""
		OOFolder = self.GetProgramFilesDirectory() + "\\Ocean Optics\\OmniDriver\\OOI_HOME\\" 
		self.OOCommon = ctypes.cdll.LoadLibrary(OOFolder + "common32.dll")
		self.OODriver = ctypes.cdll.LoadLibrary(OOFolder + "OmniDriver32.dll")
		
		self.wrapperHandle = self.OODriver.Wrapper_Create()
		
		# Get build number
		print "OceanOptics Build number: %d" %  self.OODriver.Wrapper_getBuildNumber( self.wrapperHandle )

		# Get API value
		apiVersion = self.OOCommon.JString_Create()
		self.OODriver.Wrapper_getApiVersion(self.wrapperHandle,apiVersion)
		print "OO API version %s" % ctypes.c_char_p(self.OOCommon.JString_getASCII(apiVersion)).value
		self.OOCommon.JString_Destroy(apiVersion)

		# How many spectrometers are attached
		numberOfSpectrometersAttached = self.OODriver.Wrapper_openAllSpectrometers(self.wrapperHandle)
		if numberOfSpectrometersAttached == 0 :
			raise ValueError("No OceanOptics spectrometer connected!")
		if numberOfSpectrometersAttached > 1 :
			print "More then one OceanOptics spectrometer connected. The first one will be used."
		self.spectrometerIndex = 0

		# Print spectrometer info
		firmwareVersion = self.OOCommon.JString_Create()
		serialNumber = self.OOCommon.JString_Create()
		spectrometerName = self.OOCommon.JString_Create()
		self.OODriver.Wrapper_getSerialNumber(self.wrapperHandle,self.spectrometerIndex,serialNumber)
		self.OODriver.Wrapper_getName(self.wrapperHandle,self.spectrometerIndex,spectrometerName)
		self.OODriver.Wrapper_getFirmwareVersion(self.wrapperHandle,self.spectrometerIndex,firmwareVersion)
		
		print "OO Spectrometer type: %s   serial number: %s  firmware version: %s\n" % (
			ctypes.c_char_p( self.OOCommon.JString_getASCII(spectrometerName) ).value,
			ctypes.c_char_p( self.OOCommon.JString_getASCII(serialNumber) ).value,
			ctypes.c_char_p( self.OOCommon.JString_getASCII(firmwareVersion) ).value
		)

		self.OOCommon.JString_Destroy(firmwareVersion)
		self.OOCommon.JString_Destroy(serialNumber)
		self.OOCommon.JString_Destroy(spectrometerName)
		
		# Create the buffer for saving wavelengths and spectrum
		self.bufferHandle = self.OOCommon.DoubleArray_Create()
		# Fill the array to find out the shape
		self.OODriver.Wrapper_getWavelengths(self.wrapperHandle,self.spectrometerIndex,self.bufferHandle)
		numberOfPixels = self.OOCommon.DoubleArray_getLength(self.bufferHandle)
		bufferValues = self.OOCommon.DoubleArray_getDoubleValues(self.bufferHandle)
		
		# Convert buffer to numpy array
		self.buffer = np.ctypeslib.as_array( ctypes.cast(bufferValues, ctypes.POINTER(ctypes.c_double)), (numberOfPixels,) ) 
		
		return RETURN_SUCCESS

	def StopDevice (self, arguments=None) :
		"""
		Closing spectrometer 
		"""
		# Cleaning
		self.OOCommon.DoubleArray_Destroy(self.bufferHandle)
		# Closing
		self.OODriver.Wrapper_closeAllSpectrometers(self.wrapperHandle)
		self.OODriver.Wrapper_Destroy(self.wrapperHandle)
		
		return RETURN_SUCCESS

	def SetSettings (self, arguments) :
		"""
		Set acquisition parameters
		"""
		# Set integration time
		arguments["exposure_time"] *= 1000 # Convert ms to us
		minimumAllowedIntegrationTime = self.OODriver.Wrapper_getMinimumIntegrationTime(self.wrapperHandle,self.spectrometerIndex)
		if arguments["exposure_time"] <  minimumAllowedIntegrationTime :
			raise RuntimeError("Warning: Requested exposure time is smaller then the minimum allowed value. The minimal value is %d usec." % minimumAllowedIntegrationTime)
		self.OODriver.Wrapper_setIntegrationTime(self.wrapperHandle,self.spectrometerIndex,arguments["exposure_time"])
		
		# Set averages -- number of consecutive acquisitions
		self.OODriver.Wrapper_setScansToAverage(self.wrapperHandle, self.spectrometerIndex, arguments["scans_average"])
		
		#self.OODriver.Wrapper_setBoxcarWidth(self.wrapperHandle,self.spectrometerIndex,0)
		#self.OODriver.Wrapper_setScansToAverage(self.wrapperHandle,self.spectrometerIndex,1)
		#self.OODriver.Wrapper_setCorrectForElectricalDark(self.wrapperHandle,self.spectrometerIndex,0)
		
		# Select the wavelength interval
		self.GetWavelengths()
		if  self.buffer.min() >  arguments["start_wavelength"] or self.buffer.max() <  arguments["end_wavelength"] :
			raise ValueError ("\nWavelength range (start_wavelength, end_wavelength) is incorrect, sported values are (%f, %f) nm" % ( self.buffer.min(),  self.buffer.max()) )
		
		self.spectral_interval = np.nonzero( (arguments["start_wavelength"] <= self.buffer)&(self.buffer < arguments["end_wavelength"]) )
	
		
		
		return RETURN_SUCCESS
	
	def GetWavelengths (self, arguments=None) :
		"""
		Get wavelength range
		"""
		self.OODriver.Wrapper_getWavelengths(self.wrapperHandle,self.spectrometerIndex,self.bufferHandle)
		
		try : return self.buffer[self.spectral_interval]
		except AttributeError : return self.buffer
	
	def AcquiredData (self, arguments=None) :
		"""
		Acquire spectrum from spectrometer
		"""
		self.OODriver.Wrapper_getSpectrum(self.wrapperHandle,self.spectrometerIndex,self.bufferHandle)
		
		if self.OODriver.Wrapper_isSaturated(self.wrapperHandle,self.spectrometerIndex) :
			print "Warning: OcenOptics spectrometer is saturated!"
			
		try : return self.buffer[self.spectral_interval]
		except AttributeError : return self.buffer
		
	def run (self) :
		"""
		Overloaded function provided by BasicDevice.
		"""	
		# Initialize the devices
		self.Initialize()

		BasicDevice.run(self)

		# Closing the device
		self.StopDevice()
		
########################################################################

class OceanOpticsSpectrometerTab (HardwareGUIControl) :
	"""
	This class represents a GUI controlling properties of the spectrometer.
	"""
	def __init__(self, parent, dev) :
		HardwareGUIControl.__init__(self, parent, dev)
		
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		################################################
		
		sizer.Add (wx.StaticText(self, label="Exposure time (millisecond)"), flag=wx.LEFT, border=5)
		exposure_time_ctr = wx.SpinCtrl (self, value="20", min=1, max=1e6)
		exposure_time_ctr.SetLabel("Exposure time")
		sizer.Add (exposure_time_ctr, flag=wx.EXPAND, border=5)
		
		################################################
		
		sizer.Add (wx.StaticText(self, label="Scans to average"), flag=wx.LEFT, border=5)
		scans_average_ctr = wx.SpinCtrl (self, value="1", min=1, max=1e6)
		scans_average_ctr.SetLabel("scans_average")
		sizer.Add (scans_average_ctr, flag=wx.EXPAND, border=5)
		
		################################################
		
		sizer.Add (wx.StaticText(self, label="\nStart wavelength (nm)"), flag=wx.LEFT, border=5)
		start_wavelength_ctr = wx.SpinCtrl (self, value="196", min=1, max=1e6)
		#start_wavelength_ctr = wx.SpinCtrl (self, value="580", min=1, max=1e6)
		start_wavelength_ctr.SetLabel("start_wavelength")
		sizer.Add (start_wavelength_ctr, flag=wx.EXPAND, border=5)
		
		################################################
		
		sizer.Add (wx.StaticText(self, label="End wavelength (nm)"), flag=wx.LEFT, border=5)
		end_wavelength_ctr = wx.SpinCtrl (self, value="1118", min=1, max=1e6)
		#end_wavelength_ctr = wx.SpinCtrl (self, value="680", min=1, max=1e6)
		end_wavelength_ctr.SetLabel("end_wavelength")
		sizer.Add (end_wavelength_ctr, flag=wx.EXPAND, border=5)
		
		############  Parameters for displaying spectrum #############
		# Spacer
		sizer.Add (wx.StaticText(self), flag=wx.LEFT, border=5)

		sb_sizer = wx.StaticBoxSizer( wx.StaticBox(self, label="Plotting options"),  wx.VERTICAL )
		
		fix_vert_ax_ctrl = wx.CheckBox(self, label="Fix vertical axis")
		fix_vert_ax_ctrl.__label__ = "fix_vertical_axis"
		sb_sizer.Add (fix_vert_ax_ctrl, flag=wx.EXPAND, border=5)
		
		sb_sizer.Add (wx.StaticText(self, label="max val"), flag=wx.LEFT, border=5)
		max_val_ctr = wx.SpinCtrl (self, value="50000", min=1, max=1e6)
		max_val_ctr.SetLabel("vertical_axis_max_val")
		sb_sizer.Add (max_val_ctr, flag=wx.EXPAND, border=5)
		
		sb_sizer.Add (wx.StaticText(self, label="min val"), flag=wx.LEFT, border=5)
		min_val_ctr = wx.SpinCtrl (self, value="100", min=1, max=1e6)
		min_val_ctr.SetLabel("vertical_axis_min_val")
		sb_sizer.Add (min_val_ctr, flag=wx.EXPAND, border=5)
		
		sizer.Add (sb_sizer, flag=wx.EXPAND, border=5)
		################################################
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		
		self.CreateSettingsDict()
	
	