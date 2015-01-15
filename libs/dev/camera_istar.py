########################################################################
#
#	This file contains the class for controlling 
#	the Andor iStar 334T camera as well as 
#	GUI panel for  representing its settings 
#
########################################################################

from libs.gui.hardware_control import HardwareGUIControl
from libs.dev.basic_device import BasicDevice
from libs.dev.consts import * 


import multiprocessing, ctypes, wx
import numpy as np
from libs.dev.consts import * 

########################################################################
#
#	Managing class that communicates with the process where
#	the device resides
#
########################################################################

class ManagerIStarCamera :
	"""
	Class that manages Andor IStar camera
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
		# Create the image buffer 
		self.spectrometer_buffer = IStarCamera.AllocateBuffer()
		p = IStarCamera(self.child_connection, self.spectrometer_buffer)
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
		if self.run("AcquiredData") == RETURN_FAIL : print "IStar Camera Acquisition failed"
		new_shape = self.run("GetBufferShape")
		spectrum = np.frombuffer(self.spectrometer_buffer.get_obj(), 
				dtype=np.dtype(ctypes.c_long), count=np.prod(new_shape))
				
		if new_shape[0] == 1 or new_shape[1] == 1: return spectrum
		else : return spectrum.reshape( new_shape )
		
	def GetTemperature (self) : return self.run("GetTemperature")
	
########################################################################

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
DRV_LOAD_FIRMWARE_ERROR = 20096
DRV_ERROR_NOCAMERA 	= 20990

########################################################################

# Transcript of Error messages
SpectrometerErrorMsg = { DRV_VXDNOTINSTALLED : "VxD not loaded",
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
	DRV_TEMP_NOT_STABILIZED : "Temperature reached but not stabilized", DRV_LOAD_FIRMWARE_ERROR : "Load firmware error"
}


class IStarCamera (BasicDevice):
	"""
	Control spectrometer hardware 
	"""
	def __init__ (self,  pipe, spectrometer_buffer) :
		"""
		<pipe> is for communicating.
		<spectrometer_buffer> is a multiprocessing.Array where the spectra will be saved.
		"""
		BasicDevice.__init__(self, pipe)
		# saving the buffer where the spectrum will be saved
		self.buffer = spectrometer_buffer
	
	@classmethod
	def GetAndorDriverFolder (cls) :
		"""
		Return the directory containing Andor driver
		"""
		return cls.GetProgramFilesDirectory() + "/Andor SOLIS/Drivers/"
	
	@classmethod
	def __Get_Lib_Xpixels_Ypixels__ (cls) :
		"""
		Return the reference to driver's DLL, number of pixels in X and Y dimension
		"""
		# Loading the spectrometer driver
		lib = ctypes.WinDLL (cls.GetAndorDriverFolder() + "atmcd32d.dll")
		
		# Verifying that there is a single Andor camera
		totalCameras = ctypes.c_long()
		if lib.GetAvailableCameras( ctypes.byref(totalCameras) ) != DRV_SUCCESS : 
			raise RuntimeError ("Error in GetAvailableCameras")

		if totalCameras.value > 1 : 
			raise RuntimeError ("More than one Andor camera is present")

		# Initialize the camera	
		result = lib.Initialize( cls.GetAndorDriverFolder() ) 
		if result != DRV_SUCCESS : raise RuntimeError ("Error in Initialize: %s " % SpectrometerErrorMsg[result])
		
		# Find out the number of pixels for figures
		xpixels = ctypes.c_int(); ypixels = ctypes.c_int()
		if lib.GetDetector( ctypes.byref(xpixels), ctypes.byref(ypixels) ) != DRV_SUCCESS : 
			raise RuntimeError ("Error in GetDetector")
		return lib, xpixels.value, ypixels.value
		
	@classmethod
	def AllocateBuffer (cls) :
		"""
		This static method allocates buffer that corresponds to Full Vertical Binning readout mode
		"""
		# Find out buffer's size
		lib,  xpixels, ypixels = cls.__Get_Lib_Xpixels_Ypixels__()
		# Closing camera's driver
		lib.ShutDown()
		del lib
		return multiprocessing.Array (ctypes.c_long, xpixels*ypixels )
	
	def InitializeDevice (self, arguments=None) :
		"""
		Initiate hardware
		"""
		self.lib, self.xpixels, self.ypixels = self.__Get_Lib_Xpixels_Ypixels__()
		return RETURN_SUCCESS
			
	def GetBufferShape (self, arguments=None) :
		return self.buffer_shape
	
	def SetSettings (self, settings) :
		"""
		Assign settings 
		"""		
		# check the temperature range
		mintemp = ctypes.c_int(); maxtemp = ctypes.c_int()
		result = self.lib.GetTemperatureRange( ctypes.byref(mintemp), ctypes.byref(maxtemp) )
		if result != DRV_SUCCESS : raise RuntimeError ("Error in GetTemperatureRange: %s " % SpectrometerErrorMsg[result])
		temperature = settings["temperature"]
		if temperature > maxtemp.value or temperature < mintemp.value : raise RuntimeError("Requested temperature is out of range")
		
		# Set the temperature
		if self.lib.CoolerON() != DRV_SUCCESS : raise RuntimeError ("Error in CoolerON")
		if self.lib.SetTemperature (temperature) != DRV_SUCCESS : raise RuntimeError ("Error in SetTemperature")
		
		# Set accumulate acquisition mode
		result = self.lib.SetAcquisitionMode(2)
		if result != DRV_SUCCESS : raise RuntimeError ("Error in SetAcquisitionMode: %s " % SpectrometerErrorMsg[result])
		
		accumulate_number = settings["scans_to_accumulate"]
		if self.lib.SetNumberAccumulations(accumulate_number) != DRV_SUCCESS : 
			raise RuntimeError ("Error in SetNumberAccumulations")
		
		if settings["full_vertical_binning"] :
			# Set Full Vertical Binning readout mode
			result = self.lib.SetReadMode(0)
			if result != DRV_SUCCESS : raise RuntimeError ("Error in SetReadMode: %s " % SpectrometerErrorMsg[result])
			
			self.buffer_shape = (1, self.xpixels)
			self.buffer_size = self.xpixels
		
		if settings["image"] :
			# Set Full resolution Image readout more
			result = self.lib.SetReadMode(4)
			if result != DRV_SUCCESS : raise RuntimeError ("Error in SetReadMode: %s " % SpectrometerErrorMsg[result])
			
			x_bin = settings["x_bin"]; y_bin = settings["y_bin"]
			x_min = settings["x_min"]; x_max = settings["x_max"]
			y_min = settings["y_min"]; y_max = settings["y_max"]
			
			self.buffer_shape = ( int((y_max-y_min+1)//y_bin), int((x_max-x_min+1)//x_bin) )
			self.buffer_size = self.buffer_shape[0]*self.buffer_shape[1]
		
			# Making `x_max` and `y_max` be compatible with selected binning
			y_max = self.buffer_shape[0]*y_bin + y_min - 1
			x_max = self.buffer_shape[1]*x_bin + x_min - 1
		
			result = self.lib.SetImage(x_bin,y_bin,x_min,x_max,y_min,y_max)
			if result != DRV_SUCCESS : raise RuntimeError ("Error in SetImage")
		
		# Set trigger mode
		if settings["internal"] :
			result = self.lib.SetTriggerMode(0)
			if result != DRV_SUCCESS : raise RuntimeError ("Error in SetTriggerMode: %s " % SpectrometerErrorMsg[result] )
			
			result = self.lib.SetGateMode (3) 
			if result != DRV_SUCCESS : raise RuntimeError ("Error in SetGateMode: %s " % SpectrometerErrorMsg[result] )
			
		if settings["external"] :	
			result = self.lib.SetTriggerMode(1)
			if result != DRV_SUCCESS : raise RuntimeError ("Error in SetTriggerMode: %s " % SpectrometerErrorMsg[result] )
			
			# Set the parameters of on-chip integration to eliminate the scattering
			if self.lib.SetDDGIOC(1) != DRV_SUCCESS : raise RuntimeError("Error in SetDDGIOC")
			if self.lib.SetDDGIOCTrigger(1) != DRV_SUCCESS : raise RuntimeError("Error in SetDDGIOCTrigger")
			
			result = self.lib.SetGateMode (5) 
			if result != DRV_SUCCESS : raise RuntimeError ("Error in SetGateMode: %s " % SpectrometerErrorMsg[result] )
			
			# Setting DDG gating parameters
			self.SetDDGGateTime( (settings["ddg_delay_time"], settings["ddg_width"]) )
			
			# Number of pulses to squire 
			if self.lib.SetDDGIOCNumber( ctypes.c_ulong(settings["iocnumber"]) ) != DRV_SUCCESS : raise  RuntimeError("Error in SetDDGIOCNumber")
			
		# Set gain
		if self.lib.SetMCPGain(settings["gain"]) != DRV_SUCCESS : raise RuntimeError ("Error in SetMCPGain")
		# Making sure the set gain is correct
		gain = ctypes.c_int()
		if self.lib.GetMCPGain( ctypes.byref(gain) ) != DRV_SUCCESS : raise RuntimeError ("Error in GetMCPGain")
		if gain.value != settings["gain"] : raise RuntimeError ("Error: Gain is not set")
		
		# Set exposure time (this must be set at the end)
		exposure_time = float(settings["exposure_time"])/1000.
		result = self.lib.SetExposureTime( ctypes.c_float(exposure_time) )
		if result != DRV_SUCCESS : raise RuntimeError ("Error in SetExposureTime: %s " % SpectrometerErrorMsg[result] )
		
		# Retrieve Acquisition timing and compare with requested values
		exposure = ctypes.c_float(); accumulate = ctypes.c_float(); kinetic = ctypes.c_float()
		result = self.lib.GetAcquisitionTimings( ctypes.byref(exposure), ctypes.byref(accumulate), ctypes.byref(kinetic) )
		if result != DRV_SUCCESS : raise RuntimeError ("Error in GetAcquisitionTimings: %s " % SpectrometerErrorMsg[result] )
		exposure = exposure.value; accumulate = accumulate.value; kinetic = kinetic.value
		
		if not np.isclose(exposure_time,exposure) : 
			raise RuntimeError ("Requested exposure time cannot be set. Nearest available value is %f (s)"% exposure)
		
		#if not np.isclose(accumulate_number,accumulate) : 
		#	raise RuntimeError ("Requested accumulation number cannot be set. Nearest available value is %f (s)"% accumulate)
			
		return RETURN_SUCCESS
	
	def SetDDGGateTime (self, arguments) :
		"""
		Set DDF gating parameters (in nanoseconds) 
		"""
		delay_set = ctypes.c_ulonglong(1000*arguments[0])
		width_set = ctypes.c_ulonglong(1000*arguments[1])
		if self.lib.SetDDGGateTime(delay_set, width_set)  != DRV_SUCCESS : raise RuntimeError("Error in SetDDGGateTime")
		# Check whether the values are set properly
		Delay = ctypes.c_ulonglong(); Width = ctypes.c_ulonglong()
		if self.lib.GetDDGGateTime( ctypes.byref(Delay), ctypes.byref(Width) ) != DRV_SUCCESS : raise RuntimeError("Error in GetDDGGateTime")
		if Delay.value != delay_set.value or Width.value != width_set.value : raise RuntimeError("DDG gating parameters are not properly set")
		
		return RETURN_SUCCESS
	
	def GetTemperature (self, arguments=None) :
		"""
		Get current temperature
		"""
		temperature = ctypes.c_int()
		result = self.lib.GetTemperature( ctypes.byref(temperature) ) 
		#if result != DRV_TEMP_STABILIZED : print (SpectrometerErrorMsg[result])
		return temperature.value
	
	def StopDevice(self, arguments=None) :
		"""
		Closing the camera
		"""
		if self.lib.ShutDown() != DRV_SUCCESS : print ("Error in ShutDown")
		if self.lib.CoolerOFF() != DRV_SUCCESS : print ("Error in CoolerOFF")
		return RETURN_SUCCESS

	def AcquiredData (self, arguments=None) :
		"""
		Acquire data from spectrometer
		"""
		# Check status 
		status = ctypes.c_int()
		if self.lib.GetStatus( ctypes.byref(status) ) != DRV_SUCCESS : raise RuntimeError("Error in GetStatus")
		
		if status.value != DRV_IDLE : 
			raise RuntimeError ("Status Error: %s" % SpectrometerErrorMsg[status.value])
		
		# Camera is ready to accept commands, then begin acquisition
		if self.lib.StartAcquisition() != DRV_SUCCESS : raise RuntimeError("Error in StartAcquisition")
	
		# Watling for the end of data accumulation
		status.value = DRV_ACQUIRING	
		while status.value == DRV_ACQUIRING :
			# Watling till acquisition finishes	
			result = self.lib.WaitForAcquisition()
			if result != DRV_SUCCESS : raise RuntimeError ("Error in WaitForAcquisition: %s" % SpectrometerErrorMsg[result])
			
			# Check status, once more
			self.lib.GetStatus( ctypes.byref(status) ) 
			
		# Moving the Data into the buffer
		self.buffer.acquire()
		result = self.lib.GetAcquiredData(self.buffer.get_obj(), self.buffer_size) 
		self.buffer.release()
		if result != DRV_SUCCESS : raise RuntimeError ("Error in GetAcquiredData: %s" % SpectrometerErrorMsg[result])
		
		return RETURN_SUCCESS
	
	def run (self) :
		"""
		Overloaded function provided by BasicDevice
		"""
		# Start the hardware
		self.InitializeDevice()
		
		BasicDevice.run(self)
		
		self.StopDevice()
		
		
########################################################################

class IStarCameraTab (HardwareGUIControl) :
	"""
	This class represents a GUI controlling properties of IStar camera .
	"""
	def __init__(self, parent, dev) :
		HardwareGUIControl.__init__(self, parent, dev)
		
		xpixels, ypixels = 1024, 1024
		
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		################ Readout mode settings ################
		boxsizer_ = wx.StaticBoxSizer( wx.StaticBox(self, label="Readout mode"), wx.HORIZONTAL)
		
		boxsizer_.Add ( wx.RadioButton(self, label="Full vertical binning", style=wx.RB_GROUP), flag=wx.EXPAND, border=5)
		boxsizer_.Add ( wx.RadioButton(self, label="Image"), flag=wx.EXPAND, border=5)
		
		sizer.Add (boxsizer_, flag=wx.EXPAND, border=10)
		
		################ Image dimensions ################
		boxsizer_ = wx.StaticBoxSizer( wx.StaticBox(self, label="Image resolution (image mode only)"), wx.VERTICAL)
		
		boxsizer_.Add ( wx.StaticText(self, label="X min"), flag=wx.LEFT, border=5)
		x_min_ctrl = wx.SpinCtrl (self, value="1", min=1, max=xpixels)
		x_min_ctrl.SetLabel ("x_min")
		boxsizer_.Add (x_min_ctrl, flag=wx.EXPAND, border=5)
		
		boxsizer_.Add ( wx.StaticText(self, label="X max"), flag=wx.LEFT, border=5)
		x_max_ctrl = wx.SpinCtrl (self, value=str(xpixels), min=1, max=xpixels)
		x_max_ctrl.SetLabel ("x_max")
		boxsizer_.Add (x_max_ctrl, flag=wx.EXPAND, border=5)
		
		boxsizer_.Add ( wx.StaticText(self, label="Y min"), flag=wx.LEFT, border=5)
		y_min_ctrl = wx.SpinCtrl (self, value="1", min=1, max=ypixels)
		y_min_ctrl.SetLabel ("y_min")
		boxsizer_.Add (y_min_ctrl, flag=wx.EXPAND, border=5)
		
		boxsizer_.Add ( wx.StaticText(self, label="Y max"), flag=wx.LEFT, border=5)
		y_max_ctrl = wx.SpinCtrl (self, value=str(ypixels), min=1, max=ypixels)
		y_max_ctrl.SetLabel ("y_max")
		boxsizer_.Add (y_max_ctrl, flag=wx.EXPAND, border=5)
		
		# Image binning
		boxsizer_.Add ( wx.StaticText(self, label="X binning"), flag=wx.LEFT, border=5)
		x_bin_ctrl = wx.SpinCtrl (self, value="1", min=1, max=xpixels)
		x_bin_ctrl.SetLabel ("x_bin")
		boxsizer_.Add (x_bin_ctrl, flag=wx.EXPAND, border=5)
		
		boxsizer_.Add ( wx.StaticText(self, label="Y binning"), flag=wx.LEFT, border=5)
		y_bin_ctrl = wx.SpinCtrl (self, value="1", min=1, max=ypixels)
		y_bin_ctrl.SetLabel ("y_bin")
		boxsizer_.Add (y_bin_ctrl, flag=wx.EXPAND, border=5)
		
		sizer.Add (boxsizer_, flag=wx.EXPAND, border=10)
		################################################
		
		sizer.Add (wx.StaticText(self, label="\nExposure time (ms)"), flag=wx.LEFT, border=5)
		exposure_time_ctr = wx.SpinCtrl (self, value="20", min=1, max=1e6)
		exposure_time_ctr.SetLabel( "Exposure time")
		sizer.Add (exposure_time_ctr, flag=wx.EXPAND, border=5)
		
		sizer.Add (wx.StaticText(self, label="Scans to accumulate"), flag=wx.LEFT, border=5)
		scans_to_average_ctr = wx.SpinCtrl (self, value="1", min=1, max=1e6)
		scans_to_average_ctr.__label__ = "scans_to_accumulate"
		sizer.Add (scans_to_average_ctr, flag=wx.EXPAND, border=5)
		
		################### Gain ################################
		sizer.Add (wx.StaticText(self, label="\nGain"), flag=wx.LEFT, border=5)
		gain_ctr = wx.SpinCtrl (self, value="1", min=1, max=1e6)
		gain_ctr.SetLabel ("Gain")
		sizer.Add (gain_ctr, flag=wx.EXPAND, border=5)
		
		################ Trigger modes ################
		boxsizer_ = wx.StaticBoxSizer(wx.StaticBox(self, label="Trigger mode"), wx.HORIZONTAL)	
	
		boxsizer_.Add ( wx.RadioButton(self, label="Internal", style=wx.RB_GROUP), flag=wx.EXPAND, border=5)
		boxsizer_.Add ( wx.RadioButton(self, label="External"), flag=wx.EXPAND, border=5)
		
		sizer.Add (boxsizer_, flag=wx.EXPAND, border=10)
		
		################ Parameters of DDG used only when the external trigger is on ################
		
		sizer.Add (wx.StaticText(self, label="\nDDG gate time delay (ns)"), flag=wx.LEFT, border=5)
		DDG_delay_time_ctr = wx.SpinCtrl (self, value="0", min=0, max=1e9)
		DDG_delay_time_ctr.SetLabel ("DDG delay time")
		sizer.Add (DDG_delay_time_ctr, flag=wx.EXPAND, border=10)
		
		sizer.Add (wx.StaticText(self, label="DDG gate width time (ns)"), flag=wx.LEFT, border=5)
		DDG_width_ctr = wx.SpinCtrl (self, value="4", min=0, max=1e9)
		DDG_width_ctr.SetLabel ("DDG width")
		sizer.Add (DDG_width_ctr, flag=wx.EXPAND, border=10)
		
		################ Number of pulses for the integrate on chip option ################ 
		sizer.Add (wx.StaticText(self, label="\nNumber of pulses to integrate on chip"), flag=wx.LEFT, border=5)
		IOC_num_ctr = wx.SpinCtrl (self, value="1", min=0, max=1e6)
		IOC_num_ctr.SetLabel ("IOCNumber")
		sizer.Add (IOC_num_ctr, flag=wx.EXPAND, border=10)
		
		################ Temperature ################
		sizer.Add (wx.StaticText(self, label="\nTemperature"), flag=wx.LEFT, border=5)
		temperature_ctr = wx.SpinCtrl (self, value="0", min=-50, max=50)
		temperature_ctr.SetLabel ("Temperature")
		sizer.Add (temperature_ctr, flag=wx.EXPAND, border=10)
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		
		self.CreateSettingsDict()
	
	