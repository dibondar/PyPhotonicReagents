"""

This module contains classes for controlling and GUI representation of
PicoHarp (https://www.picoquant.com/products/category/tcspc-and-time-tagging-modules/picoharp-300-stand-alone-tcspc-module-with-usb-interface)

"""

from libs.gui.hardware_control import HardwareGUIControl
from libs.dev.basic_device import BasicDevice
from libs.gui.multi_state_button import MultiStateButton

import ctypes
import wx
import multiprocessing
import functools
import numpy as np
from libs.dev.consts import * 


########################################################################
#
#	Manager class that communicates with the process where
#	the device resides
#
########################################################################

class ManagerPicoHarp :
	"""
	Class that manges the PicoHarp 
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
		Start the process controlling PicoHarp
		"""
		# Create the image buffer 
		self.histogram_buffer = PicoHarp.AllocateBuffer()
		p = PicoHarp(self.child_connection, self.histogram_buffer)
		p.start()
		return p
	
	def StartHistogramMeas (self, tacq=None) :
		"""
		Start acquiring histogram in parallel fashion
		"""
		self.lock.acquire()
		self.parent_connection.send( ("GetHistogram", tacq) )
	
	def StopHistogramMeas (self) :
		"""
		Wait till the histogram acquisition finished and return a reference to histogram 
		"""
		result = self.parent_connection.recv()
		self.lock.release()
		if result === RETURN_FAIL : 
			print "PicoHarp acquisition failed"
		return  np.frombuffer(self.histogram_buffer.get_obj(), dtype=np.uint)
		
	def GetHistogram (self, tacq=None) :
		"""
		Acquire histogram sequentially
		"""
		self.StartHistogramMeas(tacq)
		return self.StopHistogramMeas()
	
	def run(self, command, arguments=None) :
		"""
		Send the command to the spectrometer through the pipe
		"""
		self.lock.acquire()
		self.parent_connection.send( (command, arguments) )
		result = self.parent_connection.recv()
		self.lock.release()
		return result
	
	def __getattr__ (self, name) :
		"""
		Redirect all other request to the method `run`
		"""
		return functools.partial( self.run, name )
		
###########################################################################
# constants from Phdefin.h

HISTCHAN  =  65536 # number of histogram channel
FLAG_OVERFLOW = 0x0040
		
###########################################################################

class PicoHarp (BasicDevice):
	"""
	Control spectrometer hardware 
	"""
	def __init__ (self,  pipe, histogram_buffer) :
		"""
		<pipe> is for communicating.
		<histogram_buffer> is a multiprocessing.Array where histogram is saved.
		"""
		BasicDevice.__init__(self, pipe)
		# saving the buffer where the spectrum will be saved
		self.buffer = histogram_buffer
		
	@classmethod
	def AllocateBuffer (cls) :
		return multiprocessing.Array (ctypes.c_uint, HISTCHAN)
		
	def Initialize (self, arguments=None) :
		
		# Open library 
		self.lib = ctypes.WinDLL ("C:/Windows/System32/phlib.dll")
	
		# Allocate string buffer
		self.str_p =  ctypes.create_string_buffer (500)
		
		# Extract the Pico Harp library version
		lib.PH_GetLibraryVersion (str_p)
		pico_harp_lib_version_dev = "2.3"
		if str_p.value <> pico_harp_lib_version_dev :
			print "Warning: The current code was developed for PicoHarp Library version %s. The version of the currently installed library is %s" \
				% (pico_harp_lib_version_dev, str_p.value)
				
		# Pico harp device ID to be utilized
		self.DEV_INDX = ctypes.c_int(0)
		
		# Closing Pico Harp, just in case if it was already open 
		self.lib.PH_CloseDevice (self.DEV_INDX)
		
		# Open the first PicoHarp
		if self.lib.PH_OpenDevice (DEV_INDX, str_p) < 0 :
			raise RuntimeError ("The first Pico Harp could not be opened")
		
		print "PicoHarp %s opened" % str_p.value
		
		return RETURN_SUCCESS
		
	def GetWarnigns (self, arguments=None) :
		"""
		This method checks for warnings, and if there are any it prints them.
		"""
		warnings = self.lib.PH_GetWarnings (self.DEV_INDX)
		if warnings < 0 :
			raise RuntimeError ("PH_GetWarnings failed")
		elif warnings > 0 :
			# Get the warning message
			if lib.PH_GetWarningsText (self.DEV_INDX, self.str_p, warnings) < 0 :
				raise RuntimeError("PH_GetWarningsText failed")
			print self.str_p.value
			
		return self.str_p.value

	def SetSettings (self, settings) :
		"""
		Assign settings
		"""
		self.MODE_HIST = ctypes.c_int(0)

		# remember acquisition time
		self.tacq = settings["tacq"]
		
		# Initializing the Pico Harp device with device ID <DEV_INDEX> in the Histogram mode 
		if self.lib.PH_Initialize (self.DEV_INDX, self.MODE_HIST) < 0 :
			raise RuntimeError ("The first Pico Harp device could not be initiated in the Histogram mode") 
		
		if self.lib.PH_Calibrate (self.DEV_INDX) < 0 :
			raise RuntimeError ("Pico Harp calibration failed")
		
		if self.lib.PH_SetSyncDiv (self.DEV_INDX, settings["syncdiv"]) < 0 :
			raise RuntimeError ("PH_SetSyncDiv failed") 

		if self.lib.PH_SetCFDLevel (self.DEV_INDX, ctypes.c_int(0), settings["cfdlevel0"]) < 0 :
			raise RuntimeError ("PH_SetCFDLevel 0 failed")

		if self.lib.PH_SetCFDLevel (self.DEV_INDX, ctypes.c_int(1), settings["cfdlevel1"]) < 0 :
			raise RuntimeError ("PH_SetCFDLevel 1 failed") 

		if self.lib.PH_SetCFDZeroCross (self.DEV_INDX, ctypes.c_int(0), settings["cfdzerox0"]) < 0 :
			raise RuntimeError ("PH_SetCFDZeroCross 0 failed") 

		if self.lib.PH_SetCFDZeroCross (self.DEV_INDX, ctypes.c_int(1), settings["cdfzerox1"]) < 0 :
			raise RuntimeError ("PH_SetCFDZeroCross 1 failed") 

		if self.lib.PH_SetRange (self.DEV_INDX, settings["range"]) < 0 :
			raise RuntimeError ("PH_SetRange failed") 

		if self.lib.PH_SetOffset (self.DEV_INDX, settings["offset"]) < 0 :
			raise RuntimeError ("PH_SetOffset failed") 

		if self.lib.PH_SetStopOverflow (self.DEV_INDX, ctypes.c_int(1), ctypes.c_int(65535)) < 0 :
			raise RuntimeError ("PH_SetStopOverflow failed") 

		resolution = self.lib.PH_GetResolution (self.DEV_INDX)
		if resolution < 0 :
			raise RuntimeError ("PH_GetResolution failed")

		print "PicoHarp initialized with resolution is %d\n" % resolution
			
		# Note: after Init or SetSyncDiv you must allow 100 ms for valid new count rate readings
		time.sleep (0.2)
		
		print "Input0 count rate = %d / s\nInput1 count rate = %d / s\n" % self.GetCountRates()
		
		self.GetWarnigns()
		
		return RETURN_SUCCESS
		
	def GetCountRates (self, arguments=None) :
		"""
		Return count rates for both channels
		"""
		return 	self.lib.PH_GetCountRate (self.DEV_INDX, ctypes.c_int(0)),
				self.lib.PH_GetCountRate (self.DEV_INDX, ctypes.c_int(1))
		
	def GetHistogram (self, arguments=None) :
		"""
		Acquire histogram  
		"""
		# Saving the acquisition time, if provided
		if arguments :
			self.tacq = arguments[0]

		if self.lib.PH_ClearHistMem (self.DEV_INDX, ctypes.c_int(0)) < 0 : # always use Block 0 if not Routing
			raise RuntimeError ("PH_ClearHistMem failed")

		if self.lib.PH_StartMeas (self.DEV_INDX, self.tacq) < 0 :
			raise RuntimeError ("PH_StartMeas failed") 

		# Measuring for self.tacq ms
		ctcdone = 0
		while ctcdone == 0 : 
			ctcdone = lib.PH_CTCStatus (self.DEV_INDX)

		if self.lib.PH_StopMeas (self.DEV_INDX) < 0 :
			raise RuntimeError ("PH_StopMeas failed") 
				
		# Retrieving the histogram
		self.buffer.acquire()
		ret = self.lib.PH_GetBlock (self.DEV_INDX, self.buffer.get_obj(), ctypes.c_int(0)) 
		self.buffer.release()
		if ret < 0 :
			raise RuntimeError ("PH_GetBlock failed\n")
	
		flags = self.lib.PH_GetFlags (self.DEV_INDX)
		if flags < 0 : 
			raise RuntimeError ("PH_GetFlags failed")
				
		if flags & FLAG_OVERFLOW : 
			print "\n\tHistogram overflow!\n"
	
		return RETURN_SUCCESS
	
	####################################################################################
	def UpdateOffset (self, arguments) :
		"""
		Interactively update
		"""
		if self.lib.PH_SetOffset (self.DEV_INDX, arguments[0]) < 0 :
			print "Interactive update of PH_SetOffset failed" 
			return RETURN_FAIL
		else :
			return RETURN_SUCCESS

	def UpdateCFDLevel0 (self, arguments) :
		"""
		Interactively update
		"""
		if self.lib.PH_SetCFDLevel (self.DEV_INDX, ctypes.c_int(0), arguments[0]) < 0 :
			print "Interactive update of PH_SetCFDLevel 0 failed"
			return RETURN_FAIL
		else :
			return RETURN_SUCCESS
	
	def UpdateCFDLevel1 (self, arguments) :
		"""
		Interactively update
		"""
		if self.lib.PH_SetCFDLevel (self.DEV_INDX, ctypes.c_int(1), arguments[0]) < 0 :
			print "Interactive update of PH_SetCFDLevel 1 failed") 
			return RETURN_FAIL
		else :
			return RETURN_SUCCESS

	def UpdateCFDZero0 (self, arguments) :
		"""
		Interactively update
		"""
		if self.lib.PH_SetCFDZeroCross (self.DEV_INDX, ctypes.c_int(0), arguments[0]) < 0 :
			print "Interactive update of PH_SetCFDZeroCross 0 failed" 
			return RETURN_FAIL
		else :
			return RETURN_SUCCESS

	def UpdateCFDZero1 (self, arguments) :
		"""
		Interactively update
		"""
		if self.lib.PH_SetCFDZeroCross (self.DEV_INDX, ctypes.c_int(1), arguments[0]) < 0 :
			print "Interactive update of PH_SetCFDZeroCross 1 failed"
			return RETURN_FAIL
		else :
			return RETURN_SUCCESS
			
	def UpdateSyncDiv (self, arguments) :
		"""
		Interactively update
		"""
		if self.lib.PH_SetSyncDiv (self.DEV_INDX, arguments[0]) < 0 :
			print "Interactive update of PH_SetSyncDiv failed"
			return RETURN_FAIL
		else :
			# Note: after Init or SetSyncDiv you must allow 100 ms for valid new count rate readings
			time.sleep (0.1)
			return RETURN_SUCCESS
			
	def UpdateTacq (self, arguments) :
		"""
		Interactively update
		"""
		self.tacq = arguments[0]
		return RETURN_SUCCESS
		
	def UpdateRange (self, arguments) :
		"""
		Interactively update
		"""
		if self.lib.PH_SetRange (self.DEV_INDX, arguments[0]) < 0 :
			print "Interactive update of PH_SetRange failed"
			return RETURN_FAIL
		else :
			return RETURN_SUCCESS
			
########################################################################

class PicoHarpTab (HardwareGUIControl) :
	"""
	This class represents a GUI controlling properties of PicoHarp.
	"""
	def __init__(self, parent, dev) :
		HardwareGUIControl.__init__(self, parent, dev)
		
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		# Offset
		sizer.Add (wx.StaticText(self, label="Offset (ns)"), flag=wx.LEFT|wx.TOP, border=5)
		offset_ctrl = wx.SpinCtrl (self, value="0", min=0, max=1e6)
		offset_ctrl.__label__ = "offset"
		sizer.Add (offset_ctrl, flag=wx.EXPAND, border=5)
		# Enable interactive update 
		offset_ctrl.Bind (wx.EVT_SPINCTRL, 
			functools.partial(self.UpdateProperty, offset_ctrl, "UpdateOffset") 
		)
		
		# Zero level of channel 0
		sizer.Add (wx.StaticText(self, label="CFDZeroX0 (mV)"), flag=wx.LEFT, border=5)
		CFDZeroX0_ctrl = wx.SpinCtrl (self, value="10", min=0, max=1e3)
		CFDZeroX0_ctrl.__label__ = "cfdzerox0"
		sizer.Add (CFDZeroX0_ctrl, flag=wx.EXPAND, border=5)
		# Enable interactive update 
		CFDZeroX0_ctrl.Bind (wx.EVT_SPINCTRL, 
			functools.partial(self.UpdateProperty, CFDZeroX0_ctrl, "UpdateCFDZero0") 
		)
		
		# One level for channel 0
		sizer.Add (wx.StaticText(self, label="CFDLevel0 (mV)"), flag=wx.LEFT, border=5)
		CFDLevel0_ctrl = wx.SpinCtrl (self, value="80", min=0, max=1e3)
		CFDLevel0_ctrl.__label__ = "cfdlevelx0"
		sizer.Add (CFDLevel0_ctrl, flag=wx.EXPAND, border=5)
		# Enable interactive update 
		CFDLevel0_ctrl.Bind (wx.EVT_SPINCTRL, 
			functools.partial(self.UpdateProperty, CFDLevel0_ctrl, "UpdateCFDLevel0") 
		)
		
		# Zero level of channel 1
		sizer.Add (wx.StaticText(self, label="CFDZeroX1 (mV)"), flag=wx.LEFT, border=5)
		CFDZeroX1_ctrl = wx.SpinCtrl (self, value="10", min=0, max=1e3)
		CFDZeroX1_ctrl.__label__ = "cfdzerox1"
		sizer.Add (CFDZeroX1_ctrl, flag=wx.EXPAND, border=5)
		# Enable interactive update 
		CFDZeroX1_ctrl.Bind (wx.EVT_SPINCTRL, 
			functools.partial(self.UpdateProperty, CFDZeroX1_ctrl, "UpdateCFDZero1") 
		)
		
		# One level for channel 1
		sizer.Add (wx.StaticText(self, label="CFDLevel1 (mV)"), flag=wx.LEFT, border=5)
		CFDLevel1_ctrl = wx.SpinCtrl (self, value="80", min=0, max=1e3)
		CFDLevel1_ctrl.__label__ = "cfdlevelx1"
		sizer.Add (CFDLevel1_ctrl, flag=wx.EXPAND, border=5)
		# Enable interactive update 
		CFDLevel1_ctrl.Bind (wx.EVT_SPINCTRL, 
			functools.partial(self.UpdateProperty, CFDLevel1_ctrl, "UpdateCFDLevel1") 
		)
		
		# Divider
		sizer.Add (wx.StaticText(self, label="Divider"), flag=wx.LEFT, border=5)
		syncdiv_ctrl = wx.SpinCtrl (self, value="8", min=1, max=1e3)
		syncdiv_ctrl.__label__ = "syncdiv"
		sizer.Add (syncdiv_ctrl, flag=wx.EXPAND, border=5)
		# Enable interactive update 
		syncdiv_ctrl.Bind (wx.EVT_SPINCTRL, 
			functools.partial(self.UpdateProperty, syncdiv_ctrl, "UpdateSyncDiv") 
		)
		
		# Range	
		sizer.Add (wx.StaticText(self, label="Range"), flag=wx.LEFT, border=5)
		range_ctrl = wx.SpinCtrl (self, value="2", min=0, max=1e6)
		range_ctrl.__label__ = "range"
		sizer.Add (range_ctrl, flag=wx.EXPAND, border=5)
		# Enable interactive update 
		range_ctrl.Bind (wx.EVT_SPINCTRL, 
			functools.partial(self.UpdateProperty, range_ctrl, "UpdateRange") 
		)
		
		# Acquisition time
		sizer.Add (wx.StaticText(self, label="Acquisition time (ms)"), flag=wx.LEFT, border=5)
		tacq_ctrl = wx.SpinCtrl (self, value="1000", min=100, max=1e6)
		tacq_ctrl.__label__ = "tacq"
		sizer.Add (tacq_ctrl, flag=wx.EXPAND, border=5)
		# Enable interactive update 
		tacq_ctrl.Bind (wx.EVT_SPINCTRL, 
			functools.partial(self.UpdateProperty, tacq_ctrl, "UpdateTacq") 
		)
		
		# Spacer
		sizer.Add (wx.StaticText(self, label=""), flag=wx.LEFT, border=5)
		
		# Button for displaying histogram
		show_histogram_button = MultiStateButton(self, 
			actions=[self.StartShowingHistogram, self.StopShowingHistogram], 
			labels=["Show histogram", "STOP histogram"] 
		)
		sizer.Add (show_histogram_button, flag=wx.EXPAND, border=5)
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		
		self.CreateSettingsDict()
		
	def StartShowingHistogram (self, event) :
		pass
		
	def StopShowingHistogram (self, event) :
		pass