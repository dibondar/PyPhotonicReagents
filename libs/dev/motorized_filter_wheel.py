##########################################################################################
#
#	This file contains classes for controlling FW102C Thorlabs  Motorized Filter Wheel
#
##########################################################################################

import wx
import functools
import ctypes
import multiprocessing
 
from libs.gui.hardware_control import HardwareGUIControl
from libs.dev.basic_device import BasicDevice 
from libs.dev.consts import * 

########################################################################
#
#	Manager of the wheel
#
########################################################################

class ManagerFilterWheel :
	"""
	Class that manges W102C Thorlabs  Motorized Filter Wheel
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
		Start the process controlling the filter wheel
		"""
		p = FilterWheel(self.child_connection)
		p.start()
		return p
	
	def exit(self) :
		"""
		Close the process
		"""
		self.StopDevice()
		return self.Exit()
	
	def run(self, command, *args) :
		"""
		Send the command to the filter wheel through the pipe
		"""
		if len(args) == 0 : args = None
		elif len(args) == 1 : args = args[0]
		
		self.lock.acquire()
		self.parent_connection.send( (command, args) )
		result = self.parent_connection.recv()
		self.lock.release()
		return result
	
	def __getattr__ (self, name) :
		"""
		Redirect all other request to the method `run`
		"""
		return functools.partial( self.run, name )

########################################################################
#
#	Process where motorized filter wheel resides 
#
########################################################################

class FilterWheel (BasicDevice):
	"""
	Control motorized filter wheel
	"""
	def Initialize (self, settings) :
	
		# Opening the driver
		driver_path = self.GetProgramFilesDirectory() + "\\Thorlabs\\FW102C\\Sample\\msvc\\uart_library.dll"
		self.lib = ctypes.cdll.LoadLibrary (driver_path)

		if self.lib.fnUART_LIBRARY_open (settings["port"], 115200) != 0 :
			raise ValueError( "FilterWheel Error: device is not initialize" )
		
		return RETURN_SUCCESS
		
	def StopDevice (self, arg=None) :
		"""
		Close the device
		"""
		try :
			self.lib.fnUART_LIBRARY_close()
			del self.lib
		except AttributeError : pass
		
		return RETURN_SUCCESS
	
	def GetNumFilters (self, arg=None) :
		"""
		Return maximum number of filters used
		"""
		return 6
		
	def SetFilter (self, filter_num) :	
		"""
		Go to filter specified by its number
		"""
		assert ( 1 <= filter_num ) and ( filter_num <= self.GetNumFilters() )
		
		if self.lib.fnUART_LIBRARY_Set( ctypes.c_char_p("pos=%d\r" % filter_num ), 0 ) != 0 :
			raise RuntimeError("FilterWheel Error: cannot go to filter %n" %  filter_num )
			
		return RETURN_SUCCESS
		
########################################################################

class FilterWheelTab (HardwareGUIControl) :
	"""
	This class represents a GUI controlling properties of the filter wheel.
	"""
	def __init__(self, parent, dev) :
		HardwareGUIControl.__init__(self, parent, dev)
		
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		################################################
		# Specify the communication port name
		sizer.Add (wx.StaticText(self, label="Communication port"), flag=wx.LEFT, border=5)
		port_name = wx.SpinCtrl (self, value="18")
		port_name.__label__ = "port"
		sizer.Add (port_name, flag=wx.EXPAND, border=5)
		
		################## Select filter ###################
		sizer.Add( wx.StaticText(self, label="\nSelect filter"), flag=wx.LEFT, border=5)
		self.filter_num_ctrl = wx.SpinCtrl (self, min=1, max=6, value="1")
		self.filter_num_ctrl.__label__ = "selected_filter"
		
		def GoToFilter (event) :
			if self.dev.Initialize( self.GetSettings() ) == RETURN_FAIL : return
			self.dev.SetFilter( self.filter_num_ctrl.GetValue() )
		
		self.filter_num_ctrl.Bind(wx.EVT_SPINCTRL, GoToFilter)
		sizer.Add (self.filter_num_ctrl, flag=wx.EXPAND, border=5)
		################################################
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		
		self.CreateSettingsDict()
	
		