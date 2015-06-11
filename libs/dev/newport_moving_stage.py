"""
This module contains classes for controlling and GUI representation of
Newport moving stage (http://www.newport.com/MFA-Series-Miniature-Steel-Linear-Stages/300762/1033/info.aspx)
"""

from libs.gui.hardware_control import HardwareGUIControl
from libs.dev.basic_device import BasicDevice
from libs.dev.basic_manager import BasicManager

import wx
import serial
import functools
import multiprocessing
from libs.dev.consts import * 

########################################################################
#
#	Manager class that communicates with the process where
#	the device resides
#
########################################################################

class ManagerNewportMovingStage (BasicManager):
	"""
	Class that manges Newport moving stage 
	"""
	def start(self) :
		"""
		Start the process controlling Newport moving stage 
		"""
		p = NewportMovingStage(self.child_connection, self.histogram_buffer)
		p.start()
		return p
		
###########################################################################


class NewportMovingStage (BasicDevice):
	"""
	Control spectrometer hardware 
	"""
	def Initialize (self, settings) :
		# Close the port if it is already used
		try : del self.serial_port
		except AttributeError : pass 
	
		# Start the communication port
		self.serial_port = serial.Serial (port=settings["port"], 
			baudrate=19200,  bytesize=8, parity=serial.PARITY_NONE, stopbits=1, timeout=1.)
		
		return RETURN_SUCCESS
		
	def GetCurrentPosition (self, arguments=None) :
		"""
		Return current position of the moving stage
		"""
		# issue command
		self.serial_port.write("1TP?;2TP?;3TP?\r")

		# read the response of the moving stage
		response = self.serial_port.readlines()
		
		try :
			return tuple( float(_) for _ in response ) 
		except ValueError :
			print "Warning: Moving stage response %s count not be converted to int" % str(response) 
			return RETURN_FAIL
		
	def IsMoving (self, arguments=None) :
		"""
		Check whether motion is completed
		"""
		# issue command
		self.serial_port.write("1MD?;2MD?;3MD?\r")
		
		# read the response of the moving stage
		response = self.serial_port.readlines()
		
		try :
			return all( int(_) for _ in response ) 
		except ValueError :
			return False 
		
	def WaitUntillStops (self, arguments=None) :
		"""
		Wait till motion is over 
		"""
		# issue command
		self.serial_port.write("1WS0;2WS0;3WS0;1MD?;2MD?;3MD?\r")
		
		# read the response of the moving stage
		response = self.serial_port.readlines()
		
		try :
			if not all( int(_) for _ in response ) :
				print "Warning: Moving stage may not finished moving"
		except ValueError :
			print "Warning: Moving stage response %s count not be converted to int" % str(response) 
		
		return RETURN_SUCCESS

	def MoveAbs (self, position) :
		"""
		Moving to new absolute position
		"""
		self.serial_port.write("1PA%.8e;2PA%.8e;3PA%.8e\r" % position)
		return RETURN_SUCCESS
	
	def MoveAbsWait (self, position) :
		"""
		Moving to new absolute position and wait till motion stops
		"""
		self.MoveAbs(position)
		return self.WaitUntillStops()
		
	def MoveRel (self, position) :
		"""
		Moving to new relative position
		"""
		self.serial_port.write("1PR%.8e;2PR%.8e;3PR%.8e;\r" % position)
		return RETURN_SUCCESS
	
	def MoveRelWait (self, position) :
		"""
		Moving to new relative position and wait till motion stops
		"""
		self.MoveRel(position)
		return self.WaitUntillStops()
		
###########################################################################
			
class NewportMovingStageTab (HardwareGUIControl) :
	"""
	This class represents a GUI controlling properties of Newport moving stage.
	"""
	def __init__(self, parent, dev=None) :
		HardwareGUIControl.__init__(self, parent, dev)
		
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		# Specify the communication port name
		sizer.Add (wx.StaticText(self, label="Communication port"), flag=wx.LEFT, border=5)
		port_name = wx.TextCtrl (self, value="COM15")
		port_name.__label__ = "port"
		sizer.Add (port_name, flag=wx.EXPAND, border=5)
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		
		self.CreateSettingsDict()