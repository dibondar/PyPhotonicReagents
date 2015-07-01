"""
This module contains classes for controlling and GUI representation of
Newport moving stage (http://www.newport.com/MFA-Series-Miniature-Steel-Linear-Stages/300762/1033/info.aspx)
"""

from libs.gui.hardware_control import HardwareGUIControl
from libs.gui.multi_state_button import MultiStateButton
from libs.dev.basic_device import BasicDevice
from libs.dev.basic_manager import BasicManager

import wx
import serial
import functools
import numpy as np
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
		p = NewportMovingStage(self.child_connection)
		p.start()
		return p
		
###########################################################################


class NewportMovingStage (BasicDevice):
	"""
	Control Newport moving stage
	"""
	def SetSettings (self, settings) :
		# Close the port if it is already used
		try : 
			del self.serial_port
		except AttributeError : 
			pass 
	
		print "cool"
		# Start the communication port
		self.serial_port = serial.Serial (port=settings["port"], 
			baudrate=19200,  bytesize=8, parity=serial.PARITY_NONE, stopbits=1, timeout=1.
		)
		
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

class JoystickBinding (multiprocessing.Process) :
	"""
	Class realizing binding of joystick to moving stage 
	"""
	def __init__ (self, MovingStage, joystick_event, *args, **kwargs) :
		"""
		MovingStage -- manager of moving stage
		joystick_event - multiprocessing.Event to indicate when the joystic is on
		"""
		self.MovingStage = MovingStage
		self.joystick_event = joystick_event
		multiprocessing.Process.__init__(self, *args, **kwargs)
	
	def Initialize (self) :
		"""
		Initialize joystick via pygame
		"""
		# Initializing joystick
		import pygame
		self.pygame = pygame
		self.pygame.init()
		self.pygame.joystick.init()
		
		if self.pygame.joystick.get_count() == 0 or not self.pygame.joystick.get_init() :
			raise RuntimeError("No joystick is detected!")
		
		if self.pygame.joystick.get_count() <> 1 :
			print "Multiple joysticks are detected: the first one will be utilized"
	
		self.joystick = pygame.joystick.Joystick(0)
		self.joystick.init ()
		print 'Joystick "%s" was initialized' % self.joystick.get_name() 
		
	def UpdateMovingStagePos (self) :
		"""
		Update poisonings of moving stages based on the positions of joystick
		"""
		# Values of the joystick control so that motion is not detected
		motion_cut_off = 0.1
		
		# reading off the position of the joystick 
		# Note that the labels are transposed to align the joystick with the moving stage
		possition = np.array( [ self.joystick.get_axis(n) for n in [0,1,3] ] )
		
		# Find which axis are below cut off
		indx = ( np.abs(possition) < motion_cut_off )
		
		if indx.sum() < indx.size  :  
			# User moved joystick
				
			# Check whether the moving stage does not move
			if self.MovingStage.IsMoving() :
				
				# The rudder control is used to scale the step size
				scale = 0.05*(1. - self.joystick.get_axis(2))
				possition *= scale
					
				# Do not change values of axis if it is bellow cut off
				possition[indx] = 0
			
				# If joystick trigger button is on, then move z-axis only
				if self.joystick.get_button(0) :
					possition[:2] = 0
					possition[2] *= 0.1
				else : 
					# otherwise move x-y only
					possition[2] = 0
			
				# The moving stage is ready to be moved
				self.MovingStage.MoveRel( tuple(possition) )
		
			# Clearing the event stack
			self.pygame.event.clear ()	
	
	def  Stop (self) :
		self.joystick.quit ()
		self.pygame.joystick.quit()
		self.pygame.quit ()
	
	def run (self) :
		"""
		Overloaded function provided by multiprocessing.Process. Called upon start() signal 
		"""
		# Initialize joystick
		self.Initialize()

		# Continue linking joystick until user quits
		while self.joystick_event.is_set() :
			self.UpdateMovingStagePos()
		
		# Releasing resources
		self.Stop()		
		
###########################################################################

class NewportMovingStageTab (HardwareGUIControl) :
	"""
	This class represents a GUI controlling properties of Newport moving stage
	"""
	def __init__(self, parent, **kwards) :
		HardwareGUIControl.__init__(self, parent, manager_cls=ManagerNewportMovingStage, **kwards)
		
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		# Specify the communication port name
		sizer.Add (wx.StaticText(self, label="\nCommunication port"), flag=wx.LEFT, border=5)
		port_name = wx.TextCtrl (self, value="COM15")
		port_name.__label__ = "port"
		sizer.Add (port_name, flag=wx.EXPAND, border=5)
		
		# Spacer
		sizer.Add (wx.StaticText(self, label=""), flag=wx.LEFT, border=5)
		
		# Button to start joystick
		bind_joystick_button = MultiStateButton(self, states=[ 
				(self.StartJoystick, "Start joystick"), (self.StopJoystick, "STOP joystick")
			]
		)
		sizer.Add (bind_joystick_button, flag=wx.EXPAND, border=5)
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		
		self.CreateSettingsDict()
		
	def StartJoystick (self, event) :
		"""
		Bind joystick to moving stage
		"""
		# Decide whether the device needs to be destroy at the end of usage
		self._stop_dev = ( self.dev is None )
		
		self.StartDev()
		
		# Start joystick binding
		self.joystick_event = multiprocessing.Event()
		self.joystick_event.set ()
		JoystickBinding (self.dev, self.joystick_event).start()
		
	def StopJoystick (self, event) :
		"""
		Stop joystick-moving stage binding
		"""
		# Send the signal to stop synchronizing joystick
		self.joystick_event.clear()
	
		if self._stop_dev :
			self.StopDev()
		