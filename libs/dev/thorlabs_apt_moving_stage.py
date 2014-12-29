########################################################################
#
#	This module contains classes for controlling and displying properties 
#	of the APT Thorlabs moving stage
#
########################################################################

import multiprocessing
from consts import * 


########################################################################
#
#	Managing class that communicates with the process where
#	the device resides
#
########################################################################

class ManagerThorlabsAPTMovingStage :
	"""
	Class that manges the moving stage
	"""
	def __init__ (self, SerialNumber) :
		# Create the lock for device 
		self.lock = multiprocessing.Lock()
		# Create a pipe for communication
		self.parent_connection, self.child_connection = multiprocessing.Pipe()
		# Saving the serial number
		self.SerialNumber = SerialNumber
		
	def __del__ (self) :
		self.parent_connection.close()
		self.child_connection.close()
		
	def start(self) :
		"""
		Start the process controlling the moving stage
		"""
		p = multiprocessing.Process(target=MovingStageProc, args=(self.SerialNumber,self.child_connection))
		p.start()
		return p
		
	def run(self, command, arguments=None) :
		"""
		Send the command to a moving stage through the pipe
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
		
	def AbsMove(self, position) :
		"""
		Move the moving stage without waiting the completion of the motion.
		"""
		return self.run( "AbsMove", position )
		
	def AbsMoveWait(self, position) :
		"""
		Move the moving stage and wait till completion
		"""
		return self.run( "AbsMoveWait", position )

########################################################################
#
#	Process where the device resides
#
########################################################################

def MovingStageProc(SerialNumber, pipe) :
	"""
	This function will be run as a separate process
	"""

	import wx.lib.activex, threading
	
	class ThorlabsAPTMovingStage (wx.lib.activex.ActiveXCtrl) :
		"""
		Control moving stage
		"""
		def __init__ (self, SerialNumber, pipe) :	
			self.frame = wx.Frame (None, title="Moving stage SN: %d" % SerialNumber)
			panel = wx.Panel (self.frame)
			wx.lib.activex.ActiveXCtrl.__init__ (self, panel, 'MGMOTOR.MGMotorCtrl.1', size=self.frame.GetSizeTuple(), name='Moving stage')
			# Initializing the device by specifying its serial number 
			self.ctrl.StartCtrl()
			self.ctrl.HWSerialNum = SerialNumber
			self.ctrl.Identify()
			# Building simple gui
			sizer = wx.BoxSizer ()
			sizer.Add (self, flag=wx.EXPAND)
			panel.SetSizer (sizer)
			self.frame.Show()
			#self.frame.Hide()
			# saving the pipe
			self.pipe = pipe
			# starting thread for checking commands sent from other processes
			threading.Thread(target=self.CheckCommands).start()
		
		def AbsMove (self, possition) :
			self.ctrl.SetAbsMovePos(0, possition)
			self.ctrl.MoveAbsolute(0,0)
			return RETURN_SUCCESS
			
		def AbsMoveWait (self, position) :
			self.AbsMove(position)
			# Wait till movement finishes
			status = 0xFF
			while ((status >> 5)&0x1)|((status >> 6)&0x1) : status = self.ctrl.GetStatusBits_Bits(0)
			return RETURN_SUCCESS
			
		def CheckCommands (self) :
			# Checking the pipe for requested commands	
			# until "Exit" is sent
			for command, arguments in iter(self.pipe.recv, ("Exit",None)) :
				# Run requested command
				try : 
					try : result = getattr(self, command)(arguments)		
					except RuntimeError, e :
						result = RETURN_FAIL; print e.message;
				except AttributeError :
					print "\n\nPulse Shaper Error : Unrecognised command: " + command
					result = RETURN_FAIL
				
				# returning the result
				self.pipe.send(result)
				
			# Closing the process	
			wx.CallAfter(self.frame.Close)
			self.pipe.send(RETURN_SUCCESS)
			
	# Starting the wx application
	wx_app = wx.PySimpleApp(False)
	shutter = ThorlabsAPTMovingStage(SerialNumber, pipe)
	wx_app.MainLoop()
