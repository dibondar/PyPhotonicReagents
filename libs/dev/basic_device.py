########################################################################################
#
#	Base class for hardware
#
########################################################################################

import multiprocessing 
from libs.dev.consts import *  

class BasicDevice (multiprocessing.Process):
	"""
	Basic class representing the process where the device resides.
	Child classes should implement methods that realizes commands.
	"""
	def __init__ (self, pipe) :
		"""
		Save the pipe for communication
		"""
		multiprocessing.Process.__init__(self)
		# saving the pipe for communication with manager
		self.pipe = pipe
	
	@classmethod
	def GetProgramFilesDirectory (cls) :
		"""
		Return the program files directory for Windows systems
		"""
		import os
		try : return os.environ['PROGRAMFILES(X86)']	# for Win64
		except KeyError : return os.environ['PROGRAMFILES'] # for Win32

	def run (self) :
		"""
		Overloaded function provided by multiprocessing.Process.  Called upon start() signal
		"""	
		# Checking the pipe for requested commands	
		# until "Exit" is sent
		for command, arguments in iter(self.pipe.recv, ("Exit",None)) :
			# Find method matching the command 
			try : method = getattr(self, command)
			except AttributeError :
				print "\n\nPulse Shaper Error : Unrecognised command: " + command
				result = RETURN_FAIL
				continue
			
			# Method was found that matched the requested command 
			try : result = method(arguments) # Run requested command
			except RuntimeError, e :
				result = RETURN_FAIL; print e.message;
				
			# returning the result
			self.pipe.send(result)

		self.pipe.send(RETURN_SUCCESS)