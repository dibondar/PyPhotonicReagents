########################################################################################
#
#	Base class for hardware manager
#
########################################################################################

import functools
import multiprocessing 

class BasicManager :
	"""
	Parent class for device managers
	"""
	def __init__ (self) :
		# Create the lock for device 
		self.lock = multiprocessing.Lock()
		# Create a pipe for communication
		self.parent_connection, self.child_connection = multiprocessing.Pipe()
		
	def __del__ (self) :
		self.parent_connection.close()
		self.child_connection.close()
		
	def run(self, command, arguments=None) :
		"""
		Send the command to the device through the pipe
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