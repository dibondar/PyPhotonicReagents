"""
This module contains classes for controlling and GUI representation of
Thorlabs cameras (http://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=4024)
"""

from libs.gui.hardware_control import HardwareGUIControl
from libs.dev.basic_manager import BasicManager
from libs.dev.basic_device import BasicDevice

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

class ManagerThorlabsCamera :
	"""
	Class that manges Thorlabs cameras 
	"""	
	def start(self) :
		"""
		Start the process controlling thorlabs cameras
		"""
		p = ThorlabsCamera(self.child_connection, self.histogram_buffer)
		p.start()
		return p
		
###########################################################################