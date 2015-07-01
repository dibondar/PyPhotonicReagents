"""
This module contains classes for controlling and GUI representation of
Thorlabs cameras (http://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=4024)
"""

from libs.gui.hardware_control import HardwareGUIControl
from libs.gui.multi_state_button import MultiStateButton

from libs.dev.basic_manager import BasicManager
from libs.dev.basic_device import BasicDevice

from libs.dev.consts import * 

import wx
import ctypes
import serial
import visvis
import functools
import multiprocessing
import numpy as np

########################################################################
#
#	Manager class that communicates with the process where
#	the device resides
#
########################################################################

class ManagerThorlabsCamera (BasicManager):
	"""
	Class that manges Thorlabs cameras 
	"""	
	def start(self) :
		"""
		Start the process controlling thorlabs cameras
		"""
		self.camera_img_buffer = ThorlabsCamera.AllocateBuffer()
		p = ThorlabsCamera(self.child_connection, self.camera_img_buffer)
		p.start()
		# Initialize device
		self.Initialize()
		return p
	
	def StartImgAcquisition (self) :
		"""
		Start acquiring image in a parallel fashion
		"""
		self.lock.acquire()
		self.parent_connection.send( ("GetImage", None) )
	
	def StopImgAcquisition (self) :
		"""
		Wait till the image acquisition finished and return a reference  
		"""
		result = self.parent_connection.recv()
		self.lock.release()
		if result == RETURN_FAIL : 
			print "Image acquisition failed"
		return  ThorlabsCamera.ImgBuffer2Numpy(self.camera_img_buffer)
		
	def GetImage (self) :
		"""
		Acquire histogram sequentially
		"""
		self.StartImgAcquisition()
		return self.StopImgAcquisition()
		
###########################################################################
# constants

CAMERA_IMG_WIDTH = 1280
CAMERA_IMG_HIGHT = 1024

IS_SUCCESS 			= 0
IS_CAPTURE_RUNNING 	= 140	

###########################################################################

class ThorlabsCamera (BasicDevice):
	"""
	Control Thorlabs camera
	"""
	@classmethod
	def AllocateBuffer (cls) :
		return multiprocessing.Array (ctypes.c_uint8, CAMERA_IMG_WIDTH*CAMERA_IMG_HIGHT)

	@classmethod
	def ImgBuffer2Numpy (cls, camera_img_buffer) :
		"""
		Convert image buffer allocated by cls.AllocateBuffer into a numpy array
		"""
		img_buffer = np.frombuffer(camera_img_buffer.get_obj(), dtype=np.uint8)
		img_buffer = img_buffer.reshape( (CAMERA_IMG_HIGHT, CAMERA_IMG_WIDTH) )
		return img_buffer
		
	def __init__ (self,  pipe, camera_img_buffer) :
		"""
		<pipe> is for communicating.
		<histogram_buffer> is a multiprocessing.Array where histogram is saved.
		"""
		BasicDevice.__init__(self, pipe)
		# saving the buffer where images will be saved
		self.camera_img_buffer 	= camera_img_buffer
		self.img_buffer			= self.ImgBuffer2Numpy(self.camera_img_buffer)
		
	def Initialize (self, aguments=None) :
		"""
		Initialize Thorlabs camera
		"""
		# Load driver
		self.lib = ctypes.cdll.LoadLibrary("uc480")
		
		num_cams = ctypes.c_int()
		if self.lib.is_GetNumberOfCameras( ctypes.byref(num_cams) ) != IS_SUCCESS :
			raise RuntimeError("Error in is_GetNumberOfCameras")

		if num_cams.value > 1 :
			print "More than one camera is detect. The first one will be used."
		
		# Camera handle
		self.ph_cam = ctypes.wintypes.HDC(0)

		# Window handle
		h_wnd = ctypes.wintypes.HWND(0)

		# Initialize camera
		if self.lib.is_InitCamera( ctypes.byref(self.ph_cam), h_wnd ) != IS_SUCCESS :
			raise RuntimeError("Error in is_InitCamera")
			
		# ID of the memory for image storage
		pid = ctypes.c_int()
		# pointer to buffer
		pc_img_mem = self.img_buffer.ctypes.data_as(ctypes.POINTER(ctypes.c_char))
		
		if self.lib.is_SetAllocatedImageMem (self.ph_cam,  CAMERA_IMG_WIDTH, CAMERA_IMG_HIGHT, 
			8*self.img_buffer.dtype.itemsize, pc_img_mem, ctypes.byref(pid) 
		) != IS_SUCCESS :
			raise RuntimeError("Error in is_SetAllocatedImageMem")

		# Set colour
		IS_CM_MONO8 = 6
		if self.lib.is_SetColorMode (self.ph_cam,  IS_CM_MONO8) != IS_SUCCESS :
			raise RuntimeError("Error in is_SetColorMode")
			
		# Set active memory
		if self.lib.is_SetImageMem(self.ph_cam, pc_img_mem, pid) != IS_SUCCESS :
			raise RuntimeError("Error in is_SetImageMem")	
			
	def SetSettings(self, settings) :
		"""
		Set settings of Thorlabs camera
		"""
		# pre-calculate necessary information, if user choose to draw a gun-sight on the image
		self.draw_gun_sight = settings["draw_gun_sight"]
		if self.draw_gun_sight :
			# Draw circule
			sqr = (np.arange(self.img_buffer.shape[0])[:,np.newaxis] - self.img_buffer.shape[0]/2.)**2 \
				+ (np.arange(self.img_buffer.shape[1])[np.newaxis,:] - self.img_buffer.shape[1]/2.)**2
			radius = settings["radius"]
			self.circular_cut = np.nonzero( ( (radius-1)**2 < sqr)&(sqr < (radius+1)**2 )  )
	
	def GetImage(self, arguments=None) :
		"""
		Acquire image
		"""
		self.camera_img_buffer.acquire ()
		
		if self.lib.is_FreezeVideo(self.ph_cam, 0x0001) in [IS_SUCCESS, IS_CAPTURE_RUNNING] :
			if self.draw_gun_sight :
				# Drawing the "gun" sight
				self.img_buffer[self.img_buffer.shape[0]/2, :] = 0
				self.img_buffer[:, self.img_buffer.shape[1]/2] = 0
				# Drawing the circule 
				self.img_buffer[ self.circular_cut ] = 0	
			
			result = RETURN_SUCCESS
		else :
			result = RETURN_FAIL
		
		self.camera_img_buffer.release ()
		
		return result
		
###########################################################################

class ThorlabsCameraTab (HardwareGUIControl) :
	"""
	This class represents a GUI controlling properties of Thorlabs camera
	"""
	def __init__(self, parent, **kwards) :
		HardwareGUIControl.__init__(self, parent, manager_cls=ManagerThorlabsCamera, **kwards)
	 
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		# Spacer
		sizer.Add (wx.StaticText(self, label=""), flag=wx.LEFT, border=5)
		
		# Check-box to draw gun-sight
		draw_gun_sight_ctl = wx.CheckBox(self, label="Draw gun-sight")
		draw_gun_sight_ctl.__label__ = "draw_gun_sight"
		sizer.Add (draw_gun_sight_ctl, flag=wx.EXPAND, border=5)
		# Enable interactive update 
		draw_gun_sight_ctl.Bind(wx.EVT_CHECKBOX, self.UpdateSettings)
		
		# Radius of gun-sight
		sizer.Add (wx.StaticText(self, label="\nGun-sight radius"), flag=wx.LEFT, border=5)
		gun_sight_radius_ctrl = wx.SpinCtrl (self, value="43", min=0, 
											max=min(CAMERA_IMG_WIDTH, CAMERA_IMG_HIGHT) )
		gun_sight_radius_ctrl.__label__ = "radius"
		sizer.Add (gun_sight_radius_ctrl, flag=wx.EXPAND, border=5)
		# Enable interactive update 
		gun_sight_radius_ctrl.Bind (wx.EVT_SPINCTRL, self.UpdateSettings)
		
		# Spacer
		sizer.Add (wx.StaticText(self, label=""), flag=wx.LEFT, border=5)
		
		# Button for displaying image 
		show_img_button = MultiStateButton(self, states=[
				(self.StartShowingImg, "Show image"), (self.StopShowingImg, "STOP")
			]
		)
		sizer.Add (show_img_button, flag=wx.EXPAND, border=5)
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		
		self.CreateSettingsDict()
		
	def StartShowingImg (self, event) :
		"""
		Start continuous acquisition of image
		"""
		# Decide whether the device needs to be destroy at the end of usage
		self._stop_dev = ( self.dev is None )
		self.StartDev()
		
		# Start acquisition of image
		self.dev.StartImgAqusition()
	
		# Start displaying images
		self._abort_img = False
		wx.CallAfter(self.ShowImg)
		
	def StopShowingImg (self, event) :
		"""
		Stop showing images interactively
		"""
		self._abort_img = True
		self.dev.StopImgAqusition()
		
		if self._stop_dev :
			self.StopDev()
		
		del self._img_plot
		
	def ShowImg (self) :
		"""
		Draw image
		"""
		if self._abort_img  :
			# Exit
			return
		
		# Get image
		img = self.dev.StopImgAqusition()
		
		# Display image
		try :
			if img <> RETURN_FAIL :
				self._img_plot.SetData(img)
		except AttributeError :
			visvis.cla()
			visvis.clf()
			self._img_plot = visvis.imshow(img)
			visvis.title ('Camera view')
			ax = visvis.gca()
			ax.axis.xTicks = []
			ax.axis.yTicks = []
		
		# Start acquisition of histogram
		self.dev.StartImgAqusition()
		wx.CallAfter(self.ShowImg)