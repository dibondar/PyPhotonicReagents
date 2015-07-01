########################################################################
#
#	This file contains the abstract base class for device as well as 
#	GUI representation of device settings 
#
########################################################################

from libs.dev.consts import * 
from h5py._hl.dataset import Dataset as HD5Dataset
import wx
from wx.lib.agw.floatspin import FloatSpin as wxFloatSpin
import multiprocessing 
import numpy as np


class HardwareGUIControl (wx.Panel) :
	"""
	This abstract class represents a GUI controlling properties of a device.
	Note that the children of this class must call method <self.CreateSettingsDict()>
	after all GUI elements were added.
	"""
	def __init__(self, parent=None, dev=None, manager_cls=None, **kwargs) :
		"""
		`parent` is a wx.Notebook container containing all other controls
		`dev`	is a device whose settings are contained in the current control
		`manager_cls` is a class that can initialize `dev`
		"""
		wx.Panel.__init__(self, parent, **kwargs)
		self.parent = parent
		self.dev = dev
		self.manager_cls = manager_cls
	
	def StartDev (self) :
		"""
		This function implements the singleton patern of the device initialization 
		"""
		if self.dev is None :
			# Device has not been initialized yet
			if self.manager_cls :
				# Use the specified device manager to initialize the class
				self.dev = self.manager_cls()
				self.dev.start()
		
		self.UpdateSettings()
		
		return self.dev
		
	def StopDev (self) :
		"""
		Stop the device
		"""
		if self.dev is not None :
			self.dev.Exit()
			self.dev = None
	
	def UpdateProperty (self, control, update_func_name, *args, **kwords) :
		"""
		This method is used to intransitively update (via `control`) a particular settings  
		of the device by calling the function `update_func_name` in the device manager
		"""
		if self.dev is not None :
			getattr(self.dev, update_func_name)( control.GetValue() )
	
	def UpdateSettings (self, *args, **kwords) :
		"""
		This method is used to interactively update all hardware settings, if the device is initialized 
		"""
		if self.dev is not None :
			if self.dev.SetSettings( self.GetSettings() ) == RETURN_FAIL : 
				raise RuntimeError("Settings could not be set")
				
	def CreateSettingsDict (self) :
		"""
		Generate a dictionary with settings that bind to control elements
		"""
		self.settings_to_controls = {}
		for control in self.GetChildren() :
			if isinstance(control, (wx.SpinCtrl, wx.TextCtrl, wx.RadioButton, 
				wx.CheckBox,  wx.FilePickerCtrl, wxFloatSpin, wx.ComboBox) ) :
				# Extracting the property name
				try :  label = control.__label__ 
				except AttributeError : label = control.GetLabel()
				
				# Adapt <wx.FilePickerCtrl> to desired specification
				if isinstance(control, wx.FilePickerCtrl) :
					control.GetValue = control.GetPath; control.SetValue = control.SetPath
				
				# Generating the property name based on the label of control
				setting_key = label.replace(' ', '_').lower()
				
				# Checking for duplicates
				if setting_key in self.settings_to_controls : 
					raise ValueError ("Setting <%s> has already been added" % setting_key)
				
				# Saving binding between the setting and controls
				self.settings_to_controls[ setting_key ] = control
				
				
	def GetSettings(self) :
		"""
		Save the settings into a dictionary 
		"""
		return dict( (key, control.GetValue()) for key, control in self.settings_to_controls.items() )
	
	def SetSettings(self, settings) :
		"""
		Assign settings. This method is similar to <self.GetSettings> 
		"""
		for setting_key, value in settings.items() : 
			try : 
				control = self.settings_to_controls[setting_key]
				
				# Loading the HDF5 dataset 
				if type(value) is HD5Dataset : value=value[...]
				
				if isinstance(control, (wx.TextCtrl,wx.FilePickerCtrl,wx.ComboBox) ) : 
					control.SetValue(str(value))
				else : control.SetValue(value)
				
			except KeyError :
				print "SetSettings Warning: Parameter %s is not recognized and hence ignored" % setting_key