"""
Implementation of SettingsNotebook 
"""
import sys
import wx
import h5py
import hashlib
import tempfile

class SettingsNotebook (wx.Notebook) :
	"""
	GUI container (wx.Notebook) for 
	"""
	def __init__ (self, parent, specs=[], auto_load=True, **kwargs) :
		"""
		specs -- list of tuple specifying tab settings, e..g,
					[ (tab_class, label(string), dev_name(string)) ] 
				if dev_name = None then the specified tuple does not represent device 
		"""
		wx.Notebook.__init__(self, parent, **kwargs)
		
		# dictionary of all tabs
		self._tabs = {}
		
		# dictionary of tabs representing devices
		self._devs = {}
		
		for tab_specs in specs :
		
			# Unpacking specifications
			if len(tab_specs) <= 1 :
				raise ValueError("Minimal specification of a tab is (tab_class, tab_label,)")
				
			tab 	= tab_specs[0]
			label	= tab_specs[1]
			try :
				dev_name = tab_specs[2]
			except IndexError :
				dev_name = label

			# Initialize the GUI Tab
			tab = tab(self)
			
			# Add to the notebook
			self.AddPage(tab, label)
			
			# Add to dict of settings
			self._tabs[label] = tab
			
			# Add to dict of devices, if applicable
			if dev_name :
				self._devs[dev_name] = tab
	
		# Load settings 
		if auto_load :
			try :
				self.AutoLoad()
				print("Settings are auto-loaded")
			except IOError :
				pass
			
	def __del__ (self) :
		self.StopAllDevices()
	
	def StartDev (self, dev_name) :
		"""
		Start device specified by name
		"""
		return self._devs[dev_name].StartDev()
		
	def StopAllDevices (self, event=None, auto_save=True) :
		"""
		Close all devices that have been initialized
		"""
		for tab in self._devs.values() :
			tab.StopDev()
		
		# Auto save settings
		if auto_save :
			self.AutoSave()
			print("Settings are auto-saved")
			
		if event :
			# Destroy parent window 
			event.GetEventObject().Destroy() 
		
	@classmethod
	def GetAutoSavingFilename (cls) :
		"""
		Return file name where settings are saved automatically
		"""
		return tempfile.gettempdir() + hashlib.sha224(sys.argv[0]).hexdigest() + '.hdf5'
	
	def AutoLoad (self) :
		"""
		Load settings automatically
		"""
		self.LoadSettings( filename = self.GetAutoSavingFilename() )
	
	def AutoSave (self) :
		"""
		Save settings automatically
		"""
		self.SaveSettings( filename = self.GetAutoSavingFilename() )
	
	def LoadSettings (self, event=None, filename="", title="Open HDF5 file to load settings") :
		"""
		Load settings. This method is closely related to <self.SaveSettings>
		Return filename where settings were saved
		"""
		if not len(filename) :
			# Ask user to select the file
			openFileDialog = wx.FileDialog(self, title, "", "",
                                   "HDF5 files (*.hdf5)|*.hdf5", 
								   wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR
								 )
			# Check whether user cancelled
			if openFileDialog.ShowModal() == wx.ID_CANCEL: 
				return None
		
			filename = openFileDialog.GetPath()
		
		with h5py.File (filename, 'r') as f_settings :
			for label, tab in f_settings["settings"].items() :
				try :
					self._tabs[label].SetSettings(tab)
				except KeyError :
					print "Load Settings Error: Settings %s are ignored" %  label
				
		return filename
		
	def SaveSettings (self, event=None, filename = "", 
						default_filename = "settings.hdf5", title="Open HDF5 file to save settings") :
		"""
		Method for saving setting 
		"""
		if not len(filename) :
			# Ask user to select the file
			openFileDialog = wx.FileDialog(self, title, "", default_filename, "HDF5 files (*.hdf5)|*.hdf5", 
							wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT | wx.FD_CHANGE_DIR)
			# Check whether user cancelled
			if openFileDialog.ShowModal() == wx.ID_CANCEL: 
				return None	
			filename = openFileDialog.GetPath()
			
		with h5py.File (filename, 'a') as f_settings :

			# Crete the group if it does not exist
			try : 
				parameters_grp = f_settings["settings"] 
			except KeyError : 
				parameters_grp = f_settings.create_group("settings")
				
			# Loop over all settings tab
			for label, tab in self._tabs.items() :
				# Save all settings on a given tab
				try : 
					del parameters_grp[label]
				except KeyError : pass
				grp = parameters_grp.create_group(label)
					
				for key, value in tab.GetSettings().items() : 
					grp[key] = value
			
		# return valid file name
		return filename
		
	def GetAllSettings(self) :
		"""
		Return a dictionary of all dictionary containing settings from all tabs 
		"""
		return dict(
			(label, tab.GetSettings()) for label, tab in self._tabs.items()
		)