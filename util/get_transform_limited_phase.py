"""
Find transform limited phase
"""
# Add main directory to enable imports  
if __name__ == '__main__' :
	import os
	os.sys.path.append(os.path.abspath('..'))

import wx, h5py, operator
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from itertools import izip_longest

# Real time plotting
import visvis

# GUI components
from libs.gui.basic_window import BasicWindow
from libs.gui.deap_algorthims_tabs import GATab, PixelWiseGA

# Hardware
#from libs.dev.spectrometer_ocean_optics import ManagerOceanOpticsSpectrometer as ManagerSpectrometer
#from libs.dev.spectrometer_ocean_optics import OceanOpticsSpectrometerTab as SpectrometerTab
from libs.dev.camera_istar import ManagerIStarCamera as ManagerSpectrometer
from libs.dev.camera_istar import IStarCameraTab as SpectrometerTab

from libs.dev.pulse_shaper import ManagerShaper, PulseShaperTab, PULSESHAPER_MAX_VAL, ShaperInt

from libs.dev.consts import * 

########################################################################

class SettingsNotebook (wx.Notebook) :
	"""
	GUI for listing all settings
	"""
	def __init__(self, parent):
		wx.Notebook.__init__(self, parent, DevSpectrometer, DevPulseShaper)
		
		self.GA = GATab (self) 
		self.AddPage (self.GA, "Genetic algorithm")
		
		self.Spectrometer = SpectrometerTab(self, DevSpectrometer)
		self.AddPage (self.Spectrometer, "OO Spectrometer")
		 
		self.PulseShaper = PulseShaperTab(self, DevPulseShaper)
		self.AddPage (self.PulseShaper, "Pulse shaper settings")
		
		# Dictionary to bind names to tabs for saving and loading settings
		self.settings_to_tabs = {"GA" : self.GA, "Spectrometer" : self.Spectrometer, "PulseShaper" : self.PulseShaper }
			
########################################################################

class GetTransformLimitedPhase (BasicWindow) :

	def __init__ (self, parent) :
		# Starting spectrometer
		self.Spectrometer = ManagerSpectrometer()
		self.SpectrometerProc = self.Spectrometer.start()
		
		# Starting pulse shaper
		self.PulseShaper = ManagerShaper()
		self.PulseShaperProc = self.PulseShaper.start()
		
		# Create GUI
		dw, dh = wx.DisplaySize()
		wx.Frame.__init__ (self, parent, title="Get transform limited phase by optimizing SHG",
								size=(0.9*dw, 0.88*dh) )
		
		self.ConstructGUI ()
		self.Center()
		self.Show ()
		wx.EVT_CLOSE (self, self.on_close)
		
	def __del__ (self) :	
		# Close spectrometer
		self.Spectrometer.exit(); self.SpectrometerProc.join() 
		
		# Close pulse shaper
		self.PulseShaper.exit(); self.PulseShaperProc.join()
		
	def ConstructGUI (self) :
		""" Build GUI """
		self.panel = wx.Panel(self)
		sizer = wx.GridBagSizer ()
		
		############################ Settings Notebook ############################
		self.SettingsNotebook = SettingsNotebook(self.panel, self.Spectrometer, self.PulseShaper)
		sizer.Add(self.SettingsNotebook, pos=(0, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT , border=10)

		############################ Command panel ############################
		boxsizer = wx.BoxSizer (wx.VERTICAL)
		
		# Interactively display spectrum
		boxsizer.Add (self.CreateShowSpectrumButton(), flag=wx.EXPAND, border=5)
		
		# Separator
		boxsizer.Add (wx.StaticText(self.panel), flag=wx.EXPAND, border=5)
		
		################## Calibrate button ##################
		self.find_tl_phase_button = wx.Button (self.panel)
		self.find_tl_phase_button.Bind (wx.EVT_LEFT_DOWN, self.StartSearch_TransformLimitedPhase )
		self.find_tl_phase_button.Bind (wx.EVT_LEFT_DCLICK, self.StartSearch_TransformLimitedPhase)
		boxsizer.Add(self.find_tl_phase_button, flag=wx.EXPAND, border=5)
		# Define labels
		self.find_tl_phase_button.__start_label__ 	= "Find TL phase"
		self.find_tl_phase_button.__pause_label__ 	= "PAUSE optimization"
		self.find_tl_phase_button.__resume_label__	= "RESUME optimization"
		self.find_tl_phase_button.__stop_label__ 	= "STOP optimization"
		self.find_tl_phase_button.SetLabel (self.find_tl_phase_button.__start_label__)

		###############	 Save transform limited phase ####################
		self.save_tl_phase_button = wx.Button(self.panel, label="Save TL phase")
		self.save_tl_phase_button.Bind (wx.EVT_BUTTON, self.SaveTransformLimitedPhase)
		boxsizer.Add(self.save_tl_phase_button, flag=wx.EXPAND, border=5)
		
		# Separator
		boxsizer.Add (wx.StaticText(self.panel), flag=wx.EXPAND, border=5)
		
		# Send random phase to the pulse shaper
		boxsizer.Add (self.CreateRandomPhaseButton(), flag=wx.EXPAND, border=5)
		# Send random amplitude to the pulse shaper
		boxsizer.Add (self.CreateRandomAmplitudeButton(), flag=wx.EXPAND, border=5)
		# Send zero amplitude and zero phase to the pulse shaper
		boxsizer.Add (self.CreateZeroAmplitudeButton(), flag=wx.EXPAND, border=5)
		# Open pulse shaper equalizer
		boxsizer.Add (self.CreatePulseShaperEqualizerButton(), flag=wx.EXPAND, border=5)
		
		# Separator
		boxsizer.Add (wx.StaticText(self.panel), flag=wx.EXPAND, border=5)
		
		# Save settings
		boxsizer.Add( self.CreateSaveSettingsButton(), flag=wx.EXPAND, border=5)
		# Load settings
		boxsizer.Add( self.CreateLoadSettingsButton(), flag=wx.EXPAND, border=5)
		
		sizer.Add(boxsizer, pos=(1, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.GROW, border=10)
		########################### End of constructing panel ######################################
		self.panel.SetSizer (sizer)
		
		############################# Setting visvis #######################################
		Figure = app.GetFigureClass()
		self.fig = Figure(self)
		
		boxsizer = wx.BoxSizer (wx.HORIZONTAL)
		boxsizer.Add(self.panel, 0.5, wx.EXPAND)
		boxsizer.Add(self.fig._widget, 2, wx.EXPAND)
		
		#########################################################################################			
		self.SetSizer (boxsizer)
		self.SetAutoLayout(True)
		self.Layout() 
	
	def SaveTransformLimitedPhase (self, event) :
		"""
		self.save_tl_phase_button.Bind (wx.EVT_BUTTON, self.)
		"""
		pass
		
	def StartSearch_TransformLimitedPhase (self, event) :
		"""
		Button `self.find_tl_phase_button` was cliked. Use GA to find the transfrom limited phase
		"""
		button = self.find_tl_phase_button
		try :
			# Mouse double clicking stops scanning
			if event.GetEventType() == wx.wxEVT_LEFT_DCLICK  : button.SetLabel (button.__stop_label__)
		except AttributeError : pass
			
		if button.GetLabel() == button.__start_label__ :
			self.StopAllJobs ()
			
			# Ask user to select the file
			openFileDialog = wx.FileDialog(self, button.__start_label__, "", 
							"transfrom_limited_phase.hdf5", "HDF5 files (*.hdf5)|*.hdf5", 
							wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT | wx.FD_CHANGE_DIR)
			
			# Check whether user cancelled
			if openFileDialog.ShowModal() == wx.ID_CANCEL: return None	 
			self.optimization_log_filename = openFileDialog.GetPath()
			
			# Save the global settings
			self.SaveSettings(filename=self.optimization_log_filename, OpenDialog=False)
			
			# get spectrometer's settings
			settings = self.SettingsNotebook.Spectrometer.GetSettings()
			# Initiate spectrometer
			if self.Spectrometer.SetSettings(settings) == RETURN_FAIL : return
			# Initiate pulse shaper
			settings = self.SettingsNotebook.PulseShaper.GetSettings()
			self.PulseShaper.StopDevice()
			if self.PulseShaper.Initialize(settings) == RETURN_FAIL : return
			
			# Open the file for optimization log
			with  h5py.File (self.optimization_log_filename, 'a') as optimization_log_file :
			
				# The HDF5 group where each iteration is picked
				try : 
					optimization_iterations = optimization_log_file["optimization_iterations"]
					# loading last iteration, i.e., checkpoint
					self.current_iteration_number = max( int(key) for key in optimization_iterations.keys() )
					checkpoint = optimization_iterations[ str(self.current_iteration_number) ]
					self.current_iteration_number += 1
					# ask user whether to resume optimization
					if wx.MessageDialog(self, "Should optimization be resumed?", caption="resume optimization", 
						style=wx.YES_NO | wx.YES_DEFAULT ).ShowModal() == wx.ID_NO : raise ValueError
				except (KeyError, ValueError) :
					# Start fresh optimization, since there no iteration is pickled
					checkpoint = None
					self.current_iteration_number = 0
					# delete if it exists
					try : del optimization_log_file["optimization_iterations"]
					except KeyError : pass
					# create empty group
					optimization_log_file.create_group ("optimization_iterations")
				
				
				# Get number of optimization variables 
				self.num_pixels = self.PulseShaper.GetParamNumber()
				if self.num_pixels == RETURN_FAIL :
					raise RuntimeError ("Optimization cannot be started since calibration file was not loaded")
				
				# Setting up GA
				self.GA = PixelWiseGA(self.num_pixels, self.SettingsNotebook.GA.GetSettings(),
						self.FitnessFunction, checkpoint)
						 
				# The HDF5 group where optimization parameters are stored
				try : del optimization_log_file["optimization_settings"]
				except KeyError : pass
				optimization_settings_group = optimization_log_file.create_group ("optimization_settings")
				
				# Get wavelengths
				try :
					self.wavelengths = self.Spectrometer.GetWavelengths()
					optimization_settings_group["wavelengths"] = self.wavelengths
				except AttributeError : self.wavelengths = None
				
			# Changing the button's label 
			button.SetLabel (button.__pause_label__)
			
			# Cleaning visvis objects
			try : del self.fitness_plot
			except AttributeError : pass
			
			# Start optimization
			self.pause_optimization = False
			wx.CallAfter(self.NextIteration)
	
		elif button.GetLabel() == button.__pause_label__ :
			self.pause_optimization = True; button.SetLabel (button.__resume_label__)
		
		elif button.GetLabel() == button.__resume_label__ :
			self.pause_optimization = False
			wx.CallAfter(self.NextIteration)
			button.SetLabel (button.__pause_label__)

		elif button.GetLabel() == button.__stop_label__ :
			del self.pause_optimization
			button.SetLabel (button.__start_label__)
			
		else : raise ValueError ("Unrecognised button's label")
				
	
	def NextIteration (self) :
		"""
		Perform one optimization iteration
		"""
		# pause optimization, if requested
		try :
			if self.pause_optimization : return
		except AttributeError : 
			# Halting optimization
			return
		
		# Perform the optimization iteration
		self.GA.NextIteration()
		
		self.SaveOptimization()
		self.DisplayOptimization()
		
		# Doing the next iteration
		wx.CallAfter(self.NextIteration)
		
	def SaveOptimization (self) : 
		"""
		Save current iteration into the log file
		"""
		wx.Yield()

		# Open the file for optimization log
		with  h5py.File (self.optimization_log_filename, 'a') as optimization_log_file :
			checkpoint = optimization_log_file["optimization_iterations"]
			checkpoint = checkpoint.create_group( str(self.current_iteration_number) )
	
			# Pickle optimization 
			for key, value in self.GA.GetCheckPoint().items() :
				checkpoint[key] = value
			
			# Saving population in HDF5 format in the appropriate HDF5 group 
			individuals = checkpoint.create_group("individuals")
			for num, ind in enumerate(self.GA.optimization_pop) :
				ind_grp = individuals.create_group( str(num) )
				ind_grp["pulse_shape"] 	= ind
				ind_grp["fitness"]		= ind.fitness.values
				
			
		# Going to next iteration
		self.current_iteration_number += 1
		
	def DisplayOptimization (self) :
		"""
		Display the progress of optimization
		"""
		wx.Yield()
		
		visvis.cla(); visvis.clf()
		
		optimization_log = self.GA.GetOptimizationLog()
		
		for values, color in izip_longest( optimization_log.itervalues(), ['r', 'g', 'b', 'k'],  fillvalue='y' ) :
			try : visvis.plot ( values, lc=color ) 
			except Exception : pass
		
		visvis.xlabel ('iteration')
		visvis.ylabel ('objective function')
		visvis.legend( optimization_log.keys() ) 
		
		wx.Yield()
	
	def FitnessFunction (self, individual) :
		"""
		This is the objective function to be maximized.
		`individual` is supplied by the optimization algorithm
		"""
		phase = np.array(individual)
		self.PulseShaper.SetAmplPhase( np.ones(len(phase)), phase  )
	
		wx.Yield()
		
		# Get spectrum
		spectrum = self.Spectrometer.AcquiredData() 
		
		spectrum = gaussian_filter(spectrum, sigma=1)
		
		return spectrum.sum() - spectrum.min()*spectrum.size,  
		
	
		
#########################################################################
if __name__ == '__main__' :
	app = visvis.use('wx')
	app.Create()
	GetTransformLimitedPhase (None)
	app.Run()