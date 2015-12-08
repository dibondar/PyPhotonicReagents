"""
Minimalistic GUI
"""
# Add main directory to enable imports  
if __name__ == '__main__' :
	import os
	os.sys.path.append(os.path.abspath('..'))

import wx, h5py, time
import numpy as np
from scipy.interpolate import pchip_interpolate, PchipInterpolator
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import argrelmin

# Real time plotting
import visvis

# GUI components
from libs.gui.basic_window import BasicWindow, SaveSettings

# Hardware
from libs.dev.pulse_shaper import *

from libs.dev.camera_istar import ManagerIStarCamera as ManagerSpectrometer
from libs.dev.camera_istar import IStarCameraTab as SpectrometerTab

# Uncomment these lines if you want to use OceanOptics instead of IStar specified above
#from libs.dev.spectrometer_ocean_optics import ManagerOceanOpticsSpectrometer as ManagerSpectrometer
#from libs.dev.spectrometer_ocean_optics import OceanOpticsSpectrometerTab as SpectrometerTab

########################################################################

class SettingsNotebook (wx.Notebook) :
	"""
	GUI for listing all settings
	"""
	def __init__(self, parent, DevSpectrometer, DevPulseShaper):
		wx.Notebook.__init__(self, parent)

		self.Spectrometer = SpectrometerTab(self, DevSpectrometer)
		self.AddPage (self.Spectrometer, "Spectrometer settings")
		 
		self.PulseShaper = PulseShaperTab(self, DevPulseShaper)
		self.AddPage (self.PulseShaper, "Pulse shaper settings")

		# Dictionary to bind names to tabs for saving and loading settings
		self.settings_to_tabs = {"Spectrometer" : self.Spectrometer, 
			"PulseShaper" : self.PulseShaper}
			
########################################################################

class CalibrateShaper (BasicWindow) :

	def __init__ (self, parent) :
		# Starting spectrometer
		self.Spectrometer = ManagerSpectrometer()
		self.SpectrometerProc = self.Spectrometer.start()
		
		# Starting pulse shaper
		self.PulseShaper = ManagerShaper()
		self.PulseShaperProc = self.PulseShaper.start()
		
		# Create GUI
		dw, dh = wx.DisplaySize()
		wx.Frame.__init__ (self, parent, title="Pulse shaper calibration with Ocean Optics Spectrometer",
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

		self.ActionButton = wx.Button (self.panel)
		self.ActionButton.__start_label__ = "start action"
		self.ActionButton.__stop_label__ = "STOP action"
		self.ActionButton.SetLabel (self.ActionButton.__start_label__)
		self.Bind (wx.EVT_BUTTON, self.SartAction, self.ActionButton)
		boxsizer.Add (self.ActionButton, flag=wx.EXPAND, border=5)

		# Separator
		boxsizer.Add (wx.StaticText(self.panel), flag=wx.EXPAND, border=5)

		# Interactively display spectrum
		boxsizer.Add (self.CreateShowSpectrumButton(), flag=wx.EXPAND, border=5)

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

	def SartAction (self, event) :
		"""
		<self.ActionButton> was clicked
		"""
		button = self.ActionButton
		if button.GetLabel() == button.__start_label__ :
			self.StopAllJobs ()

			# get spectrometer's settings
			settings = self.SettingsNotebook.Spectrometer.GetSettings()
			# Initiate spectrometer
			if self.Spectrometer.SetSettings(settings) == RETURN_FAIL : return
			# Initiate pulse shaper
			settings = self.SettingsNotebook.PulseShaper.GetSettings()
			if self.PulseShaper.Initialize(settings) == RETURN_FAIL : return

			# Get wavelengths
			self.wavelengths = self.Spectrometer.GetWavelengths()

			# Start doing some actiom
			self._continue_action = True
			wx.CallAfter(self.DoAction)

			# Change button's label
			button.SetLabel (button.__stop_label__)

		elif button.GetLabel() == button.__stop_label__ :
			self._continue_action = False
			button.SetLabel (button.__start_label__)

		else : raise ValueError("Label is not recognized")


	def DoAction (self):
		"""
		This method is affiliated to the method <self.SartAction>
		"""
		####################################################
		#
		# Zak add you code here
		# e.g., in what follows we just acquire spectra
		#
		####################################################

		# If you want to apply some radom phase maske then uncomment the following:
		#amplitude_mask = self.PulseShaper.GetZeroMask() + 1
		#phase_mask = np.random.rand(max_amplitude_mask.size)
		#self.PulseShaper.SetMasks(amplitude_mask, phase_mask)

		# Getting spectrum
		spectrum = self.Spectrometer.AcquiredData()
		
		# Plot the spectra
		visvis.gca().Clear()
		visvis.plot (self.wavelengths, spectrum)
		visvis.xlabel("wavelength (nm)")
		visvis.ylabel("counts")
		visvis.title ("Spectrum from radom pulse shape")
		self.fig.DrawNow()
			
		# Going to the next iteration
		if self._continue_action:
			wx.CallAfter (self.DoAction)

		
#########################################################################
if __name__ == '__main__' :
	app = visvis.use('wx')
	app.Create()
	CalibrateShaper (None)
	app.Run()