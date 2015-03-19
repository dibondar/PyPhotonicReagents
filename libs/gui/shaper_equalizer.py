from libs.dev.consts import *

import wx, multiprocessing
import numpy as np
import functools 
import itertools

########################################################################################################

class PulseShaperEqualizerWindow (wx.Frame) :
	"""
	Control of pulse shaper using scroll bars to set an arbitrary phase and amplitude.
	The pulse shaper is assumed to be initialized.
	"""
	def __init__ (self, parent, PulseShaper) :
		"""
		PulseShaper -- the manager of initialized pulse shaper
		"""
		# Saving pulse shaper manager
		self.PulseShaper = PulseShaper
		
		# Saving number of pixels
		self.NumPixels = self.PulseShaper.GetParamNumber()
		if self.NumPixels == RETURN_FAIL :
			raise RuntimeError ("Pulse shaper equalizer cannot be started since calibration file was not loaded")
		
		# Limit number of pixels `pixel_number_bound`
		pixel_number_bound = 75
		if self.NumPixels > pixel_number_bound :
			print "Pulse Shaper Equalizer: Original number of pixels %d was restricted to %d" % (self.NumPixels, pixel_number_bound)
			self.NumPixels = pixel_number_bound
			
		dw, dh = wx.DisplaySize()
		wx.Frame.__init__ (self, parent, title="Pulse shaper Equalizer", size=(0.8*dw, 0.8*dh) )
		
		self.ConstructGUI ()
		self.Center()
		self.Show ()
		self.PulseShapeChanged(None)

	def ConstructGUI(self) :
		self.panel = wx.Panel(self)
		
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		# Initialize controls
		self.amplitude_controls = [ wx.ScrollBar (self.panel, style=wx.SB_VERTICAL) for _ in range(self.NumPixels) ]
		self.phase_controls = [ wx.ScrollBar (self.panel, style=wx.SB_VERTICAL) for _ in range(self.NumPixels) ]
		
		FormatLabel = lambda text : "".join( ( "\n%s" % character for  character in text ) )
		
		#################################### Functional buttons #######################################

		# Amplitude controls to clipboard buttons
		amplitude_shape_clipboard = wx.Button(self.panel, label="Copy amplitude mask to clipboard")
		amplitude_shape_clipboard.Bind (wx.EVT_BUTTON, 
			functools.partial(self.Controls2Clipboard, controls=self.amplitude_controls) 
		)
		
		# Phase controls to clipboard buttons
		phase_shape_clipboard = wx.Button(self.panel, label="Copy phase mask to clipboard")
		phase_shape_clipboard.Bind (wx.EVT_BUTTON, 
			functools.partial(self.Controls2Clipboard, controls=self.phase_controls) 
		)
		
		# Clipboard to amplitude controls buttons
		clipboard_amplitude_shape = wx.Button(self.panel, label="Import from clipboard to amplitude mask")
		clipboard_amplitude_shape.Bind (wx.EVT_BUTTON, 
			functools.partial(self.Clipboard2Control, controls=self.amplitude_controls) 
		)
		
		# Clipboard to phase controls buttons
		clipboard_phase_shape = wx.Button(self.panel, label="Import from clipboard to phase mask")
		clipboard_phase_shape.Bind (wx.EVT_BUTTON, 
			functools.partial(self.Clipboard2Control, controls=self.phase_controls) 
		)
		
		# Random amplitude buttons
		random_ampl = wx.Button(self.panel, label="Random amplitude")
		random_ampl.Bind (wx.EVT_BUTTON, 
			lambda _ : self.Array2Controls(np.random.rand(self.NumPixels), self.amplitude_controls)
		)
		
		# Random phase buttons
		random_phase = wx.Button(self.panel, label="Random phase")
		random_phase.Bind (wx.EVT_BUTTON, 
			lambda _ : self.Array2Controls(np.random.rand(self.NumPixels), self.phase_controls)
		)
		
		# Set zero amplitude
		zero_ampl = wx.Button(self.panel, label="Zero amplitude")
		zero_ampl.Bind (wx.EVT_BUTTON, 
			lambda _ : self.Array2Controls(np.zeros(self.NumPixels), self.amplitude_controls)
		)
		
		# Set zero phase
		zero_phase = wx.Button(self.panel, label="Zero phase")
		zero_phase.Bind (wx.EVT_BUTTON, 
			lambda _ : self.Array2Controls(np.zeros(self.NumPixels), self.phase_controls)
		) 
		
		# Set max amplitude
		max_ampl = wx.Button(self.panel, label="Max amplitude")
		max_ampl.Bind (wx.EVT_BUTTON, 
			lambda _ : self.Array2Controls(np.ones(self.NumPixels), self.amplitude_controls)
		)
		
		# Set max phase
		max_phase = wx.Button(self.panel, label="Max phase")
		max_phase.Bind (wx.EVT_BUTTON, 
			lambda _ : self.Array2Controls(np.ones(self.NumPixels), self.phase_controls)
		)
		
		#########################################################################################
		grd_sizer = wx.GridSizer(rows=3, cols=5)
		grd_sizer.AddMany( [ 
			(amplitude_shape_clipboard, 0, wx.EXPAND ), (clipboard_amplitude_shape, 0, wx.EXPAND ), 
			(random_ampl, 0, wx.EXPAND ), (zero_ampl, 0, wx.EXPAND ), (max_ampl, 0, wx.EXPAND ),
			# New row
			(phase_shape_clipboard, 0, wx.EXPAND ), (clipboard_phase_shape, 0, wx.EXPAND ), 
			(random_phase, 0, wx.EXPAND ), (zero_phase, 0, wx.EXPAND ), (max_phase, 0, wx.EXPAND )		
		] )
		
		sizer.Add (grd_sizer, 0, wx.EXPAND, border=5)
		
		#################################### Scroll controls  #######################################
		grd_sizer = wx.GridSizer(rows=2, cols=self.NumPixels+1)
		
		# Amplitude controls panel 
		grd_sizer.Add( wx.StaticText(self.panel, label=FormatLabel("AMPLITUDE"), style=wx.ALIGN_CENTRE_HORIZONTAL), 0, wx.ALIGN_CENTRE_VERTICAL|wx.ALIGN_CENTRE_HORIZONTAL )
		for control in self.amplitude_controls :
			control.SetScrollbar(0, 2, 100, 2)
			control.Bind( wx.EVT_SCROLL_CHANGED, self.PulseShapeChanged)
			grd_sizer.Add(control, 0, wx.EXPAND )
		
		## Pixel number labels 
		#grd_sizer.Add( wx.StaticText(self.panel) )
		#for num in range(1,self.NumPixels+1) :
		#	grd_sizer.Add( wx.StaticText(self.panel, label=str(num),style=wx.ALIGN_CENTRE_HORIZONTAL), 0, wx.ALIGN_CENTRE_VERTICAL|wx.ALIGN_CENTRE_HORIZONTAL ) 
		
		# Phase controls panel 
		grd_sizer.Add( wx.StaticText(self.panel, label=FormatLabel("PHASE"), style=wx.ALIGN_CENTRE_HORIZONTAL), 0, wx.ALIGN_CENTRE_VERTICAL|wx.ALIGN_CENTRE_HORIZONTAL )
		for control in self.phase_controls :
			control.SetScrollbar(100, 2, 100, 2)
			control.Bind( wx.EVT_SCROLL_CHANGED, self.PulseShapeChanged)
			grd_sizer.Add(control, 0, wx.EXPAND )
		#########################################################################################
		sizer.Add (grd_sizer, wx.EXPAND)
		
		self.panel.SetSizer(sizer)
		
	def Controls2Array (self, controls) :
		"""
		Convert potions of scrollbar controls into [0,1] array 
		"""
		return np.array( [ 1. - float(ctrl.GetThumbPosition())/(ctrl.GetRange() - ctrl.GetThumbSize()) for ctrl in controls ] )
		
	def Array2Controls (self, array, controls) :
		"""
		Convert [0,1] array into  potions of scrollbar controls.
		This function is inverse to self.Controls2Array
		"""
		if len(controls) < len(array) :
			raise ValueError("PulseShaperEqualizer Error: Array size must not exceed number of controls")
		
		if len(controls) > len(array) :
			# Convert to numpy array
			array = np.array(array, copy=False) 
			# Linearly interpolate
			x = np.linspace(0., 1., len(controls))
			xp = np.linspace(0., 1., array.size)
			array 	= np.interp(x, xp, array)
		
		# Set values of the controls
		for ctrl, val in itertools.izip(controls, array) :
			ctrl.SetThumbPosition( (ctrl.GetRange() - ctrl.GetThumbSize())*(1 - val) )
		
		# Send data to pulse shaper
		self.PulseShapeChanged(None)
		
	def PulseShapeChanged (self, event) :
		"""
		This method is called if a scrollbar controlling the pulse shaper mask is moved.  
		"""
		self.PulseShaper.SetAmplPhase( self.Controls2Array(self.amplitude_controls), self.Controls2Array(self.phase_controls) )
		
	def  Controls2Clipboard (self, event, controls) :
		"""
		Copy the values of controls into the Clipboard as text 
		"""
		if wx.TheClipboard.Open() :
			txt_data = str( self.Controls2Array(controls) )[1:-1]
			wx.TheClipboard.SetData(wx.TextDataObject(txt_data)) 
			wx.TheClipboard.Close()
		else :
			wx.MessageBox("Unable to open the clipboard", "Error")

	def Clipboard2Control (self, event, controls) :
		"""
		Copy the values in clipboard to the representation of controls
		"""
		if wx.TheClipboard.Open() :
		
			# Copy string from Clipboard
			#if wx.TheClipboard.IsSupported(wx.DF_TEXT):
			data = wx.TextDataObject()
			wx.TheClipboard.GetData(data)		
			wx.TheClipboard.Close()
			
			# Convert string to array
			data = map(np.float, data.GetText().split())
			
			# Set the controls
			self.Array2Controls(data, controls) 
		else :
			wx.MessageBox("Unable to open the clipboard", "Error")

########################################################################################################

class PulseShaperEqualizer (multiprocessing.Process):
	"""
	Process which opens pulse shaper equalizer
	"""
	def __init__ (self, PulseShaper) :
		"""
		`PulseShaper` -- pulse shaper manager 
		"""
		multiprocessing.Process.__init__(self)
		self.PulseShaper = PulseShaper

	def run (self) :
		app = wx.PySimpleApp()
		PulseShaperEqualizerWindow(None, self.PulseShaper)
		app.MainLoop()