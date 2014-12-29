from libs.dev.consts import *

import wx, multiprocessing
import numpy as np

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
		
		dw, dh = wx.DisplaySize()
		wx.Frame.__init__ (self, parent, title="Pulse shaper Equalizer", size=(0.8*dw, 0.8*dh) )
		
		self.ConstructGUI ()
		self.Center()
		self.Show ()
		self.PulseShapeChanged(None)

	def ConstructGUI(self) :
		self.panel = wx.Panel(self)
		sizer = wx.GridSizer(rows=3, cols=self.NumPixels+1)
		
		# Initialize controls
		self.amplitude_controls = [ wx.ScrollBar (self.panel, style=wx.SB_VERTICAL) for _ in range(self.NumPixels) ]
		self.phase_controls = [ wx.ScrollBar (self.panel, style=wx.SB_VERTICAL) for _ in range(self.NumPixels) ]
		
		FormatLabel = lambda text : "".join( ( "\n%s" % character for  character in text ) )
		
		sizer.Add( wx.StaticText(self.panel, label=FormatLabel("AMPLITUDE"), style=wx.ALIGN_CENTRE_HORIZONTAL), 0, wx.ALIGN_CENTRE_VERTICAL|wx.ALIGN_CENTRE_HORIZONTAL )
		for control in self.amplitude_controls :
			control.SetScrollbar(0, 2, 100, 2)
			control.Bind( wx.EVT_SCROLL_CHANGED, self.PulseShapeChanged)
			sizer.Add(control, 0, wx.EXPAND )
		
		# Number the pixels
		sizer.Add( wx.StaticText(self.panel) )
		for num in range(1,self.NumPixels+1) :
			sizer.Add( wx.StaticText(self.panel, label=str(num),style=wx.ALIGN_CENTRE_HORIZONTAL), 0, wx.ALIGN_CENTRE_VERTICAL|wx.ALIGN_CENTRE_HORIZONTAL ) 
		
		sizer.Add( wx.StaticText(self.panel, label=FormatLabel("PHASE"), style=wx.ALIGN_CENTRE_HORIZONTAL), 0, wx.ALIGN_CENTRE_VERTICAL|wx.ALIGN_CENTRE_HORIZONTAL )
		for control in self.phase_controls :
			control.SetScrollbar(100, 2, 100, 2)
			control.Bind( wx.EVT_SCROLL_CHANGED, self.PulseShapeChanged)
			sizer.Add(control, 0, wx.EXPAND )
		
		#########################################################################################
		self.panel.SetSizer(sizer)
		
	def Controls2Array (self, controls) :
		"""
		Convert potions of scrollbar controls into [0,1] array 
		"""
		return np.array( [ 1. - float(ctrl.GetThumbPosition())/(ctrl.GetRange() - ctrl.GetThumbSize()) for ctrl in controls ] )
	
	def PulseShapeChanged (self, event) :
		"""
		This method is called if a scrollbar controlling the pulse shaper mask is moved.  
		"""
		self.PulseShaper.SetAmplPhase( self.Controls2Array(self.amplitude_controls), self.Controls2Array(self.phase_controls) )

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