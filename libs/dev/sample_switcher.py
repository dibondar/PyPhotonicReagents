########################################################################
#
#	This file contains classes for controlling sample Switcher
#
########################################################################

from libs.gui.hardware_control import HardwareGUIControl
from libs.dev.basic_device import BasicDevice
from libs.dev.consts import * 

import wx
from wx.lib.agw.floatspin import FloatSpin as wxFloatSpin
import multiprocessing
import serial
import visvis

import numpy as np
from scipy.ndimage.filters import gaussian_filter

########################################################################
#
#	Manager of the sample Switcher
#
########################################################################

class ManagerSampleSwitcher :
	"""
	Class that manges the sample Switcher
	"""
	def __init__ (self) :
		# Create the lock for device 
		self.lock = multiprocessing.Lock()
		# Create a pipe for communication
		self.parent_connection, self.child_connection = multiprocessing.Pipe()
	
	def __del__ (self) :
		self.parent_connection.close()
		self.child_connection.close()	
	
	def start(self) :
		"""
		Start the process controlling the shaper
		"""
		p = SampleSwitcher(self.child_connection)
		p.start()
		return p
	
	def run(self, command, arguments=None) :
		"""
		Send the command to the shaper through the pipe
		"""
		self.lock.acquire()
		self.parent_connection.send( (command, arguments) )
		result = self.parent_connection.recv()
		self.lock.release()
		return result
		
	def exit(self) :
		"""
		Close the process
		"""
		self.StopDevice()
		return self.run("Exit")
	
	def Initialize(self, settings) :
		"""
		Initialize shaper 
		"""
		return self.run("Initialize", settings)
	
	def StopDevice(self) :
		return self.run("StopDevice")
		
	def MoveTo(self, position) :
		return self.run("MoveTo", position)
	
	def MoveToChannel(self, channel_num) :
		return self.run("MoveToChannel", channel_num)
	
	def GetCurrentPosition(self) :
		return self.run("GetCurrentPosition")
		
	def GetChannelNum(self) :
		return self.run("GetChannelNum")
		
	
########################################################################
#
#	Process where the sample switcher resides
#
########################################################################

class SampleSwitcher (BasicDevice) : 
	"""
	Control sample Switcher by moving a translation stage
	"""
	def Initialize (self, settings) :
		# Close the port if it is already used
		try : del self.serial_port
		except AttributeError : pass
	
		# Start the communication port
		self.serial_port = serial.Serial (port=settings["port"], 
			baudrate=19200,  bytesize=8, parity=serial.PARITY_NONE, stopbits=1)
			
		# Wrapper to define '\r' as the end of line
		#self.serial_port = io.BufferedRWPair(self.serial_port, self.serial_port)
		#self.serial_port = io.TextIOWrapper(self.serial_port, newline='\r', line_buffering=True)
		
		# Save the location of channels 
		try :
			self.chanel_positions = eval( "( %s, )" % settings["chanel_positions"] )
		except (NameError, SyntaxError, KeyError), msg :
			print "SampleSwitcher Error: Positions of channels could not be localized (%s)" % msg
			self.chanel_positions = ()
		
		return RETURN_SUCCESS
		
	def MoveToChannel (self, channel_num) :
		"""
		Move to the channel specified by the number `channel_num`
		"""
		return self.MoveTo( self.chanel_positions[channel_num] )
	
	def GetChannelNum (self, settings=None) :
		"""
		Get number of calibrated channels 
		"""
		return len(self.chanel_positions)
	
	def MoveTo (self, position) :
		"""
		Move a linear translation stage to a specified absolute position in mm
		"""
		self.serial_port.write ("1PA%.6e;1WS0;1MD?\r" % position)
		if not int( self.serial_port.readline() ) : 
			print "Error: Moving stage is still in motion!"
			
		return RETURN_SUCCESS
	
	def GetCurrentPosition (self, settings=None) :
		"""
		Get current location of the moving state
		"""
		self.serial_port.write  ("1WS0;1TP\r")
		return float( self.serial_port.readline() )
	
	def StopDevice (self, arguments=None) :
		try : del self.serial_port
		except AttributeError : pass
		return RETURN_SUCCESS

########################################################################

class SampleSwitcherTab (HardwareGUIControl) : 
	"""
	This class represents a GUI controlling properties of sample Switcher.
	"""
	def __init__(self, parent, dev) :
		"""
		`dev` is a sample switcher manager
		"""
		HardwareGUIControl.__init__(self, parent, dev)
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		# Specify the communication port name
		sizer.Add (wx.StaticText(self, label="Communication port"), flag=wx.LEFT, border=5)
		port_name = wx.TextCtrl (self, value="COM1")
		port_name.__label__ = "port"
		sizer.Add (port_name, flag=wx.EXPAND, border=5)
		
		# List of positions of channels
		sizer.Add (wx.StaticText(self, label="Position of channels"), flag=wx.LEFT, border=5)
		self.chanel_positions_ctrl = wx.TextCtrl (self, value="", style=wx.TE_MULTILINE|wx.EXPAND)
		self.chanel_positions_ctrl.__label__ = "chanel_positions"
		sizer.Add (self.chanel_positions_ctrl, flag=wx.EXPAND, border=5)
		
		# Spacer
		sizer.Add (wx.StaticText(self), flag=wx.LEFT, border=5)
		
		##################### Parameters for automatic calibration #######################
		
		sb_sizer = wx.StaticBoxSizer( wx.StaticBox(self, label="Automatic calibration"),  wx.VERTICAL )
		
		# Starting position for scanning
		sb_sizer.Add ( wx.StaticText(self, label="Beginning of scanning range (mm)") )
		init_scan_pos_ctrl = wxFloatSpin (self, increment=0.01, value=16, digits=3)
		init_scan_pos_ctrl.__label__ = "initial_scan_position"
		sb_sizer.Add (init_scan_pos_ctrl,  flag=wx.EXPAND, border=5)
		
		# Final position for scanning
		sb_sizer.Add ( wx.StaticText(self, label="End of scanning range (mm)") )
		fin_scan_pos_ctrl = wxFloatSpin (self, increment=0.01, value=80., digits=3)
		fin_scan_pos_ctrl.__label__ = "final_scan_position"
		sb_sizer.Add (fin_scan_pos_ctrl,  flag=wx.EXPAND, border=5)
		
		# Scanning step size
		sb_sizer.Add ( wx.StaticText(self, label="Scanning step size (mm)") )
		scan_step_size_ctrl = wxFloatSpin (self, increment=0.001, value=0.08,  digits=3, min_val=0.001, max_val=1)
		scan_step_size_ctrl.__label__ = "scan_step"
		sb_sizer.Add (scan_step_size_ctrl,  flag=wx.EXPAND, border=5)
		
		# Threshold to recognize peaks 
		sb_sizer.Add ( wx.StaticText(self, label="Background signal cut-off") )
		background_cutoff_ctrl = wxFloatSpin (self, increment=0.01, value=0.9, digits=3, min_val=0, max_val=1)
		background_cutoff_ctrl.__label__ = "background_cutoff"
		sb_sizer.Add (background_cutoff_ctrl,  flag=wx.EXPAND, border=5)
		
		# Methods to find peaks
		sb_sizer.Add ( wx.StaticText(self, label="Peak finding method") )
		self.peak_finders = { 
			"total fluoresce" 				: 
				lambda x : np.copy(x),
			"log(total fluoresce)"			: 
				lambda x : np.log10(x),
			"diff(total fluoresce)"			: 
				lambda x : gaussian_filter( np.abs(gaussian_filter(x, sigma=1, order=1)), sigma=1 ),
			"log(diff(total fluoresce))"	: 
				lambda x : gaussian_filter( np.log10(np.abs(gaussian_filter(x, sigma=1, order=1))), sigma=1 ),
			"diff(log(total fluoresce))"	: 
				lambda x : gaussian_filter( np.abs(gaussian_filter(np.log10(x), sigma=1, order=1)), sigma=1 )
		}
		choices = list(self.peak_finders.keys())
		peak_finder_ctrl = wx.ComboBox (self, choices=choices, value=choices[0], style=wx.CB_READONLY )
		peak_finder_ctrl.__label__ = "peak_finder"
		sb_sizer.Add (peak_finder_ctrl,  flag=wx.EXPAND, border=5)	
		
		# Spacer
		sb_sizer.Add (wx.StaticText(self), flag=wx.LEFT, border=5)
		
		# Calibrate button
		calibrate_button = wx.Button (self)
		calibrate_button._start_label ="Scan and calibrate"
		calibrate_button._stop_label ="STOP calibration"
		calibrate_button.SetLabel (calibrate_button._start_label)
		calibrate_button._start_method = self.CalibrateSampleSwitcher
		calibrate_button._stop_method = self.StopCalibration
		calibrate_button.Bind (wx.EVT_BUTTON, calibrate_button._start_method)
		sb_sizer.Add (calibrate_button,  flag=wx.EXPAND, border=5)	
		
		# Re-analysed calibration data
		analyse_button = wx.Button (self, label="Re-calibrate without scanning")
		analyse_button.Bind (wx.EVT_BUTTON, self.AnalyzeTotalFluorescence)
		sb_sizer.Add (analyse_button,  flag=wx.EXPAND, border=5)
		
		sizer.Add(sb_sizer, flag=wx.EXPAND, border=5)
		# Spacer
		sizer.Add (wx.StaticText(self), flag=wx.LEFT, border=5)
		##################### Parameters for automatic calibration #######################
		
		sb_sizer = wx.StaticBoxSizer( wx.StaticBox(self, label="Move to"),  wx.VERTICAL )
		
		sb_sizer.Add ( wx.StaticText(self, label="Position") )
		self.moving_stage_position_ctrl = wxFloatSpin (self, increment=0.01, value=80., digits=5)
		self.moving_stage_position_ctrl.__label__ = "move_to_position"
		sb_sizer.Add (self.moving_stage_position_ctrl,  flag=wx.EXPAND, border=5)
		
		# Update current position button
		def OnGetCurrentPosition (event) :
			if self.dev.Initialize( self.GetSettings() ) == RETURN_FAIL : return
			self.moving_stage_position_ctrl.SetValue( self.dev.GetCurrentPosition() )
		
		get_current_possition_button = wx.Button (self, label="Get current position")
		get_current_possition_button.Bind (wx.EVT_BUTTON, OnGetCurrentPosition)
		sb_sizer.Add (get_current_possition_button,  flag=wx.EXPAND, border=5)
		
		# Move-to position button
		def OnMoveTo (event) :
			if self.dev.Initialize( self.GetSettings() ) == RETURN_FAIL : return
			self.dev.MoveTo( self.moving_stage_position_ctrl.GetValue() )
			
		move_to_button = wx.Button (self, label="Go to")
		move_to_button.Bind (wx.EVT_BUTTON, OnMoveTo)
		sb_sizer.Add (move_to_button,  flag=wx.EXPAND, border=5)
		
		# Spacer
		sb_sizer.Add (wx.StaticText(self), flag=wx.LEFT, border=5)
		
		# Go to sample button
		def OnMoveToSample (event) :
			if self.dev.Initialize( self.GetSettings() ) == RETURN_FAIL : return
			button = event.GetEventObject()
			# Get current sample number
			label_split = button.GetLabel().split()
			sample_num = int(label_split[-1]) 
			# Get next sample number
			sample_num = (sample_num + 1) % self.dev.GetChannelNum()
			# Go to the sample
			self.dev.MoveToChannel( sample_num )
			# Update the label
			label_split[-1] = str(sample_num)
			button.SetLabel( " ".join(label_split) )
			
		move_to_sample = wx.Button (self, label="Go to sample 0")
		move_to_sample.Bind (wx.EVT_BUTTON, OnMoveToSample)
		sb_sizer.Add (move_to_sample,  flag=wx.EXPAND, border=5)
		
		sizer.Add(sb_sizer, flag=wx.EXPAND, border=5)
		###################################################################################
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		self.CreateSettingsDict()
	
	def StopCalibration (self, event) :
		"""
		Abort calibration
		"""
		self.need_abort = True
		# Adjusting button's settings
		button = event.GetEventObject()
		button.SetLabel (button._start_label)
		button.SetBackgroundColour('')
		button.Bind( wx.EVT_BUTTON, button._start_method)
	
	def CalibrateSampleSwitcher (self, event) :
		"""
		`calibrate_button` was clicked.
		Perform the automated calibration.
		"""
		# Initiate spectrometer
		settings = self.parent.Spectrometer.GetSettings()
		if self.parent.Spectrometer.dev.SetSettings( settings ) == RETURN_FAIL : return
		
		# Initialize sample switching device
		settings = self.GetSettings()
		if self.dev.Initialize( settings ) == RETURN_FAIL : return
		
		initial_position = min(settings["initial_scan_position"], settings["final_scan_position"])
		final_position = max(settings["initial_scan_position"], settings["final_scan_position"])
		
		# Job started: Adjusting button's settings
		button = event.GetEventObject()
		button.SetLabel (button._stop_label)
		button.SetBackgroundColour('red')
		button.Bind( wx.EVT_BUTTON, button._stop_method)
		
		# Set's get started
		self.need_abort = False
		positions = []
		total_fluorescence = []
		
		for position in np.arange( initial_position, final_position, abs(settings["scan_step"]) ) :
			# Move to
			self.dev.MoveTo (position)
			# Save total intensity
			total_fluorescence.append( self.parent.Spectrometer.dev.AcquiredData().sum() )
			positions.append(position) 
			# Perform check every 10 steps
			if len(positions) % 10 == 0 : 
				wx.Yield()
				# abort, if requested 
				if self.need_abort : return
		
		# Saving measurements
		self.total_fluorescence = np.array(total_fluorescence)
		self.positions = np.array(positions)
		
		self.AnalyzeTotalFluorescence ()
		
		# Job finished: Adjusting button's settings
		button = event.GetEventObject()
		button.SetLabel (button._start_label)
		button.SetBackgroundColour('')
		button.Bind(wx.EVT_BUTTON, button._start_method)
		
	def AnalyzeTotalFluorescence (self, event=None) :
		"""
		`analyse_button` was clicked
		"""
		# Get current settings
		settings = self.GetSettings()
	
		# Apply peak finding filter
		signal = self.peak_finders[ settings["peak_finder"] ](self.total_fluorescence)
		
		# Scale to (0,1)
		signal -= signal.min()
		signal = signal / signal.max()
		
		##########################################################################
		
		# Partition signal into segments that are above the background noise  
		#signal = gaussian_filter(total_fluorescence, sigma=0.5)
		background_cutoff = settings["background_cutoff"]
		segments = [ [] ]
		for num, is_segment in enumerate( signal > background_cutoff ) :
			if is_segment : 
				# this index is in the segment
				segments[-1].append( num )
			elif len(segments[-1]) : # this condition is not to add empty segments
				# Start new segments
				segments.append( [] )
		
		# Find peaks as weighted average of the segment
		peaks = [ np.average(self.positions[S], weights=self.total_fluorescence[S]) for S in segments if len(S) ]
		
		##########################################################################
		
		# Saving the positions
		self.chanel_positions_ctrl.SetValue( ", ".join( "%2.4f" % p for p in peaks ) )
		
		##########################################################################
		
		# Plot acquired data
		visvis.cla(); visvis.clf()
		
		visvis.plot( self.positions, signal )
		visvis.plot( peaks, background_cutoff*np.ones(len(peaks)), ls=None, ms='+', mc='r', mw=20)
		visvis.plot( [self.positions.min(), self.positions.max()], [background_cutoff, background_cutoff], lc='r', ls='--', ms=None)
		
		visvis.legend( ["measured signal", "peaks found", "background cut-off"] )
		
		visvis.ylabel( "total fluorescence") 
		visvis.xlabel( 'position (mm)')
		
		