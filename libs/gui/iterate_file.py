########################################################################
#
#	The abstract base class for iterating over a data set in the HDF5 file.
#	This class can be used, .e.g., to view scans
#
########################################################################

import wx
import h5py
import ntpath

import matplotlib
matplotlib.use('WXAgg')
from matplotlib.pyplot import cm
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

########################################################################

class IterateFile (wx.Frame) :
	"""
	The abstract base class for iterating over a data set in the HDF5 file.
	This class can be used, .e.g., to view scans.
	
	To use this abstract class the method <self.UpdateFrame(event=None)> must be defined.
	This method is called when the current frame is changed. 
	
	The optional method <self.IniLoad (hdf5_file)> will be called after the HDF5 file was 
	successfully opened to load desired data.
	
	The first step in <self.UpdateFrame> should be a call to <self.GetCurrentFrame>
	"""
	def __init__ (self, groupename, parent=None, filename = None, title=None) :
		"""
		<filename> the HDF5 file name which will be iterated over
		<groupename> the HDF5 group name where trainable data sets are stored.  
			The data set over which iterations is be performed must be labels by integers.
		"""
		
		# If file name is not specified, then ask user what file should be open 
		if filename == None :
			openFileDialog = wx.FileDialog (parent, "Chose HDF5 file for viewing", "", "", \
				"HDF5 files (*.hdf5)|*.hdf5", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR)

			# Check whether user cancelled
			if openFileDialog.ShowModal() == wx.ID_CANCEL : 
				raise ValueError ("User did not chose the file to view")
			else : filename = openFileDialog.GetPath()
		
		# Open the HDF5 file
		self.data_file = h5py.File (filename, 'r')
		
		# Loading other data, if desired
		try : self.IniLoad (self.data_file)
		except AttributeError : pass
		
		# Saving the group over which the iteration happens  
		self.iteration_group = self.data_file[groupename]
		
		# Converting the data set names in the itterable group to a list of int
		self.frames = sorted( self.iteration_group.keys(), key=int )
		
		# Construct GUI
		if title == None : 
			title = "Viewing file [%s]" % ntpath.basename(filename)

		wx.Frame.__init__ (self, parent, title=title)
		self.ConstructGUI ()
		self.Show ()
		self.Maximize (True)
		self.UpdateFrame ()
		
	def ConstructGUI (self) :
		"""
		Create GUI
		"""
		panel = wx.Panel(self)
		sizer = wx.BoxSizer (wx.VERTICAL)
		
		################################ Animation control button ###################################
		boxsizer = wx.BoxSizer (wx.HORIZONTAL)

			# Go to initial frame button
		self.initial_frame_button = wx.Button (panel, label="<<")
		def GoToInitialFrame (event) :
			self.current_frame.SetValue (0); self.UpdateFrame ()
		self.Bind (wx.EVT_BUTTON, GoToInitialFrame, self.initial_frame_button)
		boxsizer.Add(self.initial_frame_button, flag=wx.LEFT, border=5)

		# Go to the previous frame
		self.previous_frame_button = wx.Button (panel, label="<")
		def GoToPreviousFrame (event) :
			current_value = self.current_frame.GetValue()
			if current_value > 0 : 
				self.current_frame.SetValue(current_value-1); self.UpdateFrame ()
		self.Bind (wx.EVT_BUTTON, GoToPreviousFrame, self.previous_frame_button)
		boxsizer.Add(self.previous_frame_button, flag=wx.LEFT, border=5)
		
		# Variable storing current frame number
		self.current_frame = wx.SpinCtrl (panel, value="0", min=0, max=len(self.frames)-1)
		self.current_frame.Bind (wx.EVT_SPINCTRL, self.UpdateFrame)
		boxsizer.Add(self.current_frame, flag=wx.LEFT, border=5)

		# Animation button
		self.animation_button = wx.Button (panel)
		self.animation_button.__start_label = "Start animation"
		self.animation_button.__stop_label = "STOP animation"
		self.animation_button.SetLabel (self.animation_button.__start_label)
		self.Bind (wx.EVT_BUTTON, self.OnAnimation, self.animation_button)
		boxsizer.Add(self.animation_button, flag=wx.LEFT, border=5)

		# Go to the next frame button
		self.next_frame_button = wx.Button (panel, label=">")
		def GoToNextFrame (event) :
			current_value = self.current_frame.GetValue()
			if current_value < len(self.frames)-1 : 
				self.current_frame.SetValue(current_value+1); self.UpdateFrame ()
		self.Bind (wx.EVT_BUTTON, GoToNextFrame, self.next_frame_button)
		boxsizer.Add(self.next_frame_button, flag=wx.LEFT, border=5)

		# Go to the last frame button
		self.final_frame_button = wx.Button (panel, label=">>")
		def GoToLastFrame (event) :
			self.current_frame.SetValue (len(self.frames)-1); self.UpdateFrame ()
		self.Bind (wx.EVT_BUTTON, GoToLastFrame, self.final_frame_button)
		boxsizer.Add(self.final_frame_button, flag=wx.LEFT, border=5)

		sizer.Add(boxsizer, flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.GROW, border=5)
		
		################################ Canvas for plotting ##############################
		# Matplotlib canvas
		boxsizer = wx.BoxSizer (wx.VERTICAL)
		
		self.dpi = 80
		display_width, display_hight = wx.DisplaySize()
		
		self.fig = Figure((0.49*display_width/self.dpi, 0.8*display_hight/self.dpi), dpi=self.dpi)
		self.canvas = FigCanvas (panel, -1, self.fig)
		
		boxsizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
		boxsizer.Add(NavigationToolbar(self.canvas), 0, wx.EXPAND)
		
		sizer.Add(boxsizer, flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.GROW, border=5)	
		#########################################################################################			
		panel.SetSizer (sizer)
		
		
	def GetCurrentFrame (self) :
		"""
		Load the data set
		"""
		return self.iteration_group[ self.frames[self.current_frame.GetValue()] ]
	
	def OnAnimation (self, event=None) :
		"""
		<self.animation_button> was clicked
		"""
		if self.animation_button.GetLabel() == self.animation_button.__start_label :
			
			def DoAnimation () :
				wx.Yield()
				current_value = self.current_frame.GetValue()
				if current_value < len(self.frames)-1 : 
					# Continue animation
					self.current_frame.SetValue(current_value+1)
					self.UpdateFrame()
					
					# Decide whether to continue animation
					if self.animation_button.GetLabel() == self.animation_button.__stop_label :
						wx.CallAfter(DoAnimation)
				else :
					# Stop animation
					self.OnAnimation()
			
			# Initiate animation
			wx.CallAfter(DoAnimation)
			self.animation_button.SetLabel( self.animation_button.__stop_label )
		else :
			# Stop animation
			self.animation_button.SetLabel( self.animation_button.__start_label ) 
		