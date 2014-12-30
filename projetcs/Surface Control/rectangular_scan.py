from hardware_control import HardwareGUIControl
import wx, multiprocessing, itertools, h5py, time
import numpy as np

########################################################################
class ManagerRectangularScan :
	"""
	Managing class for performing rectangular scan of a sample
	"""
	def __init__ (self, settings, Spectrometer, MovingStageX, MovingStageY) :
		"""
		Start rectangular scanning
		"""
		self.settings = settings
		
		# Creating the events
		self.ScanningEvent = multiprocessing.Event()
		self.ScanningEvent.set()
		
		self.PauseScanningEvent = multiprocessing.Event()
		self.PauseScanningEvent.clear()
		
		# Starting the scanning process
		self.ScanProc = ProcRectangularScan( self.ScanningEvent, self.PauseScanningEvent, 
			settings, Spectrometer, MovingStageX, MovingStageY)
		self.ScanProc.start()
	
	def is_running(self) :
		"""
		Check whether the scanning is running
		"""
		return self.ScanningEvent.is_set()
	
	def pause(self) :
		"""
		Pause scanning 
		"""
		self.PauseScanningEvent.set()
		
	def resume(self) :
		"""
		Resume scanning 
		"""
		self.PauseScanningEvent.clear()
	
	def stop(self) :
		"""
		Stop scanning
		"""
		self.ScanningEvent.clear()
		# Post process the scan
		multiprocessing.Process(target=PostProcessRectangularScan, args=(self.settings["filename"],)).start()
		
	def __del__ (self) :
		self.ScanProc.join()
		
########################################################################
#
#	Scanning GUI settings
#
########################################################################		

class RectangularScanTab (HardwareGUIControl) :
	"""
	Setting of rectangular scanning
	"""
	def __init__ (self, parent) :
		HardwareGUIControl.__init__(self, parent)
	
		sizer = wx.BoxSizer(wx.VERTICAL)
	
		################################################
		sizer.Add (wx.StaticText(self, label="File name"), flag=wx.LEFT, border=5)
		filename_ctr = wx.TextCtrl (self, value="rectangular_scan.hdf5")
		filename_ctr.__label__ = "filename"
		sizer.Add (filename_ctr, flag=wx.EXPAND, border=10)
		
		###################### Scanning geometry parameters ######################
		sizer.Add (wx.StaticText(self, label="\nx min (mm)"), flag=wx.LEFT, border=5)
		x_min_ctr = wx.SpinCtrl (self, value="4", min=-50, max=50)
		x_min_ctr.SetLabel ("x min")
		sizer.Add (x_min_ctr, flag=wx.EXPAND, border=10)
		
		sizer.Add (wx.StaticText(self, label="x max (mm)"), flag=wx.LEFT, border=5)
		x_max_ctr = wx.SpinCtrl (self, value="6", min=-50, max=50)
		x_max_ctr.SetLabel ("x max")
		sizer.Add (x_max_ctr, flag=wx.EXPAND, border=10)
		
		sizer.Add (wx.StaticText(self, label="y min (mm)"), flag=wx.LEFT, border=5)
		y_min_ctr = wx.SpinCtrl (self, value="4", min=-50, max=50)
		y_min_ctr.SetLabel ("y min")
		sizer.Add (y_min_ctr, flag=wx.EXPAND, border=10)
		
		sizer.Add (wx.StaticText(self, label="y max (mm)"), flag=wx.LEFT, border=5)
		y_max_ctr = wx.SpinCtrl (self, value="6", min=-50, max=50)
		y_max_ctr.SetLabel ("y max")
		sizer.Add (y_max_ctr, flag=wx.EXPAND, border=10)
		
		###################### Number of points ######################
		sizer.Add (wx.StaticText(self, label="\nx points"), flag=wx.LEFT, border=5)
		x_points_ctr = wx.SpinCtrl (self, value="10", min=1, max=1e3)
		x_points_ctr.SetLabel ("x points")
		sizer.Add (x_points_ctr, flag=wx.EXPAND, border=10)
		
		sizer.Add (wx.StaticText(self, label="y points"), flag=wx.LEFT, border=5)
		y_points_ctr = wx.SpinCtrl (self, value="10", min=1, max=1e3)
		y_points_ctr.SetLabel ("y points")
		sizer.Add (y_points_ctr, flag=wx.EXPAND, border=10)
		
		###################### Fast varying axes ######################
		sb = wx.StaticBox(self, label="\n\nSlow varying axis")
		boxsizer = wx.StaticBoxSizer(sb, wx.HORIZONTAL)
		
		x_slow = wx.RadioButton(self, label="x", style=wx.RB_GROUP)
		x_slow.__label__ = "x_slow_axis"
		boxsizer.Add (x_slow, flag=wx.EXPAND, border=5)
		
		y_slow = wx.RadioButton(self, label="y")
		y_slow.__label__ = "y_slow_axis"
		boxsizer.Add (y_slow, flag=wx.EXPAND, border=5)
	
		sizer.Add (boxsizer, flag=wx.EXPAND, border=10) 
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		
		self.CreateSettingsDict()
		
########################################################################
#
#	Process performing scanning
#
########################################################################
class ProcRectangularScan (multiprocessing.Process) :
	"""
	Process which runes the scan.
	"""
	def __init__ (self, ScanningEvent, PauseScanningEvent, settings, 
			Spectrometer, MovingStageX, MovingStageY) :
		multiprocessing.Process.__init__(self)
		
		# Saving events controlling the scanning flow 
		self.ScanningEvent = ScanningEvent
		self.PauseScanningEvent = PauseScanningEvent
		# Scanning settings
		self.settings = settings
		# Saving device managers
		self.Spectrometer = Spectrometer
		self.MovingStageX = MovingStageX
		self.MovingStageY = MovingStageY
		
	def run (self) :
		"""
		Perform scanning
		"""
		# Generate axis
		x_axis = np.linspace( self.settings["x_min"], self.settings["x_max"], self.settings["x_points"])
		y_axis = np.linspace( self.settings["y_min"], self.settings["y_max"], self.settings["y_points"])
	
		# Cache the indecies of the coordinate
		indx_x = dict( (x, num) for num, x in enumerate(x_axis) )
		indx_y = dict( (y, num) for num, y in enumerate(y_axis) )
		
		# Generate rectangular mesh
		if self.settings["x_slow_axis"] : 
			Grid = itertools.product( x_axis, y_axis )
		if self.settings["y_slow_axis"] : 
			Grid = ( (x,y) for y, x in itertools.product( y_axis, x_axis ) )
		
		with h5py.File (self.settings["filename"], 'a') as hdf5_file :
			# Saving scanning parameters 
			ScaningParameters = hdf5_file.create_group("scanning_parameters") 
			ScaningParameters["wavelengths"] = self.Spectrometer.GetWavelengths()
			ScaningParameters["x_axis"] = x_axis
			ScaningParameters["y_axis"] = y_axis
			
			# Scanning and saving the results into the file
			SpectraGrope = hdf5_file.create_group("scaned_spectra") 
			
			# Progress info variables
			print "Start rectangular scanning" 
			initial_time 	= time.clock()
			completed  		= 0
			total_number = x_axis.size * y_axis.size
			
			for x, y in Grid :
				# pause scanning if user requested
				while self.PauseScanningEvent.is_set() : pass
			
				# check whether the user requested to quit
				if not self.ScanningEvent.is_set() : break
				
				# Move to a point
				self.MovingStageX.AbsMoveWait(x); self.MovingStageY.AbsMoveWait(y)
			
				# Saving the spectrum
				SpectraGrope[ "spectra_%d_%d" % (indx_x[x], indx_y[y]) ] = self.Spectrometer.AcquiredData()
		
				# Print progress info		
				completed += 1
				percentage_completed = 100.*completed / total_number
				seconds_left = (time.clock() - initial_time)*(100./percentage_completed - 1.)
				# convert to hours:min:sec
				m, s = divmod(seconds_left, 60)
				h, m = divmod(m, 60)
				print "%.2f %% completed. Time left: %d:%02d:%02d" % ( percentage_completed, h, m, s ) 

		
		# Signal that the scanning is over
		self.ScanningEvent.clear()
		
########################################################################
#
#	Pre-processing results of scanning
#
########################################################################

def PostProcessRectangularScan (filename) :
	"""
	Post process results of rectangular scanning:
	
	plot total interspecies of the scanned regions 
	"""
	import matplotlib.pyplot as plt
	from scipy.sparse import coo_matrix
	
	with h5py.File (filename, 'r') as hdf5_file :
		# Loading scanning settings
		ScaningParameters 	= hdf5_file["scanning_parameters"]
		x_axis = ScaningParameters["x_axis"][...]
		y_axis = ScaningParameters["y_axis"][...]
		
		# Getting total intensity data for creating sparse matrix
		data = (  (spectrum[...].sum(),) + tuple(map(int, key.split('_')[-2:]))
					for key,spectrum in hdf5_file["scaned_spectra"].items() )
		data, row, col = zip(*data)
			
	intensity = coo_matrix( (data, (row,col)) ).todense()
		
	plt.title("Total intensity (counts)")
	plt.imshow( intensity, origin='lower', interpolation='nearest', extent=[x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]] )
	plt.colorbar()
	plt.xlabel("x (mm)")
	plt.ylabel("y (mm)")
	plt.show()
		
if __name__ == '__main__' :
	app = wx.PySimpleApp()
	openFileDialog = wx.FileDialog (None, "Chose HDF5 file for viewing", "", "", \
		"HDF5 files (*.hdf5)|*.hdf5", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR)
	# Check whether user cancelled
	if openFileDialog.ShowModal() == wx.ID_CANCEL : raise ValueError ("User did not chose the file to view")
	del app
	
	# Run post processing
	PostProcessRectangularScan( openFileDialog.GetPath() )