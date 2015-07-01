"""
Program for setting up the scan 
"""
# Add PyPhotonicReagents directory to enable imports  
if __name__ == '__main__' :
	import os
	os.sys.path.append(os.path.abspath('..\..'))
########################################################################

import wx

# Real time plotting
import visvis

# GUI components
#from PyPhotonicReagents.libs.gui.info_tab import InfoTab
from set_scan_tab  import SetScanTab
from libs.gui.info_tab import InfoTab
from libs.gui.settings_notebook import SettingsNotebook

# Hardware
#from PyPhotonicReagents.libs.dev.gui.pico_harp import ManagerPicoHarp, PicoHarpTab
from libs.dev.pico_harp import PicoHarpTab
from libs.dev.thorlabs_camera import ThorlabsCameraTab
from libs.dev.newport_moving_stage import NewportMovingStageTab
		
########################################################################

class SetFRETmeterScan (wx.Frame): 
	"""
	GUI for setting up FRET meter scan
	"""
	def __init__ (self, parent) :
		# Create GUI
		dw, dh = wx.DisplaySize()
		wx.Frame.__init__ (self, parent, title="Set up FRETmeter scan",
								size=(0.9*dw, 0.88*dh) )
		
		self.ConstructGUI ()
		self.Center()
		#self.Maximize()
		self.Show ()
	
	def ConstructGUI (self) :
		""" Build GUI """
		self.panel = wx.Panel(self)
		sizer = wx.GridBagSizer ()
		
		############################ Settings Notebook ############################
		self.SettingsNotebook = SettingsNotebook(self.panel, specs=[
				(SetScanTab, "SetScan", None),	(InfoTab, "Info", None),
				(PicoHarpTab, "PicoHarp"), (ThorlabsCameraTab, "Camera"), 
				(NewportMovingStageTab,	"MovingStage")
			]
		)
		
		sizer.Add(self.SettingsNotebook, pos=(0, 0), span=(1, 1), flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT , border=10)

		# Make sure all devices are closed when exiting
		self.Bind(wx.EVT_CLOSE, self.SettingsNotebook.StopAllDevices)	
		
		############################ Command panel ############################
		boxsizer = wx.BoxSizer (wx.VERTICAL)
		
		# Separator
		boxsizer.Add (wx.StaticText(self.panel), flag=wx.EXPAND, border=5)
	
		
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
		
#########################################################################

if __name__ == '__main__' :
	app = visvis.use('wx')
	app.Create()
	SetFRETmeterScan (None)
	app.Run()