########################################################################
#
#	This file contains class for saving additional description
#
########################################################################

import wx
from libs.gui.hardware_control import HardwareGUIControl

class InfoTab (HardwareGUIControl) :
	"""
	This class represents a GUI for saving descriptions, 
	and it is not connected to any hardware
	"""
	def __init__(self, parent=None, dev=None) :
		HardwareGUIControl.__init__(self, parent, dev)
		
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		sizer.Add (wx.StaticText(self, label="Description"), flag=wx.LEFT, border=5)
		
		# Text entry field
		dw, dh = wx.DisplaySize()
		info_txt = wx.TextCtrl (self, style=wx.TE_MULTILINE|wx.TE_RICH2|wx.TE_AUTO_URL, size=(-1, 0.4*dh))
		info_txt.__label__ = "info_txt"
		sizer.Add (info_txt, flag=wx.EXPAND, border=5)
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		self.CreateSettingsDict()