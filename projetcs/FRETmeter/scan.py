"""
Program performing scanning
"""
import wx
import sys
import getopt

########################################################################

class FRETmeterScan (wx.Frame):
	"""
	FRETmeter scanner class
	"""
	def __init__(self, parent=None, **mysql_config):
		# save configurations needed to connect to the library
		self.mysql_config = mysql_config
	
		# Create GUI
		dw, dh = wx.DisplaySize()
		wx.Frame.__init__ (self, parent, title="FRETmeter",
								size=(0.2*dw, 0.5*dh) )
		
		self.ConstructGUI ()
		self.Center()
		#self.Maximize()
		self.Show ()
	
	def ConstructGUI (self) :
		""" Build GUI """
		self.panel = wx.Panel(self)
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		boxsizer = wx.StaticBoxSizer(wx.StaticBox(self.panel, label="My SQL settings"), wx.VERTICAL)
		
		#########################################################################################
		# My SQL settings
		#########################################################################################
		
		# Scan ID
		boxsizer.Add( wx.StaticText(self.panel, label="Scan ID"), flag=wx.LEFT, border=5 )
		scan_id = wx.SpinCtrl(self.panel, min=-1, max=1e9)
		scan_id.SetValue(int(self.mysql_config.get("scanid", -1)))
		boxsizer.Add(scan_id, flag=wx.EXPAND, border=5)
		
		# User 
		boxsizer.Add( wx.StaticText(self.panel, label="\nUser"), flag=wx.LEFT, border=5 )
		user_ctrl = wx.TextCtrl(self.panel)
		user_ctrl.SetValue(self.mysql_config.get("user", 'fretmeter'))
		boxsizer.Add(user_ctrl, flag=wx.EXPAND, border=5)
		
		# Password
		boxsizer.Add( wx.StaticText(self.panel, label="\nPassword"), flag=wx.LEFT, border=5 )
		password_ctrl = wx.TextCtrl(self.panel, style=wx.TE_PASSWORD)
		password_ctrl.SetValue(self.mysql_config.get("password", ""))
		boxsizer.Add(password_ctrl, flag=wx.EXPAND, border=5)
		
		# Host
		boxsizer.Add( wx.StaticText(self.panel, label="\nHost"), flag=wx.LEFT, border=5 )
		host_ctrl = wx.TextCtrl(self.panel)
		host_ctrl.SetValue(self.mysql_config.get("host", "localhost"))
		boxsizer.Add(host_ctrl, flag=wx.EXPAND, border=5)
		
		# Database
		boxsizer.Add( wx.StaticText(self.panel, label="\nDatabase"), flag=wx.LEFT, border=5 )
		database_ctrl = wx.TextCtrl(self.panel)
		database_ctrl.SetValue(self.mysql_config.get("database", ""))
		boxsizer.Add(database_ctrl, flag=wx.EXPAND, border=5)
		
		sizer.Add (boxsizer, flag=wx.EXPAND, border=10)
		
		#########################################################################################			
		self.panel.SetSizer(sizer)
		self.panel.SetAutoLayout(True)
		self.panel.Layout() 
		
if  __name__ == '__main__':
	# Extract the command line arguments
	mysql_config, _ = getopt.getopt(
		sys.argv[1:], "p:u:h:d:i:", 
		["password=", "user=", "host=", "database=", "scanid="]
	)
	
	# Starting the application
	wx_app = wx.App(False)
	FRETmeterScan(**dict(mysql_config))
	wx_app.MainLoop()