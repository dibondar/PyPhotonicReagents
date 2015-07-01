"""
GUI for setting up scanning 
"""
import wx
import visvis
import functools
import cPickle as pickle

from libs.gui.hardware_control import HardwareGUIControl

class SetScanTab (HardwareGUIControl) :
	"""
	This class implements tab for setting microscope scan
	"""
	def __init__(self, parent) :
		HardwareGUIControl.__init__(self, parent)
		
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		dw, dh = wx.DisplaySize()
		############################################################################
		#
		#	Geometry settings
		#
		############################################################################
		
		boxsizer = wx.StaticBoxSizer(wx.StaticBox(self, label="Scanning geometry"), wx.VERTICAL)	
	
		# Upper layer of points
		boxsizer.Add (wx.StaticText(self, label="Upper layer of points"), flag=wx.LEFT, border=5)
		self.upper_points_ctrl = wx.TextCtrl(self, style=wx.TE_MULTILINE, size=(0.1*dw, 0.15*dh))
		self.upper_points_ctrl.__label__ = "upper_points"
		boxsizer.Add(self.upper_points_ctrl)
		
		# Button to add point to upper layer
		add_upper_point_btn = wx.Button(self, label="Add point to upper layer")
		add_upper_point_btn.Bind( wx.EVT_BUTTON, 
					functools.partial(self.AddPoint, points_ctrl=self.upper_points_ctrl)
				)
		boxsizer.Add(add_upper_point_btn, flag=wx.EXPAND, border=5)
		
		# Lower layer of points
		boxsizer.Add (wx.StaticText(self, label="\nLower layer of points"), flag=wx.LEFT, border=5)
		self.lower_points_ctrl = wx.TextCtrl(self, style=wx.TE_MULTILINE, size=(0.1*dw, 0.15*dh))
		self.lower_points_ctrl.__label__ = "lower_points"
		boxsizer.Add(self.lower_points_ctrl)
		
		# Button to add point to lower layer
		add_lower_point_btn = wx.Button(self, label="Add point to lower layer")
		add_lower_point_btn.Bind( wx.EVT_BUTTON, 
					functools.partial(self.AddPoint, points_ctrl=self.lower_points_ctrl)
				)
		boxsizer.Add(add_lower_point_btn, flag=wx.EXPAND, border=5)
		
		# Initialize container for points and corresponding photos
		self.points = {}
		
		# Spacer
		boxsizer.Add (wx.StaticText(self, label=""), flag=wx.LEFT, border=5)
		
		# Delete all points button
		clear_all_points_btn = wx.Button(self, label="Delete all points")
		clear_all_points_btn.Bind (wx.EVT_BUTTON, self.ClearAllPoints)
		boxsizer.Add (clear_all_points_btn, flag=wx.EXPAND, border=5)
		
		sizer.Add (boxsizer, flag=wx.EXPAND, border=10)
		
		############################################################################
		#
		#	MySQL settings
		#
		############################################################################
		
		# Spacer
		sizer.Add (wx.StaticText(self, label=""), flag=wx.LEFT, border=5)
		
		boxsizer = wx.StaticBoxSizer(wx.StaticBox(self, label="My SQL settings"), wx.VERTICAL)
		
		# User 
		boxsizer.Add( wx.StaticText(self, label="User"), flag=wx.LEFT, border=5 )
		user_ctrl = wx.TextCtrl(self, value="fretmeter")
		user_ctrl.__label__ = "user"
		boxsizer.Add(user_ctrl, flag=wx.EXPAND, border=5)
		
		# Password
		boxsizer.Add( wx.StaticText(self, label="\nPassword"), flag=wx.LEFT, border=5 )
		password_ctrl = wx.TextCtrl(self, style=wx.TE_PASSWORD)
		password_ctrl.__label__ = "password"
		boxsizer.Add(password_ctrl, flag=wx.EXPAND, border=5)
		
		# Host
		boxsizer.Add( wx.StaticText(self, label="\nHost"), flag=wx.LEFT, border=5 )
		host_ctrl = wx.TextCtrl(self, value='localhost')
		host_ctrl.__label__ = "host"
		boxsizer.Add(host_ctrl, flag=wx.EXPAND, border=5)
		
		# Database
		boxsizer.Add( wx.StaticText(self, label="\nDatabase"), flag=wx.LEFT, border=5 )
		database_ctrl = wx.TextCtrl(self)
		database_ctrl.__label__ = "database"
		boxsizer.Add(database_ctrl, flag=wx.EXPAND, border=5)
		
		sizer.Add (boxsizer, flag=wx.EXPAND, border=10)
		
		# Spacer
		sizer.Add (wx.StaticText(self, label=""), flag=wx.LEFT, border=5)
		
		# Button to start scanning
		start_btn = wx.Button(self, label="Start scanning")
		start_btn.Bind(wx.EVT_BUTTON, self.StartScan)
		sizer.Add(start_btn, flag=wx.EXPAND, border=10)
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		
		self.CreateSettingsDict()
		
	def AddPoint (self, event, points_ctrl) :
		"""
		Adding points to both upper and lower levels
		"""
		# Check whether the devices are initialized 
		try :
			self.MovingStage
			self.Camera
		except AttributeError :
			# Initialize moving stage
			self.MovingStage = self.parent.StartDev("MovingStage")
			# Initialize camera
			#self.Camera	= self.parent.StartDev("Camera")
		
		# retrieve current position
		position = self.MovingStage.GetCurrentPosition()
		
		# Add point 
		points_ctrl.AppendText( "%s ,\n" % str(position) ) 
		
		# save the photo
		self.points[position] = self.Camera.GetImage()
		
	def ClearAllPoints (self, event=None) :
		"""
		Remove all information about points
		"""
		self.points = {}
		self.upper_points_ctrl.Clear()
		self.lower_points_ctrl.Clear()
		
	def StartScan(self, event) :
		"""
		Initialize scanning
		"""
		
		# Start connecting to data base
		import mysql.connector
		from mysql.connector import errorcode
		
		#########################################################
		#
		#	Committing settings
		#
		#########################################################
		
		# Get all configurations
		settings = self.GetSettings()
		
		# Get configurations for connecting to MySQL
		# Note database is not specified
		mysql_config = dict( 
			(k,  settings[k]) for k in ["user", "password", "host"]
		)
		
		# Establishing connection with the server
		cnx = mysql.connector.connect( **mysql_config )
		cursor = cnx.cursor()
		
		# Create the database if it does not exist
		cursor.execute(
			"CREATE DATABASE IF NOT EXISTS {}".format(settings["database"])
		)
		cnx.database = settings["database"]
		
		# Delete the table
		#cursor.execute("DROP TABLE scan_logs")
		
		# Create table where all settings are stored
		cursor.execute(
			"""
			CREATE TABLE IF NOT EXISTS scan_logs (
				id INT(11) NOT NULL PRIMARY KEY AUTO_INCREMENT,
				date_time DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP,
				settings BLOB NOT NULL  
			) Engine= InnoDB
			"""
		)
		
		# Save settings for current scanning section
		cursor.execute(
			"""
			INSERT INTO scan_logs
			(settings)
			VALUES (%s)
			""",
			( mysql.connector.Binary( 
				pickle.dumps(self.parent.GetAllSettings()) 
			), )
		)
		
		# Obtain the current scanning id
		cursor.execute("SELECT MAX(id) FROM scan_logs")
		scan_id = next(cursor)[0]
		print "\nCurrent scanning id is %d" % scan_id
	
		#########################################################
		#
		#	Committing pictures taken
		#
		#########################################################
	
		# Delete the table
		#cursor.execute("DROP TABLE bright_field_images")
		
		# Create table storing bright-field images
		cursor.execute(
			"""
			CREATE TABLE IF NOT EXISTS bright_field_images (
				id INT(11) NOT NULL PRIMARY KEY,
				x DOUBLE NOT NULL,
				y DOUBLE NOT NULL,
				z DOUBLE NOT NULL,
				image MEDIUMBLOB NOT NULL 
			) Engine= InnoDB
			"""
		)
		
		# Inserting each photo
		add_photo = """
			INSERT INTO bright_field_images
			(id, x, y, z, image)
			VALUES ( %(id)s, %(x)s, %(y)s, %(z)s, %(image)s )
		"""		
		for position, img in self.points.items() :
			photo_data = {
				"id"	: scan_id,
				"x" 	: position[0],
				"y"		: position[1],
				"z"		: position[2],
				"image" : mysql.connector.Binary(img.dumps())
			}
			cursor.execute(add_photo, photo_data)
		
		# Make sure data is committed to the database
		cnx.commit()

		cursor.close()
		cnx.close()
		
		mysql_config["scanid"] = int(scan_id)

		# Start scanning program with arguments	
		import subprocess
		subprocess.Popen(
			['python.exe', 'scan.py']
			+ ["%s=%s" % kv for kv in mysql_config.iteritems()]
		)
		