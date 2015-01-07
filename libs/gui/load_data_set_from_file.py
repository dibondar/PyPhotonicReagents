
import wx
import h5py
import ntpath
import os
from collections import defaultdict

from h5py._hl.group import Group as HD5Group

##########################################################################################
#
#	Functions for working with tree represented as list of list
#
##########################################################################################

def BuildTree (data_set, tree=[]) :
	"""
	Built Tree representation (as a nested list of list) of the HDF5 data set `data_set` recursively
	"""
	# Looping over the content of `data_set`
	for key, item in data_set.iteritems() :
		# Create new node
		if isinstance(item, HD5Group) :
			# Since this item is a grope, then recursively load  
			tree.append( (key, []) )
			BuildTree(item, tree[-1][-1] )
		else :
			# Since this item is a data set, then assign data reference 
			# that allows to load this data set	
			tree.append( key )
	return tree
		
def SortTree (tree, reverse=True) :
	"""
	Recursively sort the tree assuming that the nodes are numbers otherwise perform lexographic sort 
	"""
	key_func = lambda x : ( int(x[0]) if isinstance(x, tuple) else int(x) )
	try : 
		tree.sort(key=key_func, reverse=reverse)
	except ValueError : 
		tree.sort(reverse=reverse)
			
	for node in tree :
		if isinstance(node, tuple) : 
			SortTree( node[-1], reverse=reverse )
				
########################################################################################## 

class LoadDataSetDialog (wx.Dialog):
	"""
	This is a general dialog for loading data sets from HDF5 files. 
	"""
	def __init__(self, key=None, key_options=None, **kwargs):
		
		self.key = key
		self.key_options = key_options
		
		wx.Dialog.__init__(self, style=wx.MAXIMIZE_BOX|wx.MINIMIZE_BOX|wx.RESIZE_BORDER|wx.DEFAULT_DIALOG_STYLE, **kwargs)
		
		self.CreateGUI()
	
		# HDF5 data set from where data sets will be loaded
		self.data_to_load = {}
		self.loaded_data = []
		
	def CreateGUI(self):
		panel = wx.Panel(self)
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		######################### Adding control buttons #########################
		hsizer = wx.BoxSizer(wx.HORIZONTAL)
		
		# Load file button
		load_file_ctrl = wx.Button(panel, label="Load file...") 
		load_file_ctrl.Bind (wx.EVT_BUTTON, self.LoadFile)
		hsizer.Add (load_file_ctrl, flag=wx.CENTER, border=5)
		
		# Delete item button
		del_item_ctrl = wx.Button(panel, label="Remove") 
		del_item_ctrl.Bind (wx.EVT_BUTTON, self.DelItem)
		hsizer.Add (del_item_ctrl, flag=wx.CENTER, border=5)
		
		# Load datasets button
		load_ctrl = wx.Button(panel, label="Load data sets") 
		load_ctrl.Bind (wx.EVT_BUTTON, self.LoadDataSets)
		hsizer.Add (load_ctrl, flag=wx.CENTER, border=5)
		# This button will be treated as OK and CANCEL button of the Dialogue 
		#self.SetAffirmativeId( load_ctrl.GetId() )
		#self.SetEscapeId( load_ctrl.GetId() )
		
		sizer.Add(hsizer, flag=wx.EXPAND|wx.TOP|wx.LEFT|wx.RIGHT|wx.GROW, border=5)
		
		######################### Adding control buttons #########################
		hsizer = wx.BoxSizer(wx.HORIZONTAL)
		
		display_width, display_hight = wx.DisplaySize()
		size=(0.5*display_width, 0.8*display_hight)
		
		# Tree structure 
		self.DataTree = wx.TreeCtrl (panel, size=size)
		self.DataTree.Bind (wx.EVT_TREE_ITEM_ACTIVATED, self.AddDataSet)
		hsizer.Add (self.DataTree, flag=wx.EXPAND|wx.LEFT|wx.GROW, border=5)
		
		# List structure
		self.DataList = wx.ListBox (panel, size=size, style=wx.LB_MULTIPLE)
		self.DataList.Bind (wx.EVT_KEY_DOWN, self.OnListKey) 
		hsizer.Add (self.DataList, flag=wx.EXPAND|wx.RIGHT|wx.GROW, border=5)
		
		sizer.Add(hsizer, flag=wx.EXPAND|wx.CENTER|wx.GROW, border=5)
		###################################################################
	
		panel.SetSizer(sizer)
		sizer.Fit(panel)
		self.Maximize()
		self.Show()
	
	def LoadFile (self, event) :
		"""
		Loading HDF5 file
		"""
		openFileDialog = wx.FileDialog(self, "Open HDF5 file", "", "",
                                       "HDF5 files (*.hdf5)|*.hdf5", wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_CHANGE_DIR)
		# Check whether user cancelled
		if openFileDialog.ShowModal() == wx.ID_CANCEL: return 
	
		# Save current filename
		self.filename = os.path.abspath(openFileDialog.GetPath()) 
		
		#  Recursively construct tree representation of the HDF file 
		with h5py.File (self.filename, 'r') as F :
			tree = BuildTree( F )
			SortTree(tree)
			
		# Construct wx GUI Tree by using 
		self.DataTree.DeleteAllItems()
		root = self.DataTree.AddRoot(ntpath.basename(self.filename))
		self.BuildGUITree(tree, root)
		self.DataTree.Expand(root)
	
	def BuildGUITree (self, tree, parent_id,  abs_path='') :
		"""
		Initialize wx.TreeCtrl with `tree` recursively
		"""
		for node in tree :
			# Create new node
			if isinstance(node, tuple) : 	
				 # Since this item is a grope, then recursively load 
				key = node[0]
				self.BuildGUITree( node[-1], self.DataTree.AppendItem(parent_id, key), "%s/%s"%(abs_path, key) )
			else :
				# Since this item is a data set, then assign data reference 
				# that allows to load this data set	
				self.DataTree.AppendItem( parent_id, node, data=wx.TreeItemData("%s/%s"%(abs_path, node)) )
			
	def AddDataSet (self, event) :
		"""
		Add double clicked data set from tree item into the list
		"""
		data = self.DataTree.GetItemPyData( event.GetItem() )
		# Add to list
		if data is not None :
			# Add item into list
			label = ntpath.basename(self.filename) + data
			self.DataList.Append( label )
			# Save info for loading
			self.data_to_load[label] = (self.filename, data)
		
	def DelItem (self, event) :
		"""
		Delete selected items from the list
		"""
		selections = [ _ for _ in self.DataList.GetSelections() if _ != -1 ]
		
		# Delete selected items 
		for sel in selections :
			del self.data_to_load[ self.DataList.GetString(sel) ]
			self.DataList.Delete(sel)
			
	def OnListKey (self, event) :
		"""
		Key down event for self.DataList
		"""
		if event.GetKeyCode() == wx.WXK_DELETE :
			# Del key was pressed
			self.DelItem(event)
			
	def LoadDataSets (self, event) :
		"""
		Load data sets as specified by user
		"""
		# Sort data sets by file name
		data_groped = defaultdict(list)
		for file_name, data_name in self.data_to_load.itervalues() :
			data_groped[ file_name ].append( data_name )
			
		def FindValuesForKey (key, F, values=[]) :
			"""
			Find `key` recursively in file HDF5 file `F`.
			Values will be saved in list `values` and returned
			"""
			for k, item in F.iteritems() :
				if isinstance(item, HD5Group) :
					FindValuesForKey(key, item, values) 
				elif key in k :
					values.append( item[...] )
			return values
	
		# Loading data 
		self.loaded_data = []
		for file_name, list_data_name in data_groped.iteritems() :
			# Read specified data sets from file
			with h5py.File (file_name, 'r') as F :
				if self.key is None :
					# No labelling is required just load data
					self.loaded_data.extend( ( F[_][...] for _ in list_data_name ) )
				else :
					# Find value of key labelling the loaded data sets
					option = FindValuesForKey(self.key, F) 
					if len(option) == 1 :
						# Unique value is found 
						option = option.pop()
					else :
						# We need to ask user about that
						dlg = wx.SingleChoiceDialog (self,"Select pulse shaping options", 
								"What kind of pulse shapes are saved in file %s ?" % file_name, self.key_options )
						if dlg.ShowModal() <> wx.ID_OK : raise IOError 
						option = self.key_options[ dlg.GetSelection() ]
						
					# Marked load data set with found option
					self.loaded_data.extend( ( (option, F[_][...]) for _ in list_data_name ) )
					
		# Data loaded, close the window
		self.Close()
		
	def GetLoadedData(self) :
		"""
		This function returns loaded data
		"""
		return self.loaded_data 

##########################################################################################

class LoadPulseShapesDialog (LoadDataSetDialog):
	"""
	This is a dialogue for loading pulse shapes from HDF5 files. 
	"""
	def __init__(self, **kwargs):
	
		LoadDataSetDialog.__init__(self, "pulse_shaping_option", 
			("amplitude only", "phase only", "amplitude and phase"), **kwargs)
		

##########################################################################################
#
# Some tests
#
##########################################################################################

if __name__ == "__main__":
	app = wx.App(False)
	
	print "Tests\n"
	
	# Load pulse shapes
	dlg = LoadPulseShapesDialog(parent=None, title="Loading pulse shapes from HDF5 files")
	dlg.ShowModal()
	print "Loaded pulse shapes data:"
	print dlg.GetLoadedData()
	
	# Load anything
	dlg = LoadDataSetDialog(parent=None, title="Loading arbitrary data sets")
	dlg.ShowModal()
	print "\nLoaded unmarked data:"
	print dlg.GetLoadedData()
	
	app.MainLoop()