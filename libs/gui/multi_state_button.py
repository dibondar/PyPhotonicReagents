import wx
import itertools

class MultiStateButton (wx.Button) :
	"""
	Button that has multiple states
	"""
	def __init__ (self,  parent, actions=[], labels=[], colours=[], **kwargs) :
		"""
		Constructor:
		
			actions -- list of methods corresponding to pressing the button (to peform at different state)
			labels	-- list of names of different states
			colours -- lost of colours to shade the button at different stages
		"""
		wx.Button.__init__ (self, parent, **kwargs)
		
		if len(actions) == 0 :
			print Warning("No actions are defined. MultiStateButton will behave as an ordinary button")
			return
		
		# generate default colour scheme
		if len(colours) == 0 :
			colours = ['']
			
		if len(colours) < len(actions) :
			colours.extend( 
				itertools.islice( itertools.cycle(['red', 'blue', 'green']), len(actions)-len(colours) )
			)
		else :
			colours = colours[:len(actions)]
		
		if len(labels) < len(actions) :
			labels.extend( itertools.repeat('', len(actions)-len(labels)) )
		else :
			labels = labels[:len(actions)]
			
		# Create iterator
		self._actions = itertools.cycle( itertools.izip(actions, labels, colours) )
		
		self.SetNextStep()
		
		# Bind  
		self.Bind( wx.EVT_BUTTON, self.NextAction )
		
	def SetNextStep (self) :
		"""
		Prepare the button for the next action 
		"""
		self.CurrentAction, label, colour = self._actions.next()
		self.SetLabel(label)
		self.SetBackgroundColour(colour)
	
	def NextAction(self, event=None) :
		"""
		Button was clicked execute the current action
		"""
		self.CurrentAction(event) 
		self.SetNextStep()