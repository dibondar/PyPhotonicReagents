import wx
import itertools

class MultiStateButton (wx.Button) :
	"""
	Button that has multiple states
	"""
	def __init__ (self,  parent, states=[], **kwargs) :
		"""
		Constructor:
			state_specs -- list of tuples specifying state, e.g, [(callable function, label, colour)]
		"""
		wx.Button.__init__ (self, parent, **kwargs)
		
		if len(states) == 0 :
			print ("Warning: No actions are defined. MultiStateButton will behave as an ordinary button")
			return
		
		# generate default colour scheme
		colours = itertools.chain( [''], itertools.cycle(['', 'red', 'blue', 'green']) )
		
		actions = []
		
		for specs in states :
			assert len(specs) 
			
			if len(specs) == 1 :
				actions.append( (specs[0], '', colours.next()) )
			elif len(specs) == 2 :
				actions.append( (specs[0], specs[1], colours.next()) )
			else :
				actions.append( tuple(specs[:3]) )
			
		# Create iterator
		self._actions = itertools.cycle(actions)
		
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