########################################################################
#
# Wrappers for  algorithms in DEAP library (http://deap.readthedocs.org/en/latest/api/algo.html)
#
########################################################################

import wx
from wx.lib.agw.floatspin import FloatSpin as wxFloatSpin

from libs.gui.hardware_control import HardwareGUIControl

########################################################################

import numpy as np
import random
import array
import cPickle as pickle

from deap import algorithms
from deap import base
from deap import creator
from deap import tools

from h5py._hl.group import Group as HD5Group

########################################################################

class GATab (HardwareGUIControl) :
	"""
	GUI Tab defying set of parameters common for DEAP optimization algorithms
	"""
	def __init__ (self, parent) :
		HardwareGUIControl.__init__(self, parent)
		
		self.sizer = wx.BoxSizer(wx.VERTICAL)
		
		sizer.Add (wx.StaticText(self, label="Gaussian mutation:\nMean"), flag=wx.LEFT, border=5)
		mu_ctrl = wxFloatSpin (self, min_val=0, max_val=1, increment=0.01, value=0.0, digits=3)
		mu_ctrl.__label__ = "gaussian_mutation_mu"
		sizer.Add (mu_ctrl , flag=wx.EXPAND, border=5)
		
		sizer.Add (wx.StaticText(self, label="Sigma"), flag=wx.LEFT, border=5)
		sigma_ctrl = wxFloatSpin (self, min_val=0, max_val=1, increment=0.01, value=0.05, digits=3)
		sigma_ctrl.__label__ = "gaussian_mutation_sigma"
		sizer.Add (sigma_ctrl , flag=wx.EXPAND, border=5)
	
		sizer.Add (wx.StaticText(self, label="Independent  probability for attribute mutation"), flag=wx.LEFT, border=5)
		indpb_ctrl = wxFloatSpin (self, min_val=0, max_val=1, increment=0.01, value=1., digits=3)
		indpb_ctrl.__label__ = "mutation_indpb"
		sizer.Add (indpb_ctrl, flag=wx.EXPAND, border=5)	
			
		# Separator
		sizer.Add (wx.StaticText(self), border=5)

		# cxpb - The probability of mating two individuals 
		sizer.Add (wx.StaticText(self, label="Mating probability"), flag=wx.LEFT, border=5)
		cxpb_ctrl = wxFloatSpin (self, min_val=0, max_val=1, increment=0.01, value=0.5, digits=3)
		cxpb_ctrl.__label__ = "cxpb"
		sizer.Add (cxpb_ctrl, flag=wx.EXPAND, border=5)
		
		# mutpb - The probability of mutating an individuals
		sizer.Add (wx.StaticText(self, label="Mutating probability"), flag=wx.LEFT, border=5)
		mutpb_ctrl = wxFloatSpin (self, min_val=0, max_val=1, increment=0.01, value=0.1, digits=3)
		mutpb_ctrl.__label__ = "mutpb"
		sizer.Add (mutpb_ctrl, flag=wx.EXPAND, border=5)
		
		self.SetSizer(self.sizer)
		
		############### GUI is created, now generate settings ######################
		self.CreateSettingsDict()


##########################################################################

from itertools import repeat
from collections import Sequence

def mutBoundedGaussian(individual, mu, sigma, indpb, min_val, max_val):
	"""
	This function applies an additive gaussian mutation of mean `mu` and standard
	deviation `sigma` on the input `individual`. 
	The resulted value of `individual` stays bounded between `min_val` and `max_val`.
	This function is derived from the function  `deap.tools.mutation.mutGaussian`
	"""
	size = len(individual)
	if not isinstance(mu, Sequence):
		mu = repeat(mu, size)
	elif len(mu) < size:
		raise IndexError("mu must be at least the size of individual: %d < %d" % (len(mu), size))
		
	if not isinstance(sigma, Sequence):
		sigma = repeat(sigma, size)
	elif len(sigma) < size:
		raise IndexError("sigma must be at least the size of individual: %d < %d" % (len(sigma), size))
		
	if not isinstance(min_val, Sequence):
		min_val = repeat(min_val, size)
	elif len(min_val) < size : 
		raise IndexError("min_val must be at least the size of individual: %d < %d" % (len(min_val), size))
	
	if not isinstance(max_val, Sequence):
		max_val = repeat(max_val, size)
	elif len(max_val) < size : 
		raise IndexError("max_val must be at least the size of individual: %d < %d" % (len(max_val), size))
	
	for i, m, s, lower, upper in zip(xrange(size), mu, sigma, min_val, max_val):
		if random.random() < indpb:
			individual[i] = ( individual[i] + random.gauss(m, s) )%(upper - lower) + lower
			#individual[i] = min(max(individual[i] + random.gauss(m, s), lower), upper)
	
	return individual,
	
##########################################################################

class PixelWiseGA :
	"""
	Using DEAP library, set up a GA where each individual is an array.
	This algorithm is well suited for phase and/or amplitude optimization(e.g., finding transform limited phase).
	"""
	def __init__ (self, num_pixels, ga_settings, FitnessFunction, checkpoint=None, min_val=0, max_val=1) :
		"""
		`num_pixels` -- number of pixels in an array -- an individual member of the population
		`ga_settings` -- GA settings, as given by the tab
		`FitnessFunction` -- fitness function to be maximized
		`checkpoint` -- checkpoint to restart the optimization
		`min_val` and `max_val` define bounds for the values in the array-individual
		"""
		creator.create("FitnessMax", base.Fitness, weights=(1.0,))
		creator.create("Individual", array.array, typecode='f', fitness=creator.FitnessMax)

		self.optimization_toolbox = base.Toolbox()
			
		# Attribute generator
		self.optimization_toolbox.register("attr_float", random.rand)
						#lambda mu, sigma : random.gauss(mu,sigma) % (max_val - min_val) + min_val, 
						# lambda mu, sigma : min(max(random.gauss(mu,sigma),min_val), max_val), 
						#ga_settings["gaussian_mutation_mu"], ga_settings["gaussian_mutation_sigma"] )
		
		# Structure initializers
		self.optimization_toolbox.register("individual", 
					tools.initRepeat, creator.Individual, self.optimization_toolbox.attr_float, num_pixels)
		self.optimization_toolbox.register("population", tools.initRepeat, list, self.optimization_toolbox.individual)
		self.optimization_toolbox.register("evaluate", FitnessFunction)
		self.optimization_toolbox.register("mate", tools.cxTwoPoint)
		
		self.optimization_toolbox.register("mutate", mutBoundedGaussian, mu=ga_settings["gaussian_mutation_mu"], 
			sigma=ga_settings["gaussian_mutation_sigma"], indpb=ga_settings["mutation_indpb"],  min_val=min_val, max_val=max_val)
		self.optimization_toolbox.register("select", tools.selTournament, tournsize=ga_settings["selection_size"])
		
		# Statistics
		self.optimization_stats = tools.Statistics(lambda ind: ind.fitness.values)
		self.optimization_stats.register("avg", np.mean)
		self.optimization_stats.register("std", np.std)
		self.optimization_stats.register("min", np.min)
		self.optimization_stats.register("max", np.max)
		
		# List where log will be saved
		self.optimization_log = {"min" : [], "avg" : [], "max" : [] }
		
		if checkpoint :
			################# Load last iteration of optimization ##################
			if isinstance(checkpoint, HD5Group) :
				self.optimization_pop 		= pickle.loads( str(checkpoint["population"][...]) )
				self.optimization_hof 		= pickle.loads( str(checkpoint["halloffame"][...]) )
			else :
				self.optimization_pop 		= pickle.loads( checkpoint["population"] )
				self.optimization_hof 		= pickle.loads( checkpoint["halloffame"] )
		else : 
			######################### Start new evolution ######################
			self.optimization_pop = self.optimization_toolbox.population(n=ga_settings["population_size"])
			self.optimization_hof = tools.HallOfFame( ga_settings["halloffame_size"] )
		
		# Select the optimization algorithm and its arguments
		self.RunOptimizationIteration  = algorithms.eaSimple
		self.opt_args = (self.optimization_pop, self.optimization_toolbox)
		self.opt_kwargs = dict(cxpb=ga_settings["cxpb"], mutpb=ga_settings["mutpb"], 
					ngen=1, stats=self.optimization_stats, halloffame=self.optimization_hof, verbose=True)
	
	def NextIteration (self) :
		"""
		Perform a GA iteration
		"""
		# Perform an iteration of optimization
		self.optimization_pop, log = self.RunOptimizationIteration (*self.opt_args, **self.opt_kwargs)
		# Saving log
		log = log.pop()
		for key, value in self.optimization_log.items() : value.append( log[key] )
	
	def GetOptimizationLog (self) :
		return self.optimization_log
	
	def GetCheckPoint (self) :
		"""
		pickle the current iteration of the optimization 
		"""
		return {"population" : np.string_(pickle.dumps(self.optimization_pop)), 
				"halloffame" : np.string_(pickle.dumps(self.optimization_hof)) }

##########################################################################