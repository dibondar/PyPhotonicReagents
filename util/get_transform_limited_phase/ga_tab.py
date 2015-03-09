########################################################################
#
#	GA for optimal dynamic discrimination experiment
# 
########################################################################

import wx
from wx.lib.agw.floatspin import FloatSpin as wxFloatSpin
import h5py
from h5py._hl.group import Group as HD5Group

import numpy as np
from scipy.optimize import nnls
from scipy.ndimage.filters import gaussian_filter
import random
import multiprocessing
from functools import partial
from itertools import repeat, combinations, izip, izip_longest, chain, islice
from collections import Counter, Sequence
import array
import cPickle as pickle
from operator import itemgetter 

from deap import algorithms
from deap import base
from deap import creator
from deap import tools

import visvis 

from libs.gui.hardware_control import HardwareGUIControl
from libs.gui.basic_window import SaveSettings 
from libs.gui.load_data_set_from_file import LoadPulseShapesDialog

from libs.dev.consts import * 

########################################################################

##########################################################################
#
#	Axillary functions
#
##########################################################################

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
			# individual[i] = min(max( individual[i] + random.gauss(m, s), lower), upper)
	
	return individual,
	
	
########################################################################
		
class _SpectrumPostProcess :
	"""
	Container for prepossessing spectrum
	"""
	################ Methods of post processing spectra ################
	@classmethod
	def VertBinned (cls, spectrum) :
		if len(spectrum.shape) == 2 : return spectrum.sum(axis=0)
		else : return spectrum
	
	@classmethod
	def Smooth (cls, spectrum) :
		return gaussian_filter(spectrum, sigma=3)
	
	@classmethod
	def Sum (cls, spectrum) :
		return np.array( [spectrum.sum()] )
	
	@classmethod
	def VertBinnedSmooth (cls, spectrum) :
		return cls.Smooth( cls.VertBinned(spectrum) )
	
	@classmethod
	def SmoothSum (cls, spectrum) :
		return cls.Sum( cls.VertBinnedSmooth(spectrum) )
	
	################ Service implementations ################
	@classmethod
	def _GetDict (cls) :
		return { 	"vertical binned spectrum" 			: cls.VertBinned,
					"smooth vertical binned spectrum" 	: cls.VertBinnedSmooth,
					"total intensity" 					: cls.Sum,
					"smooth total intensity"			: cls.SmoothSum }
	
	@classmethod
	def Select (cls, key) :
		return cls._GetDict()[key]
	
	@classmethod
	def GetChoices (cls) :
		return cls._GetDict().keys()

########################################################################

class GA_Tab (HardwareGUIControl) :
	"""
	GUI for GA algorithm 
	"""
	def __init__ (self, parent) :
	
		HardwareGUIControl.__init__(self, parent)
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		# Separator
		sizer.Add (wx.StaticText(self), border=5)
		
		# Pulse shaping options
		sizer.Add ( wx.StaticText(self, label="Pulse shaping options") )
		pulse_shaping_options = {
			"amplitude only"		: 
				( 1, lambda ind : ( ind, np.zeros(len(ind)) ) ),
			"phase only"			: 
				( 1, lambda ind : ( np.ones(len(ind)), ind ) ),
			"amplitude and phase"	: 
				( 2, lambda ind : ( ind[:len(ind)/2], ind[len(ind)/2:] ) )
		} 
		self.Ind_length, self.ind2pulse_shape = zip(*pulse_shaping_options.values())
		# Multiplicative factor to determine the number of optimization variables  
		# controls the size of a single individual, which is equal to ind_len*self.num_pixels
		self.Ind_length 		= dict( zip(pulse_shaping_options.keys(), self.Ind_length) )
		# Function converting GA individual to pulse shape 
		self.ind2pulse_shape	= dict( zip(pulse_shaping_options.keys(), self.ind2pulse_shape) )
		
		choices = pulse_shaping_options.keys()
		pulse_shaping_ctrl = wx.ComboBox (self, choices=choices, value=choices[0], style=wx.CB_READONLY )
		pulse_shaping_ctrl.__label__ = "pulse_shaping_option"
		sizer.Add (pulse_shaping_ctrl,  flag=wx.EXPAND, border=5)	
		
		# Spectrum post-processing options
		sizer.Add ( wx.StaticText(self, label="Spectrum post-processing options") )
		choices = _SpectrumPostProcess.GetChoices()
		spectrum_options_ctrl = wx.ComboBox (self, choices=choices, value=choices[0], style=wx.CB_READONLY )
		spectrum_options_ctrl.__label__ = "spectrum_postprocess"
		sizer.Add (spectrum_options_ctrl,  flag=wx.EXPAND, border=5)	
		
		
		# Separator
		#sizer.Add (wx.StaticText(self), border=5)
		
		
		# Size of pixel bundle
		sizer.Add (wx.StaticText(self, label="\nPixels to bundle"), flag=wx.EXPAND, border=5)
		pixel_bundle_width = wx.SpinCtrl(self, value="1", min=1, max=640)
		pixel_bundle_width.__label__ = "pixel_bundle_width"
		sizer.Add (pixel_bundle_width, flag=wx.EXPAND, border=5)
		
		# Population size
		sizer.Add (wx.StaticText(self, label="Population size"), flag=wx.LEFT, border=5)
		population_size_ctrl  = wx.SpinCtrl (self, value="10", min=1, max=100000)
		population_size_ctrl.__label__ = "population_size"
		sizer.Add (population_size_ctrl , flag=wx.EXPAND, border=5)
		
		# Separator
		sizer.Add (wx.StaticText(self), border=5)
			
		########################################################################
		# Setting characterizing mutation
		########################################################################
	
		# sizer.Add (wx.StaticText(self, label="Gaussian mutation:"), flag=wx.LEFT, border=5)
		
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
		
		# Separator
		sizer.Add (wx.StaticText(self), border=5)
		
		# Record background signal
		background_signal_button = wx.Button (self, label="Record background")
		background_signal_button.Bind ( wx.EVT_BUTTON, self.RecordBackground )
		sizer.Add (background_signal_button, flag=wx.EXPAND, border=5)
		
		
		# Coherently suppress fluorescence button
		self.max_shg_button = wx.Button (self)
		self.max_shg_button._start_label = "Start GA to maximize SHG" 
		self.max_shg_button._start_method = partial(self.StartOptimization, DoGA=self.DoSingleIterGA_SHGMaximize)
		self.max_shg_button._stop_label = "STOP GA to maximize SHG"
		self.max_shg_button.SetLabel (self.max_shg_button._start_label)
		self.max_shg_button.Bind( wx.EVT_BUTTON, self.max_shg_button._start_method )
		sizer.Add(self.max_shg_button, flag=wx.EXPAND, border=5)
		
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		self.CreateSettingsDict()

	def RecordBackground (self, event=None) :
		"""
		Record background spectrum
		"""
		# Create pseudonyms 
		self.DevSpectrometer 	= self.parent.Spectrometer.dev
		
		# Initiate spectrometer
		settings = self.parent.Spectrometer.GetSettings()
		if self.DevSpectrometer.SetSettings(settings) == RETURN_FAIL : return
		
		# Record background 
		self.background_signal = self.DevSpectrometer.AcquiredData().astype(np.int)
	
	def RecordReferenceSignal (self) :
		"""
		Record fluorescence from a reference pulse
		"""
		# sent reference mask
		ref_pulse = np.ones(self.num_params)
		self.DevPulseShaper.SetAmplPhase( *self.Ind2PulseShape(ref_pulse) )
		
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return
		
		# Get spectrum
		spectrum = self.GetSampleSpectrum()
		self.reference_signal = spectrum
			
		try :
			self.log_reference_signal.append( spectrum.sum() )
		except KeyError :
			self.log_reference_signal = [ spectrum.sum() ]
			
	def CheckBackground (self) :
		"""
		Check whether the background signal is recorded and ready to be used. 
		"""
		try : 
			self.background_signal 
		except AttributeError :
			def SetBackgroundZero () :
				self.background_signal = 0
				
			options = { "record background now" : self.RecordBackground, 
						"continue optimization without recording background" : SetBackgroundZero }
						
			dlg = wx.SingleChoiceDialog (self, 'Background sygnal has not been recoreded. Select one of the following option', 
				'Background signal not found', options.keys(), wx.CHOICEDLG_STYLE ) 
			
			if dlg.ShowModal() == wx.ID_OK :
				options[ dlg.GetStringSelection() ]()
			else :
				# user cancel
				return
	
	def ResetLogs (self) : 
		"""
		Reset all logs when new experiment begins
		"""
		self.optimization_log 		= {}
		self.reference_signal 		= []
		self.log_reference_signal	= []
		
	def StartOptimization (self, event, DoGA) :
		"""
		Initiate GA. `self.odd_ga_button` was clinked
		`DoSingleIterationGA` is a function that performs the GA
		"""
		# Create pseudonyms of necessary devices 
		self.DevSpectrometer 	= self.parent.Spectrometer.dev
		self.DevPulseShaper		= self.parent.PulseShaper.dev
			
		# Save global settings and get the name of log file
		self.optimization_log_filename = SaveSettings(SettingsNotebook=self.parent, 
								title="Select file to save GA progress", filename="odd_experiment.hdf5")
		if self.optimization_log_filename is None : return
		
		####################### Initiate devices #############################
		
		# Initiate spectrometer
		settings = self.parent.Spectrometer.GetSettings()
		if self.DevSpectrometer.SetSettings(settings) == RETURN_FAIL : return
		
		# Initiate pulse shaper
		settings = self.parent.PulseShaper.GetSettings()
		if self.DevPulseShaper.Initialize(settings) == RETURN_FAIL : return
		
		# Get number of optimization variables 
		self.num_pixels = self.DevPulseShaper.GetParamNumber()
		if self.num_pixels == RETURN_FAIL :
			raise RuntimeError ("Optimization cannot be started since calibration file was not loaded")
		
		# Check whether the background signal array is present
		self.CheckBackground()
		
		ga_settings = self.GetSettings()
		
		#####################################################################
		
		# Open the file for optimization log to find checkpoint
		with  h5py.File (self.optimization_log_filename, 'a') as optimization_log_file :
			
			# The HDF5 group where each iteration is picked
			try : 
				optimization_iterations = optimization_log_file["optimization_iterations"]
				# loading last iteration, i.e., checkpoint
				self.current_iteration_number = max( int(key) for key in optimization_iterations.keys() )
				checkpoint = optimization_iterations[ str(self.current_iteration_number) ]
				self.current_iteration_number += 1
				# ask user whether to resume optimization
				if wx.MessageDialog(self, "Should optimization be resumed?", caption="resume optimization", 
					style=wx.YES_NO | wx.YES_DEFAULT ).ShowModal() == wx.ID_NO : raise ValueError
			except (KeyError, ValueError) :
				# Start fresh optimization, since there no iteration is pickled
				checkpoint = None
				self.current_iteration_number = 0
				# delete if it exists
				try : del optimization_log_file["optimization_iterations"]
				except KeyError : pass
				# create empty group
				optimization_log_file.create_group ("optimization_iterations")
		
		####################### Setting up GA #######################
		min_val = 0; max_val = 1
		
		creator.create("FitnessMax", base.Fitness, weights=(1.0,))
		creator.create("Individual", array.array, typecode='f', fitness=creator.FitnessMax)

		self.optimization_toolbox = base.Toolbox()
			
		# Attribute generator
		self.optimization_toolbox.register("attr_float", #random.random)
						lambda mu, sigma : min( max(random.gauss(mu,sigma), min_val), max_val), 
						#lambda mu, sigma : random.gauss(mu,sigma) % (max_val - min_val) + min_val , 
						ga_settings["gaussian_mutation_mu"], ga_settings["gaussian_mutation_sigma"] )
						#0.5*(max_val + min_val), ga_settings["gaussian_mutation_sigma"] )
		
		# Save the number of optimization variable
		self.num_params = self.Ind_length[ga_settings["pulse_shaping_option"]] \
							*(self.num_pixels / ga_settings["pixel_bundle_width"])
		
		# Structure remeasurers
		self.optimization_toolbox.register("individual", 
			tools.initRepeat, creator.Individual, self.optimization_toolbox.attr_float, self.num_params)
		
		# Define the function converting GA ind into the pulse shaper
		self.Ind2PulseShape = self.ind2pulse_shape[ ga_settings["pulse_shaping_option"] ]
		
		# Set post-processing spectrum function 
		self.SpectrumPostProcess = _SpectrumPostProcess.Select( ga_settings["spectrum_postprocess"] )
		
		self.optimization_toolbox.register("population", tools.initRepeat, list, self.optimization_toolbox.individual)
		self.optimization_toolbox.register("mate", tools.cxTwoPoint)
		
		self.optimization_toolbox.register("mutate", mutBoundedGaussian, mu=ga_settings["gaussian_mutation_mu"], 
			sigma=ga_settings["gaussian_mutation_sigma"], indpb=ga_settings["mutation_indpb"],  min_val=min_val, max_val=max_val)
		self.optimization_toolbox.register("select", tools.selBest)
		
		# Statistics
		self.optimization_stats = tools.Statistics(lambda ind: ind.fitness.values)
		self.optimization_stats.register("avg", np.mean)
		#self.optimization_stats.register("std", np.std)
		self.optimization_stats.register("min", np.min)
		self.optimization_stats.register("max", np.max)
		
		# Prepare logs
		self.ResetLogs()
		
		# Parameters for genetic algorithm
		self.cxpb 	= ga_settings["cxpb"]
		self.mutpb	= ga_settings["mutpb"]
		self.population_size = ga_settings["population_size"]
		
		if checkpoint :
			################# Load last iteration of optimization ##################
			if isinstance(checkpoint, HD5Group) :
				self.optimization_pop 		= pickle.loads( str(checkpoint["pickled_population"][...]) )
			else :
				self.optimization_pop 		= pickle.loads( checkpoint["pickled_population"] )
		else : 
			######################### Start new evolution ######################
			self.optimization_pop = self.optimization_toolbox.population(n=self.population_size)
			
		#####################################################################
		self.need_abort = False
		
		# Adjusting button's settings
		button = event.GetEventObject()
		button.SetLabel (button._stop_label)
		button.SetBackgroundColour('red')
		button.Bind( wx.EVT_BUTTON, self.StopOptimization)
		
		# Start doing a single iteration GA
		wx.CallAfter (DoGA) 
	
	def GetSampleSpectrum (self) :
		"""
		Measure sample fluorescence spectra 
		"""
		# Get spectra
		spectrum = self.DevSpectrometer.AcquiredData().astype(np.int)
			
		# Subtract the background
		spectrum  -= self.background_signal
		
		return self.SpectrumPostProcess(spectrum)
	
	def DoSingleIterGA_SHGMaximize (self, remeasure=True) :
		"""
		Perform a single iteration of GA to find phase of laser pulse 
		that coherently suppresses single sample fluorescence
		"""
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return
		
			
		if remeasure :		
			# Record spectrum
			self.optimization_pop = self.MeasureSpectra(self.optimization_pop)
			
			wx.Yield()
			# abort, if requested 
			if self.need_abort : return
			
			# Calculate the fitness the population
			for ind in self.optimization_pop :
				ind.fitness.values = (ind.spectra.sum(),)
			
			self.optimization_pop[:] = self.optimization_toolbox.select(
										self.optimization_pop, self.population_size
				)
		else :
			# Generate offspring
			offspring = algorithms.varAnd(self.optimization_pop, self.optimization_toolbox, self.cxpb, self.mutpb)
			
			# Measure spectra for offspring only
			offspring = self.MeasureSpectra(offspring)
			
			wx.Yield()
			# abort, if requested 
			if self.need_abort : return
			
			# Calculate the fitness  of offspring
			for ind in offspring :
				ind.fitness.values = (ind.spectra.sum(),)
				
			self.optimization_pop[:] = self.optimization_toolbox.select( 
				self.optimization_pop + offspring, int(1.2*self.population_size)
				)
				
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return
		
		# Saving statistics for plotting purposes
		for key, value in self.optimization_stats.compile(self.optimization_pop).items() :
			try : 
				self.optimization_log[key].append(value)
			except KeyError : 
				self.optimization_log[key] = [value]
		
		self.SaveOptimization()
		self.DisplayOptimization() 
	
		# Going to the next iteration
		wx.CallAfter (self.DoSingleIterGA_SHGMaximize,  remeasure=(not remeasure)) 
	
	def MeasureSpectra (self, individuals) :
		"""
		Measure spectra from population of pulse shapes in `individuals` 
		"""
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return individuals
		
		self.RecordReferenceSignal()
			
		for ind in individuals :
			# Set the pulse shaper
			self.DevPulseShaper.SetAmplPhase( *self.Ind2PulseShape(np.array(ind)) )
				
			wx.Yield()
			# abort, if requested 
			if self.need_abort : return individuals
				
			# Save spectra
			ind.spectra = self.GetSampleSpectrum()
		
		return individuals
	
	def SaveOptimization (self) : 
		"""
		Save current iteration into the log file
		"""
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return

		# Open the file for optimization log
		with  h5py.File (self.optimization_log_filename, 'a') as optimization_log_file :
			#################### Saving current optimization iteration ####################
			checkpoint = optimization_log_file["optimization_iterations"]
			checkpoint = checkpoint.create_group( str(self.current_iteration_number) )
	
			# Saving the population  
			checkpoint["pickled_population"] = np.string_( pickle.dumps(self.optimization_pop) )
				
			# Saving population in HDF5 format in the appropriate HDF5 group 
			individuals = checkpoint.create_group("individuals")
			for num, ind in enumerate(self.optimization_pop) :
				ind_grp = individuals.create_group( str(num) )
				ind_grp["pulse_shape"] 	= ind
				ind_grp["fitness"]		= ind.fitness.values
				
				# Saving multiple spectra data in a separate group
				spectra_grp = ind_grp.create_group("spectra")
				#for channel, spectrum in ind.spectra.iteritems() :
				#	spectra_grp[ str(channel) ] = spectrum
				spectra_grp[ "0" ] = ind.spectra
				
			#################### Saving log information ####################
			# Delete existing summary group
			try : del optimization_log_file["optimization_summary"] 
			except KeyError : pass
			
			# Create new group
			summary_group = optimization_log_file.create_group("optimization_summary")
			for key, val in self.optimization_log.iteritems() :
				summary_group[key] = val
			
			# Save the reference signal summary 
			try : del optimization_log_file["reference_signal_summary"] 
			except KeyError : pass
			
			summary_group = optimization_log_file.create_group("reference_signal_summary")
			#for key, val in self.log_reference_signal.iteritems() :
			#	summary_group[ "channel_%d" % key ] = val
			summary_group["chanel_0"] = self.log_reference_signal
			
			##########################################################################
				
		# Going to next iteration
		self.current_iteration_number += 1
	
	def DisplayOptimization (self) :
		"""
		Display the progress of optimization
		"""
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return
	
		def GetValueColourIter (d) :
			return izip_longest( d.itervalues(), ['r', 'g', 'b', 'k'],  fillvalue='y' )
	
		visvis.cla(); visvis.clf()
		visvis.subplot(211)
		
		# Plot optimization statistics
		for values, colour in GetValueColourIter(self.optimization_log) :
			try : visvis.plot ( values, lc=colour ) 
			except Exception : pass
		
		visvis.xlabel ('iteration')
		visvis.ylabel ('Objective function')
		visvis.legend( self.optimization_log.keys() )
		
		# Display reference signal
		visvis.subplot(212)
			
		# Plot reference signal
		try : visvis.plot ( self.log_reference_signal ) 
		except Exception : pass
		
		visvis.xlabel ('iteration')
		visvis.ylabel ("Signal from reference pulse")
			
	def StopOptimization (self, event) :
		"""
		Stop GA
		"""
		self.need_abort = True 
			
		# Adjusting button's settings
		button = event.GetEventObject()
		button.SetLabel (button._start_label)
		button.SetBackgroundColour('')
		button.Bind( wx.EVT_BUTTON, button._start_method)
	