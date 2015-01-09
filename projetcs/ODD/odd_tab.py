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

def RankingIndividuals (inds) :
	"""
	 This function is used to calculate the fitness function in the ODD experiment.
	 It ranks the tuple of individuals `inds` 
	 It is assumed that ind.fitness.values corresponds to to the number of the individual `ind` in
	 the population.
	"""
	# Iterator generating fluorescence matrix for the current pulse combination 
	fluorescence_matrix = izip( *chain(*[ind.spectra.itervalues() for ind in inds]) )
	N = len(inds)
	fluorescence_matrix = ( np.array(F).reshape((N,N)) for F in fluorescence_matrix )
		
	# calculate the rank (i.e., fitness) of the pulse tuple, 
	# and the pixel number in spectra 
	rank, pix_num = max( ( abs(np.linalg.det(F)), pix_num ) 
							for pix_num, F in enumerate(fluorescence_matrix) )
									
	# Find tuple of incidences numbering the utilized individuals
	ind_tuple = tuple( ind.pix_num for ind in inds )
									
	return rank, ind_tuple, pix_num

########################################################################

class ODD_Tab (HardwareGUIControl) :
	"""
	GUI for ODD algorithm 
	"""
	def __init__ (self, parent) :
	
		HardwareGUIControl.__init__(self, parent)
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		# List of positions of channels
		sizer.Add (wx.StaticText(self, label="Channel number with pure samples for learning"), flag=wx.LEFT, border=5)
		self.chanel_odd_experiment_ctrl = wx.TextCtrl (self, value="", style=wx.TE_MULTILINE|wx.EXPAND)
		self.chanel_odd_experiment_ctrl.__label__ = "channels"
		sizer.Add (self.chanel_odd_experiment_ctrl, flag=wx.EXPAND, border=5)
		
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
		
		# Fitness function options
		sizer.Add ( wx.StaticText(self, label="Fitness function") )
		self.fitness_options = {
			"determinant"		:	lambda F : np.abs(np.linalg.det(F)),
			"min eigenvalue" 	:	lambda F : np.abs(np.linalg.eigvals(F)).min(),
			"max eigenvalue"	:	lambda F : np.abs(np.linalg.eigvals(F)).max()
		}
		choices = self.fitness_options.keys()
		fitness_options_ctrl = wx.ComboBox (self, choices=choices, value=choices[1], style=wx.CB_READONLY )
		fitness_options_ctrl.__label__ = "fitness_options"
		sizer.Add (fitness_options_ctrl,  flag=wx.EXPAND, border=5)
		
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
		
		# ODD GA button
		self.odd_ga_button = wx.Button (self)
		self.odd_ga_button._start_label = "Start ODD optimization"
		self.odd_ga_button._start_method = partial(self.StartOptimization, DoGA=self.DoSingleIterGA_ODD)
		self.odd_ga_button._stop_label = "STOP ODD optimization"
		self.odd_ga_button.SetLabel (self.odd_ga_button._start_label)
		self.odd_ga_button.Bind( wx.EVT_BUTTON, self.odd_ga_button._start_method)
		sizer.Add(self.odd_ga_button, flag=wx.EXPAND, border=5)
		
		# Coherently suppress fluorescence button
		self.sup_fluor_button = wx.Button (self)
		self.sup_fluor_button._start_label = "Start GA to suppress fluorescence" 
		self.sup_fluor_button._start_method = partial(self.StartOptimization, DoGA=self.DoSingleIterGA_CoherentSuppression)
		self.sup_fluor_button._stop_label = "STOP GA to suppress fluorescence"
		self.sup_fluor_button.SetLabel (self.sup_fluor_button._start_label)
		self.sup_fluor_button.Bind( wx.EVT_BUTTON, self.sup_fluor_button._start_method )
		sizer.Add(self.sup_fluor_button, flag=wx.EXPAND, border=5)
		
		###################### Measuring concentrations ######################
		sb_sizer = wx.StaticBoxSizer( wx.StaticBox(self, label="Concentration measurements"),  wx.VERTICAL )
		
		# List of positions of channels_mixtures
		sb_sizer.Add (wx.StaticText(self, label="Channel number to measure concentrations"), flag=wx.LEFT, border=5)
		self.channel_mixtues_ctrl = wx.TextCtrl (self, value="", style=wx.TE_MULTILINE|wx.EXPAND)
		self.channel_mixtues_ctrl.__label__ = "channels_mixtures"
		sb_sizer.Add (self.channel_mixtues_ctrl, flag=wx.EXPAND, border=5)
		
		# ODD GA file
		sb_sizer.Add (wx.StaticText(self, label="ODD GA file"), flag=wx.LEFT, border=5)
		odd_ga_file_ctr = wx.FilePickerCtrl(self, message="Chose file containing optimization data...")
		odd_ga_file_ctr.__label__ = "odd_ga_file_name"
		sb_sizer.Add (odd_ga_file_ctr, flag=wx.EXPAND, border=5)
		
		# Number of pulses to be used
		sb_sizer.Add (wx.StaticText(self, label="Number of pulses to interrogate"), flag=wx.LEFT, border=5)
		num_pulses_ctrl  = wx.SpinCtrl (self, value="10", min=1, max=100000)
		num_pulses_ctrl.__label__ = "num_pulses"
		sb_sizer.Add (num_pulses_ctrl , flag=wx.EXPAND, border=5)
		
		# Measure concentration button (OLD way)
		measure_concentr_button = wx.Button (self)  
		measure_concentr_button._start_label 	= "OLD Measure concentration"
		measure_concentr_button._start_method 	= self.Old_MeasureConcentration
		measure_concentr_button._stop_label 	= "STOP measuring concentration"
		measure_concentr_button._stop_method 	= self.Stop_MeasureConcentration
		measure_concentr_button.SetLabel (measure_concentr_button._start_label)
		measure_concentr_button.Bind (wx.EVT_BUTTON, measure_concentr_button._start_method)
		sb_sizer.Add(measure_concentr_button, flag=wx.EXPAND, border=5) 
		
		# Measure concentration button (NEW way)
		measure_concentr_button = wx.Button (self)  
		measure_concentr_button._start_label 	= "Measure concentration"
		measure_concentr_button._start_method 	= self.MeasureConcentration
		measure_concentr_button._stop_label 	= "STOP measuring concentration"
		measure_concentr_button._stop_method 	= self.Stop_MeasureConcentration
		measure_concentr_button.SetLabel (measure_concentr_button._start_label)
		measure_concentr_button.Bind (wx.EVT_BUTTON, measure_concentr_button._start_method)
		sb_sizer.Add(measure_concentr_button, flag=wx.EXPAND, border=5) 
		
		sizer.Add (sb_sizer, flag=wx.EXPAND, border=5)
		######################################################################
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		self.CreateSettingsDict()

	def RecordBackground (self, event=None) :
		"""
		Record background spectrum
		"""
		# Create pseudonyms 
		self.DevSpectrometer 	= self.parent.Spectrometer.dev
		self.DevSampleSwitcher	= self.parent.SampleSwitcher.dev
		
		# Initiate spectrometer
		settings = self.parent.Spectrometer.GetSettings()
		if self.DevSpectrometer.SetSettings(settings) == RETURN_FAIL : return
		
		# Initiate sample switcher
		settings = self.parent.SampleSwitcher.GetSettings()
		if self.DevSampleSwitcher.Initialize(settings) == RETURN_FAIL : return
		
		# Saving the name of channels 
		ga_settings = self.GetSettings()
		self.channels = sorted(eval( "(%s,)" % ga_settings["channels"] ))
		if self.DevSampleSwitcher.GetChannelNum()-1 < max(self.channels) :
			raise ValueError ("Error: Some channels specified are not accessible by sample switcher.")
			
		# Record background for each available channels 
		self.background_signal = {}
		for channel in self.channels :
			self.DevSampleSwitcher.MoveToChannel(channel)
			self.background_signal[ channel ] = self.DevSpectrometer.AcquiredData().astype(np.int)
	
	def RecordReferenceSignal (self, channel) :
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
		spectrum = self.GetSampleSpectrum(channel)
		self.reference_signal[ channel ] = spectrum
			
		try :
			self.log_reference_signal[ channel ].append( spectrum.sum() )
		except KeyError :
			self.log_reference_signal[ channel ] = [ spectrum.sum() ]
			
	def CheckBackground (self) :
		"""
		Check whether the background signal is recorded and ready to be used. 
		"""
		try : 
			# The background signal must be consistent with self.channels
			for channel in self.channels :
				if channel not in self.background_signal : raise AttributeError
		except AttributeError :
			def SetBackgroundZero () :
				self.background_signal = dict( (channel, 0) for channel in self.channels ) 
				
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
		self.optimization_log 		= { }
		self.reference_signal 		= { }
		self.log_reference_signal	= { }
		self.emission_spectra 		= { }
		self.N_emission_spectra		= { }
		
	def StartOptimization (self, event, DoGA) :
		"""
		Initiate GA. `self.odd_ga_button` was clinked
		`DoSingleIterationGA` is a function that performs the GA
		"""
		# Create pseudonyms of necessary devices 
		self.DevSpectrometer 	= self.parent.Spectrometer.dev
		self.DevPulseShaper		= self.parent.PulseShaper.dev
		self.DevSampleSwitcher	= self.parent.SampleSwitcher.dev
			
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
		
		# Initiate sample switcher
		settings = self.parent.SampleSwitcher.GetSettings()
		if self.DevSampleSwitcher.Initialize(settings) == RETURN_FAIL : return
		
		# Saving the name of channels 
		ga_settings = self.GetSettings()
		self.channels = sorted(eval( "(%s,)" % ga_settings["channels"] ))
		if self.DevSampleSwitcher.GetChannelNum()-1 < max(self.channels) :
			raise ValueError ("Error: Some channels specified are not accessible by sample switcher.")
		
		# Check whether the background signal array is present
		self.CheckBackground()
		
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
						lambda mu, sigma : random.gauss(mu,sigma) % (max_val - min_val) + min_val , 
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
		
		# Define fitness function
		self.FitnessFunction = self.fitness_options[ ga_settings["fitness_options"] ] 
		
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
	
	def GetSampleSpectrum (self, channel) :
		"""
		Measure sample fluorescence spectra 
		"""
		# Get spectra
		spectrum = self.DevSpectrometer.AcquiredData().astype(np.int)
		
		# The following block is to obtain a super long time averaged 
		# emission spectra of molecules in each channel
		
		# The mean is calculated iteratively 
		# see, e.g., http://www.heikohoffmann.de/htmlthesis/node134.html
		try :
			self.N_emission_spectra[channel] += 1
			self.emission_spectra[channel] += ( spectrum - self.emission_spectra[channel] )/ self.N_emission_spectra[channel]
		except KeyError :
			self.emission_spectra[channel] = spectrum.astype(np.float)
			self.N_emission_spectra[channel] = 1
			
		# Subtract the background
		spectrum  -= self.background_signal[channel]
		
		return self.SpectrumPostProcess(spectrum)
	
	def DoSingleIterGA_CoherentSuppression (self, remeasure=True) :
		"""
		Perform a single iteration of GA to find phase of laser pulse 
		that coherently suppresses single sample fluorescence
		"""
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return
		
		# Consistency check
		assert len(self.channels) == 1, "Only one channel must be specified" 
		channel = self.channels[0]
		
		# Move to a selected channel 
		self.DevSampleSwitcher.MoveToChannel(channel)
			
		if remeasure :		
			# Record spectrum
			self.optimization_pop = self.MeasureSpectra(self.optimization_pop)
			
			wx.Yield()
			# abort, if requested 
			if self.need_abort : return
			
			# Calculate the fitness the population
			for ind in self.optimization_pop :
				ind.fitness.values = (-ind.spectra[channel].sum(),)
			
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
				ind.fitness.values = (-ind.spectra[channel].sum(),)
				
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
		wx.CallAfter (self.DoSingleIterGA_CoherentSuppression,  remeasure=(not remeasure)) 
	
	def DoSingleIterGA_ODD (self, remeasure=True) :
		"""
		Perform a single iteration of GA to find pulses for discrimination
		"""
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return
		
		if remeasure :
			# Measure the population
			self.optimization_pop = self.MeasureSpectra(self.optimization_pop)
			self.optimization_pop = self.CalculateFitness(self.optimization_pop)
			
			self.optimization_pop[:] = self.optimization_toolbox.select(
										self.optimization_pop, self.population_size
				)
		else :
			# Generate offspring
			offspring = algorithms.varAnd(self.optimization_pop, self.optimization_toolbox, self.cxpb, self.mutpb)
			
			# Measure spectra for offspring only
			offspring = self.MeasureSpectra(offspring)
			# Selection
			self.optimization_pop[:] = self.optimization_toolbox.select(
				self.CalculateFitness(self.optimization_pop + offspring), int(1.2*self.population_size)
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
		wx.CallAfter (self.DoSingleIterGA_ODD, remeasure=(not remeasure)) 
			
	def MeasureSpectra (self, individuals) :
		"""
		Measure fluorescence spectra emitted by each sample channel 
		from population of pulse shapes in `individuals` 
		"""
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return individuals
		
		# Cleaning up the previously collected spectra
		for ind in individuals :
			ind.spectra = dict()
			
		for channel in self.channels :
			# Move to a channel containing pure sample 
			self.DevSampleSwitcher.MoveToChannel(channel)
		
			self.RecordReferenceSignal(channel)
			
			for ind in individuals :
				# Set the pulse shaper
				self.DevPulseShaper.SetAmplPhase( *self.Ind2PulseShape(np.array(ind)) )
				
				wx.Yield()
				# abort, if requested 
				if self.need_abort : return individuals
				
				# Save spectra
				ind.spectra[channel] = self.GetSampleSpectrum(channel)
		
		return individuals
	
	def CalculateFitness (self, individuals) :
		"""
		Calculate fitness function for `individuals`
		"""
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return individuals
		
		# List storing tuple ranking
		pulse_tuple_ranks = []
		
		N = len(self.channels)
		
		# Loop over all possible combination of laser pulses
		for ind_tuple in combinations( range(len(individuals)), N ) :
			
			# Iterator generating fluorescence matrix for the current pulse combination 
			fluorescence_matrix = izip( *chain(*[individuals[indx].spectra.itervalues() for indx in ind_tuple]) )
			fluorescence_matrix = ( np.array(F).reshape((N,N)) for F in fluorescence_matrix )
			
			# calculate the rank (i.e., fitness) of the pulse tuple 
			rank = max( self.FitnessFunction(F) for F in fluorescence_matrix )
		
			pulse_tuple_ranks.append( ( rank, ind_tuple) )
			
			wx.Yield()
			# abort, if requested 
			if self.need_abort : return individuals
		
		# Sort tuples by rank
		pulse_tuple_ranks.sort( reverse=True )
		
		# Indexes of individuals for which fitness value has been assigned
		marked_ind = set()
		
		# Loop over ranked tuples
		for rank, ind_tuple in pulse_tuple_ranks :
		
			ind_tuple = set(ind_tuple)
			for indx in ind_tuple.difference( marked_ind ) :
				individuals[indx].fitness.values = (rank,)
			
			marked_ind.update( ind_tuple )
			
			# Stop looping if all fitness functions have been calculated
			if len(marked_ind) == len(individuals): break
		
		return individuals
		
	def TOBE_FIXED_CalculateFitness (self, individuals) :
		"""
		Calculate fitness function for `individuals`
		"""
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return individuals
		
		# Number individuals in the population
		# this is required by the function RankingIndividuals
		for num, ind in enumerate(individuals) :
			ind.pix_num = num
		
		# Starting multiprocessor mapper
		multiprocessing.freeze_support()
		pool = multiprocessing.Pool( processes=(multiprocessing.cpu_count()-1) )
		
		# Iterator for combinatorial analysis
		combinatorial_iter = combinations( individuals, len(self.channels) )
		
		############ Compute pulse ranking using multiprocessing pooling ###########
		pulse_tuple_ranked = []
		while True :
			# Regarding `islice` see 
			# http://stackoverflow.com/questions/5318936/python-multiprocessing-pool-lazy-iteration
			calculated = pool.map(RankingIndividuals, islice(combinatorial_iter, 100))
			
			# Check whether combinatorial iterations are finished
			if len(calculated) == 0 : break
			
			wx.Yield()
			# abort, if requested 
			if self.need_abort : return individuals
			
			########### Process calculated data ###########
			# Sort all tuples by rank
			pulse_tuple_ranked.extend( calculated )
			pulse_tuple_ranked.sort( reverse=True )
			
			# Find the minimal value of the index in `pulse_tuple_ranked` list 
			# up to which all the individuals are accounted for
			cut_off = 0
			
			marked_ind = set()
			for _, ind_tuple, _ in pulse_tuple_ranked :
				marked_ind.update( ind_tuple )
				# Check whether all individuals are accounted for in `marked_ind`
				if len(marked_ind) == len(individuals) : break
				else : cut_off += 1
			
			# Keep only minimal number of tuples necessary
			del pulse_tuple_ranked[ (cut_off+1): ]
			
		####################### Converting ranks into fitness #######################
		
		# Indexes of individuals for which fitness value has been assigned
		marked_ind = set()
		
		# Loop over ranked tuples
		for rank, ind_tuple, pix_num in pulse_tuple_ranks :
		
			ind_tuple = set(ind_tuple)
			for indx in ind_tuple.difference( marked_ind ) :
				# Assign fitness
				individuals[indx].fitness.values = (rank,)
				# Save the optimal pixel number 
				individuals[indx].pix_num = pix_num		
				
			marked_ind.update( ind_tuple )
		
		#####################################################################
		return individuals	
		
	def __CalculateFitness__ (self, individuals) :
		"""
		Calculate fitness function for `individuals` using multiprocessor map
		"""
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return individuals
		
		####################### Save `individuals` into a file #######################
		# so that the multiprocess mapper could access them
		individuals_filename = '_tmp_idividuals_.hdf5'

		with h5py.File (individuals_filename, 'w') as individuals_file :
			for indx, ind in enumerate(individuals) :
				# Create individual for each individual's spectra dict
				ind_grp = individuals_file.create_group( str(indx) ) 
				
				# Save spectra 
				for channel, spectrum in ind.spectra.iteritems() :
					ind_grp[ str(channel) ] = spectrum
		
		####################### Mapper #######################
		def Ranking (ind_tuple) :
			"""
			This function ranks of the tuple `ind_tuple`.
			Returns (tuple's rank, tuple, and spectral pixel used for rank calculations)
			"""
			# Load data about the tuple
			with h5py.File (individuals_filename, 'r') as individuals_file :
				# Load all spectra for the specified tuple
				spectra = [ spectrum[...] for spectrum in 
					chain( *[ individuals_file[str(indx)].itervalues() for indx in ind_tuple ] )
				]
			
			N = int(np.sqrt( len(spectra) ))
			
			# Iterator generating fluorescence matrix for the given tuple
			fluorescence_matrix = ( np.array(F).reshape((N,N)) for F in izip(*spectra) )
			
			# calculate the rank (i.e., fitness) of the pulse tuple, 
			# and the pixel number in spectra 
			rank, pix_num = max( ( abs(np.linalg.det(F)), pix_num ) 
									for pix_num, F in enumerate(fluorescence_matrix) )
			
			return rank, ind_tuple, pix_num
		
		# Starting multiprocessor mapper
		pool = multiprocessing.Pool( processes=(multiprocessing.cpu_count()-1) )  
		# Iterator for combinatorial analysis
		combinatorial_iter = combinations( range(len(individuals)), len(self.channels) )
		pulse_tuple_ranked = pool.map(Ranking, combinatorial_iter)
		#pulse_tuple_ranked = pool.map(Ranking, islice(combinatorial_iter, 1000))
		
		####################### Converting ranks into fitness #######################
		# Sort tuples by rank
		pulse_tuple_ranked.sort()
		
		# Indexes of individuals for which fitness value has been assigned
		marked_ind = set()
		
		# Loop over ranked tuples until all fitness functions are assigned
		while len(marked_ind) < len(individuals) :
			rank, ind_tuple, pix_num = pulse_tuple_ranked.pop()
		
			ind_tuple = set(ind_tuple)
			for indx in ind_tuple.difference( marked_ind ) :
				# Assign fitness
				individuals[indx].fitness.values = (rank,)
				# Save the optimal pixel number 
				# Save the optimal pixel number 
				individuals[indx].pix_num = pix_num
			marked_ind.update( ind_tuple )
			
		
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
				for channel, spectrum in ind.spectra.iteritems() :
					spectra_grp[ str(channel) ] = spectrum
			
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
			for key, val in self.log_reference_signal.iteritems() :
				summary_group[ "channel_%d" % key ] = val
				
			# Save super long time averages of emission spectra  
			try : del optimization_log_file["emission_spectra"] 
			except KeyError : pass
			
			emission_spectra_grp = optimization_log_file.create_group("emission_spectra")
			for key, val in self.emission_spectra.iteritems() :
				emission_spectra_grp[ "channel_%d" % key ] = val
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
		for values, colour in GetValueColourIter(self.log_reference_signal) :
			try : visvis.plot ( values, lc=colour ) 
			except Exception : pass
		
		visvis.xlabel ('iteration')
		visvis.ylabel ("Signal from reference pulse")
		visvis.legend( ["channel %d" % x for x in self.log_reference_signal.keys()] )

			
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
		
	def Old_MeasureConcentration (self, event) :
		"""
		Perform concentration measurement 
		"""
		# Create pseudonyms of necessary devices 
		self.DevSpectrometer 	= self.parent.Spectrometer.dev
		self.DevPulseShaper		= self.parent.PulseShaper.dev
		self.DevSampleSwitcher	= self.parent.SampleSwitcher.dev
		
		# Save global settings and get the name of log file
		measurement_filename = SaveSettings(SettingsNotebook=self.parent, 
								title="Select file to save concentration measurements", filename="concentrations.hdf5")
		if measurement_filename is None : return
		
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
		
		# Initiate sample switcher
		settings = self.parent.SampleSwitcher.GetSettings()
		if self.DevSampleSwitcher.Initialize(settings) == RETURN_FAIL : return
		
		settings = self.GetSettings()
		
		# Load channels with mixtures 
		channels_mixtures 	= sorted( eval( "(%s,)" % settings["channels_mixtures"] ) )
		channels_pure		= sorted( eval( "(%s,)" % settings["channels"] ) )
		
		# Saving the name of all channels containing pure and mixed samples
		self.channels = sorted(set( channels_pure + channels_mixtures  ))
		if self.DevSampleSwitcher.GetChannelNum()-1 < max(self.channels) :
			raise ValueError ("Error: Some channels specified are not accessible by sample switcher.")
		
		# Check whether the background signal array is present
		self.CheckBackground()
		
		# Adjusting button's settings
		self.need_abort = False
		
		button = event.GetEventObject()
		button.SetLabel (button._stop_label)
		button.Bind( wx.EVT_BUTTON, button._stop_method)
		
		############### Load optimized pulse shape ##########################
		with h5py.File (settings["odd_ga_file_name"], 'r') as odd_ga_file :
		
			ga_settings_grp = odd_ga_file["settings/ODD_GA/"]
			
			# Define the function converting GA ind into the pulse shaper
			self.Ind2PulseShape = self.ind2pulse_shape[ str(ga_settings_grp["pulse_shaping_option"][...]) ] 
		
			# Save the number of optimization variable
			self.num_params = self.Ind_length[ str(ga_settings_grp["pulse_shaping_option"][...]) ] \
							*(self.num_pixels / ga_settings_grp["pixel_bundle_width"][...])
		
			# Select pulse shapes to be used for measuring concentration
			optimization_iterations = odd_ga_file["optimization_iterations"]
			choices = sorted( optimization_iterations.keys(), key=int, reverse=True )
			dlg = wx.MultiChoiceDialog( self, "Pick GA iterations to load pulse shapes from",
                                   "Photonic reagents", choices)							   
			if dlg.ShowModal() != wx.ID_OK :
				# user cancel
				return
			
			# Load choices
			choices = [ optimization_iterations["%s/individuals" % choices[x]].itervalues() for x in dlg.GetSelections() ]
			
			# Loading laser pulses from selected GA iteration 
			loaded_pulse_shapes = [ ( x["fitness"][...], x["pulse_shape"][...] ) for x in chain(*choices) ]
			
		# Sort loaded pulses by fitness function
		loaded_pulse_shapes.sort( key=itemgetter(0) )
				
		# Find unique pulse shapes
		pulse_shapes = set()
		while len(pulse_shapes) < settings["num_pulses"] and len(loaded_pulse_shapes) :
			# ignore the fitness
			_, ps = loaded_pulse_shapes.pop()
			pulse_shapes.update( ( tuple(ps), ) )
		
		# Set post-processing spectrum function 
		self.SpectrumPostProcess = _SpectrumPostProcess.Select( settings["spectrum_postprocess"] )
		
		# This lines allow us to reuse the rest of the code for ODD GA 
		creator.create("Individual", np.ndarray, spectra=dict)
		pulse_shapes = map( creator.Individual, pulse_shapes ) 

		################ Concentration determination ################
		
		# Measure spectra from all channels 
		self.ResetLogs()
		self.MeasureSpectra( pulse_shapes )
		
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return 
		
		# Form fluorescent matrix of pure samples 
		fluoresence_matrix = np.vstack( [
			np.hstack( [pulse.spectra[channel].flat for pulse in pulse_shapes] ) for channel in channels_pure  
		] ).T
		
		# This is for debugging purposes
		
		self.MeasureSpectra( pulse_shapes )
		
		# Calculate concentrations and save the results 
		with h5py.File (measurement_filename, 'a') as measurement_file :
			
			# Create grope where results of concentration measurements will be saved 
			try : del measurement_file["concentration_measurements"]
			except KeyError : pass
			concentration_measurements_grp = measurement_file.create_group("concentration_measurements")
			
			# Perform calculations
			for channel in channels_mixtures :
				# Create group for individual channel containing mixtures 
				channel_grp = concentration_measurements_grp.create_group("channel_%d" % channel)
			
				# Grope fluorescence spectral data for the studied mixture 
				mixture_fluorescence = np.hstack( [ pulse.spectra[channel].flat for pulse in pulse_shapes ] )
			
				# least square solution using all measured fluorescence data
				channel_grp["least square concentrations"] = \
					np.linalg.lstsq( fluoresence_matrix, mixture_fluorescence )[0]
				
				# non-negative least square solution
				channel_grp["non-negative least square concentrations"] = \
					nnls( fluoresence_matrix, mixture_fluorescence )[0]
					
				# print the results
				print "\nResults of concentration determination for channel %d" % channel
				for method, results in channel_grp.iteritems() :
					print "\t %s: \n\t\t %s" % (method, "; ".join("%.2e" % r for r in results[...]) )
		
		
		# Adjusting button's settings
		button = event.GetEventObject()
		button.SetLabel (button._start_label)
		button.Bind( wx.EVT_BUTTON, button._start_method)
	
	def MeasureConcentration (self, event) :
		"""
		Perform concentration measurement 
		"""
		# Create pseudonyms of necessary devices 
		self.DevSpectrometer 	= self.parent.Spectrometer.dev
		self.DevPulseShaper		= self.parent.PulseShaper.dev
		self.DevSampleSwitcher	= self.parent.SampleSwitcher.dev
		
		# Save global settings and get the name of log file
		measurement_filename = SaveSettings(SettingsNotebook=self.parent, 
								title="Select file to save concentration measurements", filename="concentrations.hdf5")
		if measurement_filename is None : return
		
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
		
		# Initiate sample switcher
		settings = self.parent.SampleSwitcher.GetSettings()
		if self.DevSampleSwitcher.Initialize(settings) == RETURN_FAIL : return
		
		settings = self.GetSettings()
		
		# Load channels with mixtures 
		channels_mixtures 	= sorted( eval( "(%s,)" % settings["channels_mixtures"] ) )
		channels_pure		= sorted( eval( "(%s,)" % settings["channels"] ) )
		
		# Saving the name of all channels containing pure and mixed samples
		self.channels = sorted(set( channels_pure + channels_mixtures  ))
		if self.DevSampleSwitcher.GetChannelNum()-1 < max(self.channels) :
			raise ValueError ("Error: Some channels specified are not accessible by sample switcher.")
		
		# Check whether the background signal array is present
		self.CheckBackground()
		
		# Adjusting button's settings
		self.need_abort = False
		
		button = event.GetEventObject()
		button.SetLabel (button._stop_label)
		button.Bind( wx.EVT_BUTTON, button._stop_method)
		
		########################## Load pulse shapes ##########################
		dlg = LoadPulseShapesDialog(parent=self, title="Select pulse shapes for measuring concentrations")
		dlg.ShowModal()
		loaded_pulse_shapes = dlg.GetLoadedData()
		dlg.Destroy()
		if len(loaded_pulse_shapes) == 0 : 
			# No pulses loaded, then exit
			return
		
		########################## Analyze pulse shapes ##########################
		
		# Find out what kind of pulse shaping is more popular
		options = Counter( ( str(key) for key, _ in loaded_pulse_shapes ) )
		
		if len(options) > 1 :
			print "Warning: Different pulse shaping options have are in the loaded pulses. We select the most popular. The other will be ignored."
		shaping = options.most_common(1)[0]
		
		# Keep pulse shapes only with the right shaping 
		pulse_shapes = [ pulse for shaping_option, pulse in loaded_pulse_shapes if str(shaping_option) == shaping ]
		
		# Define the function converting GA ind into the pulse shaper
		self.Ind2PulseShape = self.ind2pulse_shape[ shaping ] 
		
		############################################################
		
		# Set post-processing spectrum function 
		self.SpectrumPostProcess = _SpectrumPostProcess.Select( settings["spectrum_postprocess"] )
		
		# This lines allow us to reuse the rest of the code for ODD GA 
		creator.create("Individual", np.ndarray, spectra=dict)
		pulse_shapes = map( creator.Individual, pulse_shapes ) 

		################ Concentration determination ################
		
		# Measure spectra from all channels 
		self.ResetLogs()
		self.MeasureSpectra( pulse_shapes )
		
		wx.Yield()
		# abort, if requested 
		if self.need_abort : return 
		
		# Form fluorescent matrix of pure samples 
		fluoresence_matrix = np.vstack( [
			np.hstack( [pulse.spectra[channel].flat for pulse in pulse_shapes] ) for channel in channels_pure  
		] ).T
		
		# This is for debugging purposes
		
		self.MeasureSpectra( pulse_shapes )
		
		# Calculate concentrations and save the results 
		with h5py.File (measurement_filename, 'a') as measurement_file :
			
			# Create grope where results of concentration measurements will be saved 
			try : del measurement_file["concentration_measurements"]
			except KeyError : pass
			concentration_measurements_grp = measurement_file.create_group("concentration_measurements")
			
			# Perform calculations
			for channel in channels_mixtures :
				# Create group for individual channel containing mixtures 
				channel_grp = concentration_measurements_grp.create_group("channel_%d" % channel)
			
				# Grope fluorescence spectral data for the studied mixture 
				mixture_fluorescence = np.hstack( [ pulse.spectra[channel].flat for pulse in pulse_shapes ] )
			
				# least square solution using all measured fluorescence data
				channel_grp["least square concentrations"] = \
					np.linalg.lstsq( fluoresence_matrix, mixture_fluorescence )[0]
				
				# non-negative least square solution
				channel_grp["non-negative least square concentrations"] = \
					nnls( fluoresence_matrix, mixture_fluorescence )[0]
					
				# print the results
				print "\nResults of concentration determination for channel %d" % channel
				for method, results in channel_grp.iteritems() :
					print "\t %s: \n\t\t %s" % (method, "; ".join("%.2e" % r for r in results[...]) )
		
		
		# Adjusting button's settings
		button = event.GetEventObject()
		button.SetLabel (button._start_label)
		button.Bind( wx.EVT_BUTTON, button._start_method)
		
	def Stop_MeasureConcentration (self, event) :
		"""
		Stop measuring concentration 
		"""
		self.need_abort = True
		
		# Adjusting button's settings
		button = event.GetEventObject()
		button.SetLabel (button._start_label)
		button.Bind( wx.EVT_BUTTON, button._start_method)
		