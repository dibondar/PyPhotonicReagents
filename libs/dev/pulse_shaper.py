########################################################################
#
#	This file contains classes for controlling CRi SLM pulse shaper
#
########################################################################

from libs.gui.hardware_control import HardwareGUIControl
from libs.gui.load_data_set_from_file import LoadPulseShapesDialog
from libs.dev.basic_device import BasicDevice 
from libs.dev.consts import * 

import numpy as np
from scipy.interpolate import PchipInterpolator
import wx, serial, multiprocessing, h5py, time
from itertools import izip
from collections import Counter

########################################################################
#
#	Manager of the shaper
#
########################################################################

# Int type used to denote the shaper voltages
# Note that masks be np.zeros(640, dtype=ShaperInt)
ShaperInt = np.dtype('<u2')

# Maximal voltage 
PULSESHAPER_MAX_VAL = 4095 



class ManagerShaper :
	"""
	Class that manges the CRi SLM shaper
	"""
	def __init__ (self) :
		# Create the lock for device 
		self.lock = multiprocessing.Lock()
		# Create a pipe for communication
		self.parent_connection, self.child_connection = multiprocessing.Pipe()
	
	def __del__ (self) :
		self.parent_connection.close()
		self.child_connection.close()	
	
	def start(self) :
		"""
		Start the process controlling the shaper
		"""
		p = PulseShaper(self.child_connection)
		p.start()
		return p
	
	def run(self, command, arguments=None) :
		"""
		Send the command to the shaper through the pipe
		"""
		self.lock.acquire()
		self.parent_connection.send( (command, arguments) )
		result = self.parent_connection.recv()
		self.lock.release()
		return result
		
	def exit(self) :
		"""
		Close the process
		"""
		self.StopDevice()
		return self.run("Exit")
	
	def Initialize(self, settings) :
		"""
		Initialize shaper 
		"""
		return self.run("Initialize", settings)
	
	def StopDevice(self) :
		return self.run("StopDevice")
	
	def SetMasks (self, master_mask, slave_mask) :
		return self.run ("SetMasks", (master_mask, slave_mask) )
	
	def SetAmplPhase (self, amplitude, phase) :
		"""
		Set amplitudes and phases based on the calibration data.
		Calibration file must be loaded.
		"""
		return self.run ("SetAmplPhase", (amplitude, phase) )
		
	def GetUnwrappedAmplPhase (self, amplitude, phase) :
		"""
		Get actual values of phases (in rads) and amplitude (0, 1)
		"""
		return self.run ("GetUnwrappedAmplPhase", (amplitude, phase) )
	
	def GetParamNumber (self) :
		"""
		Return number of logical pixels of the shaper as a result of calibration.
		Calibration file must be loaded.
		"""
		return self.run ("GetParamNumber")
	
	def SetUniformMasks (self, vol_master_mask, vol_slave_mask) :
		"""
		Apply uniform voltage onto masks of the pulse shaper
		"""
		return self.run ("SetUniformMasks", (vol_master_mask, vol_slave_mask))
	
	def GetPixelNumber (self) :
		"""
		Returns number of pixels of the pulse shaper
		"""
		return self.run ("GetPixelNumber")
		
	def GetZeroMask (self) :
		"""
		Returns mask with zero values
		"""
		return np.zeros( self.GetPixelNumber(), dtype=ShaperInt )

########################################################################
#
#	Classes regulating usage of calibration data
#
########################################################################

def Order (P,V) :
	""" Order both `P` and `V` by `P` """
	indx = np.argsort(P)
	return P[indx], V[indx]
	
def GetTotalPhaseBounds (min_phase_master_mask, max_phase_master_mask, 
						min_phase_slave_mask, max_phase_slave_mask) :
	"""
	Get min and max total phase that can be generated given the pulse shaper calibration data
	"""
	phase_min = 2*np.vstack( (	min_phase_master_mask, 			min_phase_slave_mask, 
								min_phase_master_mask - np.pi/4, min_phase_slave_mask + np.pi/4. )
							).max(axis=0)
		
	phase_max = 2*np.vstack( (	max_phase_master_mask, 			max_phase_slave_mask, 
								max_phase_master_mask - np.pi/4, max_phase_slave_mask + np.pi/4. )
							).min(axis=0)
							
	return phase_min, phase_max

def FindBrokenPixels (min_phase_master_mask, max_phase_master_mask, 
					 min_phase_slave_mask, max_phase_slave_mask) :
	"""
	Return positions of  "broken" pixels on the criterion that 
	the phases of slave and master masks should be approximately equal (i.e, with 10%)
	"""
	RelVariation = lambda A, B : np.abs(A - B) / np.maximum(A, B)
	return np.nonzero( ( RelVariation(min_phase_master_mask, min_phase_slave_mask) > 0.1 )| 
					   ( RelVariation(max_phase_master_mask, max_phase_slave_mask) > 0.1 ) )[0]
	
################### The following class is under development ###################									 
class AmplPhase2ShaperMasks_Pixelwise :
	"""
	Class for converting amplitudes and phases into the pulse shaper masks 
	where the calibration data is used for each pixels separately
	"""
	
	def __init__ (self, calibration_filename, pulse_shaper_resolution, transform_limited_phase) :
		"""
		Load calibration information:
			<calibration_filename> -- path to calibration file, which is generated by the corresponding calibration program.
			<pulse_shaper_resolution> -- number of pixels in the shaper
		"""
		#super(AmplPhase2ShaperMasks_Basic, self).__init_()
		raise NotImplementedError("This code is depreciated")
			
		self.pulse_shaper_resolution = pulse_shaper_resolution
		
		# Loading the calibration file
		with h5py.File (calibration_filename, 'r') as calibration_file : 
			initial_pixel 			= calibration_file["settings/CalibrateShaper/initial_pixel"][...]
			final_pixel 			= calibration_file["settings/CalibrateShaper/final_pixel"][...]
			
			# load calibration information in tabulated form
			calibration_data = [
				(	int(key.split('_')[-1]), 
					pixel["voltage_master_mask"][...], 
					pixel["phase_master_mask"][...], 
					pixel["voltage_slave_mask"][...], 
					pixel["phase_slave_mask"][...],
					pixel["pixel_bound"][...]	) 
				for key, pixel in calibration_file["calibrated_pixels/pixels"].items() 
			]
			
			# Define the constant -- number of parameters, which is the number of logical pixels that can be manipulated
			self.num_parameters = len(calibration_data)	
			"""
			# Loading the transform limited phase
			try : 
				transform_limited_phase = calibration_file["transform_limited_phase"][...]
				if transform_limited_phase.size != self.num_parameters : 
					print "Warning: transform limited phase saved in the calibration file is invalid" 
					raise KeyError 
			except KeyError : transform_limited_phase = np.zeros(self.num_parameters)
			"""
		# sort by pixel number
		calibration_data.sort()
			
		# Unpacking data
		_, self.voltage_master_mask, self.phase_master_mask, self.voltage_slave_mask, \
			self.phase_slave_mask, self.pixels_edges = zip(*calibration_data)

		self.pixels_edges = np.array(self.pixels_edges)
		self.pixels_edges = np.vstack( (self.pixels_edges.min(axis=0), self.pixels_edges.max(axis=0)) )
		
		# Apply the transform limited phase correction
		transform_limited_phase *= 0.5
		for P, P_TL in izip(self.phase_master_mask, transform_limited_phase) : P += P_TL
		for P, P_TL in izip(self.phase_slave_mask, transform_limited_phase)  : P += P_TL
		
		def MinMaxPhase (voltage, phase) :
			""" Extract min and max phase that corresponds to voltages that are bellow maximal value """
			valid_indx = [ np.nonzero( (V > 0)&(V < PULSESHAPER_MAX_VAL) ) for V in voltage ]
			gen = ( P[indx].max() for P, indx in izip(phase,valid_indx) )
			max_p = np.fromiter( gen, np.float, len(phase) ) 
			gen = ( P[indx].min() for P, indx in izip(phase,valid_indx) )
			min_p = np.fromiter( gen, np.float, len(phase) ) 
			return min_p, max_p
	
		min_phase_master_mask, max_phase_master_mask = MinMaxPhase(self.voltage_master_mask, self.phase_master_mask)
		min_phase_slave_mask, max_phase_slave_mask = MinMaxPhase(self.voltage_slave_mask, self.phase_slave_mask)
	
	
		broken_pixels_indx = FindBrokenPixels(min_phase_master_mask, max_phase_master_mask, 
											min_phase_slave_mask, max_phase_slave_mask)
		if broken_pixels_indx.size :
			print "Pulse Shaper Calibration Warning: The following pixels may not be well calibrated: "
			print broken_pixels_indx
			print "\n"
		
		# Find range of accessible phases
		#self.phase_max = 2*np.minimum( max_phase_master_mask - np.pi/2., max_phase_slave_mask )
		#self.phase_min = 2*np.maximum( min_phase_master_mask, min_phase_slave_mask + np.pi/2.)
		self.phase_min, self.phase_max = GetTotalPhaseBounds (min_phase_master_mask, max_phase_master_mask, 
											min_phase_slave_mask, max_phase_slave_mask)
		
		
		# Calibration functions: Voltage of the mask as a function of phase
		self.voltage_function_master_mask = [ PchipInterpolator(*Order(P,V)) for P,V in izip(self.phase_master_mask,self.voltage_master_mask) ]
		self.voltage_function_slave_mask = [ PchipInterpolator(*Order(P,V)) for P,V in izip(self.phase_slave_mask,self.voltage_slave_mask) ]
	
	def __len__ (self) :
		"""
		Return number of logical pixels
		"""
		return self.num_parameters

	def ValidateAmplPhase (self, amplitude, phase, copy=True) :
		"""
		`amplitude` and `phase` array consistency check and convert `phase` into radians.
		`copy` - controls whether a duplicate of `amplitude` and `phase` needs to be created
		"""
		# make sure that the arguments are numpy arrays
		amplitude 	= np.array(amplitude, copy=copy)
		phase		= np.array(phase, copy=copy)
		
		if amplitude.size != self.num_parameters :
			print "PulseShaping Warning: The array of amplitudes has wrong size. Adjusting the size"
			amplitude = np.append( amplitude[:self.num_parameters], np.zeros(max(0,self.num_parameters-amplitude.size)) )
		
		if phase.size != self.num_parameters :
			print "PulseShaping Warning: The array of phases has wrong size. Adjusting the size"
			phase = np.append( phase[:self.num_parameters], np.zeros(max(0,self.num_parameters-phase.size)) )
		
		if amplitude.max() > 1. or amplitude.min() < 0. :
			print "PulseShaping Warning: Amplitude must take the values from 0 to 1. Renormalizing amplitude to enforce this requirement"
			amplitude -= amplitude.min(); amplitude /= amplitude.max()
		
		if phase.max() > 1. or phase.min() < 0. :
			print "PulseShaping Warning: Phase must take the values from 0 to 1. Renormalizing phase to enforce this requirement"
			phase -= phase.min(); phase /= phase.max()
		
		# Scaling `phase` to fit the accessible range from `self.phase_min` to `self.phase_max`
		phase *= (self.phase_max - self.phase_min)
		phase += self.phase_min
		
		return amplitude, phase
	
	def __call__ (self, amplitude, phase, copy=True) :
		"""
		Convert amplitude and phase into the pulse shaper masks
		"""
		amplitude, phase = self.ValidateAmplPhase(amplitude, phase, copy)
		
		# Convert the phase amplitude information to the phases of the master and slave masks
		# by using the following equations (see the comments inside the constructor)
		#
		# 		amplitude = np.cos( phase_master_mask - phase_slave_mask  )**2
		# 		phase = phase_master_mask + phase_slave_mask
		
		tmp = np.arccos(np.sqrt(amplitude)) 
		phase_master_mask = 0.5*( phase + tmp )
		phase_slave_mask = 0.5*(  phase - tmp )
	
		# Creating empty masks
		master_mask = np.zeros(self.pulse_shaper_resolution, dtype=ShaperInt)
		slave_mask = np.zeros_like(master_mask)
		
		# Calculate voltages and check for consistency
		voltage_master_mask = np.abs( self.GetMasterMaskVoltage(phase_master_mask) )
		indx1 = np.nonzero(voltage_master_mask > PULSESHAPER_MAX_VAL)[0]
		voltage_master_mask[ indx1 ] = PULSESHAPER_MAX_VAL
		
		voltage_slave_mask	=  np.abs( self.GetSlaveMaskVoltage(phase_slave_mask) )
		indx2 = np.nonzero(voltage_slave_mask > PULSESHAPER_MAX_VAL)[0]
		voltage_slave_mask[ indx2 ] = PULSESHAPER_MAX_VAL
		
		if indx1.size or indx2.size :
			print "\nPulse Shaping Warning: Voltages at the following pixels is out of range: "
			print np.unique( np.append(indx1, indx2) )
		
		# Filling the masks
		for bound, MV, SV in izip(self.pixels_edges, voltage_master_mask, voltage_slave_mask ) :
			master_mask[ bound[0]:bound[1] ] 	= MV
			slave_mask[ bound[0]:bound[1] ] 	= SV
			
		return master_mask, slave_mask
			
	def GetMasterMaskVoltage (self, phase_master_mask) :
		return np.fromiter( 
				( V(Phi) for V,Phi in izip(self.voltage_function_master_mask, phase_master_mask) ), 
			ShaperInt, len(phase_master_mask) )
	
	def GetSlaveMaskVoltage (self, phase_slave_mask) :
		return np.fromiter( 
				( V(Phi) for V,Phi in izip(self.voltage_function_slave_mask, phase_slave_mask) ), 
			ShaperInt, len(phase_slave_mask) )
		
class AmplPhase2ShaperMasks_SurfaceCalibration :
	"""
	Class for converting amplitudes and phases into the pulse shaper masks 
	using surface calibration
	"""		
	def __init__ (self, calibration_filename, pulse_shaper_resolution, transform_limited_phase, ampl_correction) :
		"""
		Load calibration information:
			<calibration_filename> -- path to calibration file, which is generated by the corresponding calibration program.
			<pulse_shaper_resolution> -- number of pixels in the shaper
			<transform_limited_phase> -- additive correction to the phase to achieve the transform limited pulse
			<ampl_correction>	-- multiplicative amplitude correction to achieve spectrally rectangular pulse shapes
		"""		
		self.pulse_shaper_resolution = pulse_shaper_resolution
		
		# Loading the calibration file
		with h5py.File (calibration_filename, 'r') as calibration_file : 
			self.initial_pixel 			= int(calibration_file["settings/CalibrateShaper/initial_pixel"][...])
			self.final_pixel 			= int(calibration_file["settings/CalibrateShaper/final_pixel"][...])
			
			# Number of parameters to vary
			self.num_parameters 		= self.final_pixel - self.initial_pixel
			
			#################### Master mask settings #########################
			master_mask_data			= calibration_file["calibrated_surface/master_mask"]
			
			self.master_mask_multiplier = master_mask_data["multiplier"][...]
			self.master_mask_offset		= master_mask_data["offset"][...]
			self.master_mask_calibration = PchipInterpolator( 	
								*Order(master_mask_data["calibration_curve_phase"][...],
								master_mask_data["calibration_curve_voltage"][...])
							)
			# The inverse function to self.master_mask_calibration
			self.voltage2phase_master_mask  = PchipInterpolator( 	
								*Order(master_mask_data["calibration_curve_voltage"][...],
								master_mask_data["calibration_curve_phase"][...])
							)
			
			#################### Slave mask settings #########################
			slave_mask_data		= calibration_file["calibrated_surface/slave_mask"]
			
			self.slave_mask_multiplier 	= slave_mask_data["multiplier"][...]
			self.slave_mask_offset		= slave_mask_data["offset"][...]
			self.slave_mask_calibration = PchipInterpolator( 	
								*Order(slave_mask_data["calibration_curve_phase"][...],
								slave_mask_data["calibration_curve_voltage"][...])
							)
			# The inverse function to self.slave_mask_calibration 
			self.voltage2phase_slave_mask = PchipInterpolator( 	
								*Order(slave_mask_data["calibration_curve_voltage"][...],
								slave_mask_data["calibration_curve_phase"][...])
							)
		
		# Verifications
		assert self.master_mask_offset.size == self.slave_mask_offset.size
		assert self.master_mask_multiplier.size == self.slave_mask_multiplier.size
		assert self.master_mask_multiplier.size == self.master_mask_offset.size
		if self.pulse_shaper_resolution != self.master_mask_multiplier.size :
			raise ValueError ("Calibration does not match pulse shaper. Recalibrate current pulse shaper")
			
		##################################################################
		
		@np.vectorize
		def TotalPhaseVariation (min_mm, max_mm, min_sm, max_sm, Condition) :
			"""
			Determinate the range of total phase variation with arbitrary amplitudes 
			in each pixel of pulse shaper 
			"""
			# Get phase ranges for both masks
			bins = 200
			phase_mm = np.linspace(min_mm, max_mm, bins)[:,np.newaxis]
			phase_sm = np.linspace(min_sm, max_sm, bins)[np.newaxis,:]
			
			diff_phase 	= phase_mm - phase_sm
			sum_phase 	= phase_mm + phase_sm
			
			# get histogram for total phase
			hist, bin_edges = np.histogram( 
				np.extract(Condition(diff_phase), sum_phase), bins=bins 
			)
	
			# Find incidences of most common values in histogram
			indx = np.nonzero( hist == Counter(hist).most_common(1)[0][0] )[0]

			# Return min and max phase
			return bin_edges[indx.min()+1], bin_edges[indx.max()-1]
	
		# Get bounds in phase variation
		min_phase_master_mask 	= self.master_mask_offset
		max_phase_master_mask 	= self.master_mask_offset + self.master_mask_multiplier
		min_phase_slave_mask	= self.slave_mask_offset
		max_phase_slave_mask	= self.slave_mask_offset + self.slave_mask_multiplier
		
		################## "plus" case ##################
		min_phase_plus, max_phase_plus = TotalPhaseVariation(
			min_phase_master_mask, max_phase_master_mask, 
			min_phase_slave_mask, max_phase_slave_mask,
			lambda diff_phase : (diff_phase >= 0)&(diff_phase <= 0.5*np.pi)
		)
		################### "minus" case ##################
		min_phase_minus, max_phase_minus = TotalPhaseVariation(
			min_phase_master_mask, max_phase_master_mask, 
			min_phase_slave_mask, max_phase_slave_mask,
			lambda diff_phase : (diff_phase <= 0)&(diff_phase >= -0.5*np.pi)
		)
		
		cond =( max_phase_minus - min_phase_minus > max_phase_plus - min_phase_plus ) 
		if cond.sum() > cond.size/2 :
			# On average, "minus" case gives more phase variation
			self.phase_min = min_phase_minus
			self.phase_max = max_phase_minus
			self.phase_sign = -1
		else :
			# On average, "plus" case gives more phase variation
			self.phase_min = min_phase_plus
			self.phase_max = max_phase_plus
			self.phase_sign = +1
	
		############### Limit phase variation to [0, 2*pi] ###############
		self.phase_min = self.phase_min.max()
		self.phase_max = self.phase_min + 2*np.pi
			
		################ Consistency check ####################
		# that the total phase chosen is reachable
		assert np.all(self.phase_min >= min_phase_master_mask + min_phase_slave_mask), \
			"Minimal total phase cannot be reached"
		assert np.all(self.phase_max <= max_phase_master_mask + max_phase_slave_mask), \
			"Maximal total phase cannot be reached"
		
		#######################################################
		# Set transform limited phase and amplitude correction
		if len(transform_limited_phase) :
			self.transform_limited_phase = self.ValidateArray(transform_limited_phase)
		else : self.transform_limited_phase = None
		
		if len(ampl_correction) :
			self.ampl_correction = self.ValidateArray(ampl_correction)
		else : self.ampl_correction = None
		
	def __len__ (self) :
		"""
		Return number of sencitive pixels
		"""
		return self.num_parameters
	
	def ValidateArray (self, A) :
		"""
		A combined function for checking the amplitude and phases masks
		"""
		if A.size != self.pulse_shaper_resolution :
			if A.size > self.num_parameters :
				raise ValueError ("PulseShaper Error: the size of amplitude or phase are too large")
				
			if A.size < self.num_parameters :
				# Stretching masks are needed because pixels are bundled
				x = np.linspace(0.,1.,self.num_parameters)
				xp = np.linspace(0.,1.,A.size)
				A 	= np.interp(x, xp, A)
					
			# Padding masks
			pad_width = ( self.initial_pixel, self.pulse_shaper_resolution-self.final_pixel )
			A = np.pad( A, pad_width, 'constant', constant_values=(A[0], A[-1]) )
				
		return A
			
	def ValidateAmplPhase (self, amplitude, phase, copy=True) :
		"""
		`amplitude` and `phase` array consistency check and convert `phase` into radians.
		`copy` - controls whether a duplicate of `amplitude` and `phase` needs to be created
		"""				
		##########################################################
		# make sure that the arguments are numpy arrays
		amplitude 	= np.array(amplitude, copy=copy)
		phase		= np.array(phase, copy=copy)
		
		# Check amplitude
		amplitude = self.ValidateArray(amplitude)
		# Check phase
		phase = self.ValidateArray(phase)
	
		# Add transform limited phase and apply wrapping
		if self.transform_limited_phase is not None :
			phase += self.transform_limited_phase
		
		# Multiply the amplitude correction
		if self.ampl_correction is not None :
			amplitude *= self.ampl_correction
		
		# Enforcing (0,1) value range
		np.clip(amplitude, 0, 1, out=amplitude)  
		phase %= 1.
	
		# Scaling `phase` to fit the accessible range from `self.phase_min` to `self.phase_max`
		phase *= (self.phase_max - self.phase_min)
		phase += self.phase_min
	
		return amplitude, phase
		
	def GetUnwrappedAmplPhase (self, amplitude, phase) :
		"""
		Get actual values of phases (in rads) and amplitude (0, 1).
		This function is inverse to self.__call__
		"""
		# Get mask voltages 
		master_mask, slave_mask = self.__call__(amplitude, phase)
		
		# Obtain corresponding slave phases 
		phase_slave_mask = self.voltage2phase_slave_mask(slave_mask)
		phase_slave_mask *= self.slave_mask_multiplier
		phase_slave_mask += self.slave_mask_offset
		
		# Obtain corresponding master phases 
		phase_master_mask = self.voltage2phase_master_mask(master_mask)
		phase_master_mask *= self.master_mask_multiplier
		phase_master_mask += self.master_mask_offset
		
		# Recover amplitude and phase 
		amplitude = np.cos( phase_master_mask - phase_slave_mask  )**2
		phase = phase_master_mask + phase_slave_mask
		
		return amplitude, phase
		
	def __call__ (self, amplitude, phase, copy=True) :
		"""
		Convert amplitude and phase into the pulse shaper masks
		"""
		amplitude, phase = self.ValidateAmplPhase(amplitude, phase, copy)
		
		# Convert the phase amplitude information to the phases of the master and slave masks
		# by using the following equations (see the comments inside the constructor)
		#
		# 		amplitude = np.cos( phase_master_mask - phase_slave_mask  )**2
		# 		phase = phase_master_mask + phase_slave_mask
		tmp = np.arccos(np.sqrt(amplitude))
		tmp *= self.phase_sign
		phase_master_mask = 0.5*( phase + tmp )
		phase_slave_mask = 0.5*(  phase - tmp )
	
		# Get voltages for master mask
		phase_master_mask -= self.master_mask_offset
		phase_master_mask /= self.master_mask_multiplier
		np.clip(phase_master_mask, 0, 1, out=phase_master_mask)
		#if not np.all(phase_master_mask >= 0.) : print "Problem master < 0"
		#if not np.all(phase_master_mask <= 1.) : print "Problem master > 1"
		master_mask = self.master_mask_calibration(phase_master_mask)
		np.round(master_mask, out=master_mask)
		master_mask = master_mask.astype(ShaperInt)
		
		# Get voltages for slave mask
		phase_slave_mask -= self.slave_mask_offset
		phase_slave_mask /= self.slave_mask_multiplier
		#if not np.all(phase_slave_mask >= 0.) : print "Problem slave < 0"
		#if not np.all(phase_slave_mask <= 1.) : print "Problem slave  > 1"
		np.clip(phase_slave_mask, 0, 1, out=phase_slave_mask)
		slave_mask = self.slave_mask_calibration(phase_slave_mask)
		np.round(slave_mask, out=slave_mask)
		slave_mask = slave_mask.astype(ShaperInt)
		
		return master_mask, slave_mask
		
		
class AmplPhase2ShaperMasks :
	"""
	Class facilitating the selection of methods for converting amplitudes and phases into shaper masks
	"""
	@classmethod
	def Choices (cls) :
		return [ "surface calibration", "pixel-wise calibration" ]
	
	@classmethod
	def Initialize (cls, option, *args, **kwargs) :
		"""
		Return an instance of class specified by string `option`. 
		The method `cls.Choices()` returns list of valid values for `option`.
		"""
		class_choices = dict( zip(cls.Choices(), 
					[AmplPhase2ShaperMasks_SurfaceCalibration, AmplPhase2ShaperMasks_Pixelwise]
				) )
		return class_choices[ option ](*args, **kwargs)
	
########################################################################
#
#	Process where the shaper resides
#
########################################################################


class PulseShaper (BasicDevice):
	"""
	Control pulse shaper
	"""
	
	def GetPixelNumber (self, arguments=None) :
		"""
		Returns number of pixels of the pulse shaper
		"""
		return self.npixles
	
	def ValidateMask (self, mask) :
		"""
		Verify that <mask> is a valid mask and convert the mask to a string
		suitable for sending to the shaper
		"""
		#print mask
		#print "================================================================"
		
		
		if not isinstance(mask, np.ndarray) :
			raise TypeError ("ValidateMask Error: Mask must be a numpy array")
		
		if len(mask) != self.npixles :
			raise ValueError ("ValidateMask Error: Mask's size is %d, but should be %d" % (len(mask), self.npixles ))
		
		if mask.dtype != ShaperInt :
			raise TypeError ("ValidateMask Error: Mask must be of ShaperInt type")
		
		if ( mask > PULSESHAPER_MAX_VAL ).sum() != 0 :
			raise ValueError ("ValidateMask Error: Mask must not contain values higher than %d" % PULSESHAPER_MAX_VAL)
		
		return (chr(x) for x in np.frombuffer(mask.data, dtype=np.uint8))
			
	def SetMasks (self, arguments) :
		"""
		Set voltages for masks
		"""
		# Extracting and validating masks
		master_mask, slave_mask = arguments
		master_mask = self.ValidateMask(master_mask)
		slave_mask	= self.ValidateMask(slave_mask)
		
		# Select the master mask (only for dual mask systems)
		self.serial_port.write("M0\r")
		
		# Select the defying frame
		self.serial_port.write("F0\r")
		
		# Initializing the block transfer 
		self.serial_port.write("B1\r")
		
		# Transferring the block 
		self.serial_port.write( "".join(master_mask) )

		# Clean the buffer
		self.serial_port.readlines()		
		
		######## Check that the sent data to the pulse shaper coincides with returned ######
		#self.serial_port.write("B?\r")
		#retuned = "".join(self.serial_port.readlines())[3:]
		#if result == retuned : print "Send and returned masks are the same"
		#else : print "Error: Send and returned masks are different"
		############################################3
		
		# Repeating the same steps as above for the slave mask
		self.serial_port.write("M1\rF0\rB1\r")
		
		self.serial_port.write( "".join(slave_mask) )
		
		# Select active frame
		self.serial_port.write("P0\r")
		
		# Clean the buffer
		self.serial_port.readlines()
		
		# Waiting for updating
		time.sleep( 1.e-3*self.update_time_delay )
		
		return RETURN_SUCCESS
	
	def SetUniformMasks (self, arguments) :
		"""
		Apply uniform masks
		"""
		vol_master_mask, vol_slave_mask = arguments
		ones = np.ones(self.npixles, dtype=ShaperInt)
		return self.SetMasks( (vol_master_mask*ones, vol_slave_mask*ones) )
		
	def SetAmplPhase (self, arguments) :
		try : self.AmplPhase2ShaperMasks 
		except AttributeError :
			print "PulseShaping Warning: Calibration file was not loaded"
			return RETURN_FAIL
	
		return self.SetMasks( self.AmplPhase2ShaperMasks(*arguments, copy=False) ) 
			
	def GetParamNumber (self, arguments=None) :
		"""
		Return number of logical pixels of the shaper as a result of calibration 
		"""
		try : return len(self.AmplPhase2ShaperMasks)
		except AttributeError :
			print "PulseShaping Warning: Calibration file was not loaded"
			return RETURN_FAIL
	
	def GetUnwrappedAmplPhase (self, arguments) :
		"""
		Get actual values of phases (in rads) and amplitude (0, 1)
		"""
		try : self.AmplPhase2ShaperMasks 
		except AttributeError :
			print "PulseShaping Warning: Calibration file was not loaded"
			return RETURN_FAIL
			
		return  self.AmplPhase2ShaperMasks.GetUnwrappedAmplPhase(*arguments)
		
	def Initialize (self, settings) :
		# Close the port if it is already used
		try : del self.serial_port
		except AttributeError : pass 
	
		# Start the communication port
		# Timout is chosen 50ms
		self.serial_port = serial.Serial (port=settings["port"], baudrate=460800,  bytesize=8, 
			parity=serial.PARITY_NONE, stopbits=1, timeout=0.05)
		# Wrapper to define '\r' as the end of line
		#self.serial_port = io.BufferedRWPair(self.serial_port, self.serial_port)
		#self.serial_port = io.TextIOWrapper(self.serial_port, newline='\r', line_buffering=True)
	
		# get shaper info
		self.serial_port.write ("v?\r")
		#print self.serial_port.readlines()[-1][20:-2]
		info = self.serial_port.readlines()[0].split(' ')
		
		# Delete non numerical values
		def IsInt (x) :
			"""
			Is x convertible to int?
			"""
			try : int(x) 
			except ValueError : return False
			return True
		info = filter(IsInt, info)
		
		# Convert to numbers
		info = map(int, info)
		
		# Saving the characteristics
		kind = {0 : "Phase", 1 : "Amplitude", 2 : "Dual"}[ info.pop() ]	
		info = (kind,) + tuple(info) 
		self.kind, self.firmware_rev, self.hardware_rev, self.serial_number, self.npixles = info
		
		if self.kind != "Dual" : raise ValueError ("Pulse shaper error: This is not a dual pulse shaper")
		
		# Printing the characteristics
		print "\n%s pulse shaper initialized: firmware rev %d; hardware rev %d; serial number %d; # of pixels %d\n" % info

		# Saving time delay 
		self.update_time_delay = settings["update_time_delay"]
		
		# Get transform limited phase
		transform_limited_phase = np.array( map(np.float, settings["transform_limited_phase"].split() ) )
		
		# Get the amplitude correction
		ampl_correction = np.array( map(np.float, settings["ampl_correction"].split() ) )
		
		# Loading the calibration data
		try : 
			# Select the conversion class
			self.AmplPhase2ShaperMasks = AmplPhase2ShaperMasks.Initialize(  
				settings["amplphase2shapermasks"], settings["calibration_file_name"], 
				self.npixles, transform_limited_phase, ampl_correction
			)
		except (IOError, KeyError, ValueError), msg :
			print "Warning: Pulse shaper calibration file could not be loaded: %s" % msg  
			
		return RETURN_SUCCESS
	
	def StopDevice (self, arguments=None) :
		try : del self.serial_port
		except AttributeError : pass
		return RETURN_SUCCESS
		
		
########################################################################

class PulseShaperTab (HardwareGUIControl) :
	"""
	This class represents a GUI controlling properties of pulse shaper.
	"""
	def __init__(self, parent, dev) :
		HardwareGUIControl.__init__(self, parent, dev)
		
		sizer = wx.BoxSizer(wx.VERTICAL)
		
		# Specify the communication port name
		sizer.Add (wx.StaticText(self, label="Communication port"), flag=wx.LEFT, border=5)
		port_name = wx.TextCtrl (self, value="COM15")
		port_name.__label__ = "port"
		sizer.Add (port_name, flag=wx.EXPAND, border=5)
		
		# How long to wait after new pulse shape was sent
		sizer.Add (wx.StaticText(self, label="\nWait to update (ms)"), flag=wx.LEFT, border=5)
		update_time_delay_ctr = wx.SpinCtrl (self, value="50", min=0, max=10000)
		update_time_delay_ctr.__label__ = "update_time_delay"
		sizer.Add (update_time_delay_ctr, flag=wx.EXPAND, border=5)	
		
		# Calibration file
		sizer.Add (wx.StaticText(self, label="\nCalibration file"), flag=wx.LEFT, border=5)
		calibration_file_ctr = wx.FilePickerCtrl(self, message="Chose pulse shaper calibration file...")
		calibration_file_ctr.__label__ = "calibration_file_name"
		sizer.Add (calibration_file_ctr, flag=wx.EXPAND, border=5)		
		
		# Select how to utilize calibration information
		sizer.Add (wx.StaticText(self, label="\nAmplitude/phase to voltage conversion"), flag=wx.LEFT, border=5)
		choices = AmplPhase2ShaperMasks.Choices()
		amplphase_convertion_ctrl = wx.ComboBox(self, choices=choices, value=choices[0], style=wx.CB_READONLY )
		amplphase_convertion_ctrl.__label__ = "amplphase2shapermasks"
		sizer.Add (amplphase_convertion_ctrl, flag=wx.EXPAND, border=5)
		
		# Transform limited phase 
		sizer.Add (wx.StaticText(self, label="\nTrasform limited phase (multi-line text)"), flag=wx.LEFT, border=5)
		self.transform_lim_phase_ctrl = wx.TextCtrl (self, style=wx.TE_MULTILINE)
		self.transform_lim_phase_ctrl.__label__ = "transform_limited_phase"
		sizer.Add (self.transform_lim_phase_ctrl, flag=wx.EXPAND, border=5)
		
		# Button to load transform limited phase
		load_tl_ctrl = wx.Button (self, label="Load transform limited phase")
		load_tl_ctrl.Bind (wx.EVT_BUTTON, self.LoadTransformLimited)
		sizer.Add (load_tl_ctrl, flag=wx.EXPAND, border=5)
		
		# Amplitude correction
		sizer.Add (wx.StaticText(self, label="\nAmplitude correction (multi-line text)"), flag=wx.LEFT, border=5)
		self.ampl_correction_ctrl = wx.TextCtrl (self, style=wx.TE_MULTILINE)
		self.ampl_correction_ctrl.__label__ = "ampl_correction"
		sizer.Add (self.ampl_correction_ctrl, flag=wx.EXPAND, border=5)
		
		self.SetSizer(sizer)
		############### GUI is created, now generate settings ######################
		self.CreateSettingsDict()
	
	def LoadTransformLimited (self, event) :
		"""
		Loading transform limited phase from HDFD5 file
		"""
		dlg = LoadPulseShapesDialog (parent=self, title="Load transform limited phase from HDF5 file")
		dlg.ShowModal()
		tl_phase = dlg.GetLoadedData()
	
		if len(tl_phase) == 0 : return
		if len(tl_phase) > 1 :
			print "PulseShaper Warning: More then one pulse shape is selected, we will load the last one"
		
		tl_phase = tl_phase.pop()
		
		if str(tl_phase[0]) == "phase only" :
			self.transform_lim_phase_ctrl.SetValue( " ".join( str(_) for _ in tl_phase[-1] ) )
		else :
			print "PuseShaper Error: Transform limited pulse must be phase-only shaped"
		
		dlg.Destroy()