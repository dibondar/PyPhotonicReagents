"""
Assume we want to solve a linear system of equations:
	Ax = b
where A and b are contaminated by multiplicative noise. 

How does the noise influence the accuracy of the solution?
"""

import numpy as np
import scipy.stats as stats
from scipy.optimize import nnls
from numpy.linalg import lstsq
from functools import partial

class CAnalysis :
	
	def __init__ (self, MulNoiseGen, LinEqsSolver) :
		# Save the function that generates the multiplicative noise
		self.MulNoiseGen = MulNoiseGen
		
		# Save the solver of linear equations
		self.LinEqsSolver = LinEqsSolver

	def GetXExact (self) :
		# Generate solution x
		x_exact = np.random.rand( self.A_exact.shape[0] )
		
		# x_exact should imitate the concentrations, i.e.,
		# 	x_exat > 0 and sum(x_exact) = 1
		x_exact /= x_exact.sum()
		
		return x_exact
		
	def SingleSampleAnalysis (self, x_exact) :
		"""
		Draw one sample of noise and perform the analysis
		"""
		b_exact = self.A_exact.dot( x_exact )
		
		# Generate noisy linear equation
		noise = self.MulNoiseGen( size=self.A_exact.size ).reshape(self.A_exact.shape)
		A_noisy = self.A_exact * ( 1. + noise )
		# Enforce positivity 
		np.abs(A_noisy, out=A_noisy)
		
		noise = self.MulNoiseGen( size=b_exact.size )
		b_noisy = b_exact *( 1. + noise )
		# Enforce positivity 
		np.abs(b_noisy, out=b_noisy)
		
		# Solve the noisy linear equation 
		x_noisy = self.LinEqsSolver( A_noisy, b_noisy )[0]
		
		return np.abs( x_exact - x_noisy ).max()
		
	def __call__ (self, A_exact, x_sample_size, noise_sample_size) :
		"""
		A_exact 			-- exact matrix of the linear equations
		x_sample_size 		-- how many solutions x to be drawn
		noise_sample_size 	-- how many revelations of the multiplicative noise to be drawn
		"""
		self.A_exact = A_exact
		
		gen_x_exact = ( self.GetXExact() for _ in xrange(x_sample_size) )
		
		stats = np.vstack( [
			np.fromiter( (self.SingleSampleAnalysis(x) for _ in xrange(noise_sample_size)), np.float) 
				for x in gen_x_exact 
		] )
		
		return stats
		
		
if __name__ == '__main__' :
	
	Analys = CAnalysis( partial(np.random.normal, scale=0.1),  nnls )
	Analys = partial(Analys, x_sample_size=1, noise_sample_size=30000)
	
	#multiprocessing.freeze_support()
	#pool = multiprocessing.Pool( processes=(multiprocessing.cpu_count()-1) )
	
	#dim_dependence = pool.map(analys, [ np.random.rand(N,N) for N in [2,3,4,5] ] )
	
	import matplotlib.pyplot as plt
	
	for S, N in zip( [221,222,223,224], [2,3,4,5] ) :
		plt.subplot(S)
		r = np.ravel(Analys(np.random.rand(N,N)))
		stats.probplot(r, dist=stats.cauchy, plot=plt)
		#plt.hist(r, bins=100 )
		plt.title("Dimension %d" % N)
		
		print "STD for N = %d is %.2e" % (N, np.std(r))
		#print stats.anderson(r, dist='gumbel')
	plt.show()
	
	
	