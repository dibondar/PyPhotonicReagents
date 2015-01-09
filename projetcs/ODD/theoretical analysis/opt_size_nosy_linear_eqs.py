"""
Optimal dimension for noisy linear equations to get accurate concentrations
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import lstsq
from scipy.optimize import nnls
from functools import partial

########################################################
#
#	Parameters
#
########################################################

#np.random.seed(121)

# number of unknowns 
N_min, N_max = 2, 6

# number of equations
M_min, M_max = 2, 100

# multiplicative noise
MulNoiseGen = partial(np.random.normal, scale=0.005)

# Linear equations solver
LinEqsSolver = nnls

########################################################
#
#	
#
########################################################

# Generate unknowns 
x_exact = np.abs(np.random.rand(N_max))
x_exact /= x_exact.sum()

# Generate matrix
A_exact = np.abs( np.random.rand(M_max, N_max) )

# Find the rhs of the linear equations
b_exact = A_exact.dot( x_exact )


########################################################
#
#	Analysis
#
########################################################

errors = np.zeros( (N_max-N_min+1, M_max-M_min+1), dtype=np.float )

for n in np.arange(N_min, N_max+1) :
	for m in np.arange(M_min, M_max+1) :
		# Selecting subsystem with 
		# `n` unknowns and `m` equations 
		A = np.copy( A_exact[:m, :n] )
		b = np.copy( b_exact[:m] )
		
		# Contaminating system with the multiplicative noise
		A *= (1. + MulNoiseGen(size=A.size).reshape(A.shape) )
		b *= (1. + MulNoiseGen(size=b.size).reshape(b.shape) )
		#A += MulNoiseGen(size=A.size).reshape(A.shape) 
		#b += MulNoiseGen(size=b.size).reshape(b.shape) 
		
		# Enforce positivity
		A = np.abs(A)
		b = np.abs(b)
		
		# Solve noisy linear equation
		x = LinEqsSolver(A, b)[0]
		
		errors[n-N_min, m-M_min] =  np.max( np.abs( x - x_exact[:n] ) / x_exact[:n] )

print ( np.abs(LinEqsSolver(A_exact, b_exact)[0] - x_exact)/x_exact ).max()
########################################################
#
#	Plotting
#
########################################################

plt.subplot(211)
plt.imshow( np.log(errors), extent=(M_min, M_max, N_min, N_max), interpolation ='nearest', origin='lower' )
#plt.imshow( np.log10(errors),  )
plt.xlabel("Number of equations")
plt.ylabel("Number of unknowns")
plt.title("Log of maximum relative error")
plt.colorbar()

plt.subplot(212)
plt.plot( np.arange(M_min, M_max+1), np.log10(errors[-1,:]) )
plt.title ("Log of maximum relative error for *** unknowns")

print errors.mean()

plt.show()
		