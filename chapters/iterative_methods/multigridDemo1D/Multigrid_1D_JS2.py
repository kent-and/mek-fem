from __future__ import division

import	numpy as np
import scipy.linalg as sp

def exact( x0, x1, N ):
	x = np.linspace( x0, x1, N + 1 )
	return np.sin(np.pi*x)

def A_matrix(N):
	"""
		Coefficient matrix
	"""
	
	c = np.zeros(N+1)
	c[1]=-1
	c[0]=2

	r=np.zeros(N+1)
	r[1]=-1

	A = sp.toeplitz(c,r)
	#FIXME this looks ugly
	A[0][0] = 1; A[-1][-1] = 1;
	A[0][1] = 0; A[1][0] = 0; 
	A[-1][-2] = 0; A[-2][-1] = 0; 
	return A
def f(x):
	"""
		source function
	"""
	return np.pi**2*np.sin(np.pi*x)
	
def b_vector( N, x0, x1 ):
	"""
		rhs vector
	"""
	h = (x1 - x0)/float(N)
	x = np.linspace( x0, x1, N + 1 )
	b = h**2*f(x)
	b[0], b[N] = 0, 0
	return b

def S(A, b, u0 = 0, n_iter = 100):
	"""
	Richardson iteration, Jacobi preconditioner
	
		A = FDM coefficient matrix N+1xN+1
		b = h**2*f
		u0 = initial guess (N+1 array)
		n = number of iterations, default = 100
	"""
	if type(u0) == int:
		u0 = np.random.rand(len(b))
	u_1 = u0	
	n = 0
#	print u_1.shape, A.shape
	while n < n_iter:
		n += 1
		D = np.diag(np.diag(A))
		u = u_1 - np.dot(sp.inv(D), np.dot(A, u_1) - b)
		u_1 = u
		
	return u

def res_op(u):
	"""
		restriction operator
	"""
	
	U = np.zeros( (len(u) - 1) / 2 + 1 )
	U = np.zeros(((len(u) - 2)/2) + 2)
	U[0] = u[0]; U[-1] = u[-1];#Keeping the boundarys
	for j in range( 1, len(U)-1):
  		U[j] = 1/4*( u[2*j - 1] + 2*u[2*j] + u[2*j + 1])
  	return U
	
def int_op(U):
	"""
		interpolator operator
	"""
	#print "\n"
	print "started interpolating"
	#print "\n"
	u = np.zeros( (len(U) - 2) * 2 + 2)
	u[0] = U[0]; u[-1] = U[-1];#Keeping the boundarys
	for j in range(1, len(U) - 1):
		u[2*j] = U[j]
	for j in range(0, len(U) - 1 ):
		u[2*j + 1] = 1/2*( U[j] + U[j + 1] )
	return u
	
def run(N, N_start, x0, x1, interpol, U = 0, max_grid = 16):
	"""
		N must be of size 2**k
	"""		
	A = A_matrix( N )
	b = b_vector( N, x0, x1 )
	#print "before smoothing: A, b, U", len(A), len(b)
	u = S(A, b, U)
	u[0]  = 0; u[-1] = 0 #FIXME, this is a quick fix
	#print "after smoothing"
#	print u
	run_it = False
	#print 'N is now:', N
	if not interpol:
		
		if N > max_grid:
			N = (N-1)/2 + 1
		  	#print u
			U = res_op( u )
			#print "N > max grid", N, U
			run_it = True
			interpol = False
			
		if N == max_grid + 1:
			#print " NNN", N
			A_next = A_matrix( (N+1)/2 )
			A = A_matrix( N )
			b = b_vector( N, x0, x1 )
			#print "At lowest level"
			#print res_op( np.dot( A, U ) - b )
			#print len(A), len(U)
			##print A_next
			#print '________________________________'
			#print sp.inv(A_next)
			#print '________________________________'
			#print len(int_op( np.dot(sp.inv( A_next ), res_op( np.dot( A, U ) - b )) ))
			u_solv = U - int_op( np.dot(sp.inv( A_next ), res_op( np.dot( A, U ) - b )) )
			u = u_solv
			#u = U # for kjoring uten solver
			#print "lowest level, N is:",N
			interpol = True
	  		
	if interpol:
		if N < N_start:
			N = (N-1)*2 + 1
			U = int_op( u )
			run_it = True
			
	if run_it:
		#print "calling myself, N:", N, 'U is',U
		u, A, b = run( N, N_start, x0, x1, interpol, U )
			
	return u, A, b
 	
if __name__ == '__main__':
 	u,A,b = run(129,129,0,1, False)
	print "-------------------------------------"
	print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
	print "-------------------------------------"
	#print u
	#print '-------'
	u_e = exact(0,1,129)
	error = abs(u-u_e)
	print u
	print u_e
	#print error
	#u,A,b = run(129, 129, 0, 1, False, U = u, max_grid = 16)
	#error = abs(u-u_e)
	#print error
	#print len(u), len(exact(0,1,17))

