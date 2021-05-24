from __future__ import division

import	numpy as np
import scipy.linalg as sp
import matplotlib.pyplot as plt

def exact( x0, x1, N ):
	x = np.linspace( x0, x1, N + 1 )
	return np.sin(np.pi*x)

def A_matrix(N):
	"""
		Coefficient matrix
	"""
	c = np.zeros( N + 1 )
	c[1]=-1
	c[0]=2

	r=np.zeros( N + 1 )
	r[1]=-1
	A = sp.toeplitz(c,r)
	A[0,:] = np.zeros( N + 1 )
	A[-1,:] = np.zeros( N + 1 )
	A[0,0] = A[-1,-1] = 1
#	A[0,1] = A[-1,-2] = - 1
	A[1,0] = A[-2,-1] = 0
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
	h = (x1 - x0)/(N)
	x = np.linspace( x0, x1, N + 1 )
	b = h**2*f(x)
	b[0], b[N] = 0, 0
	return b

def S(A, b, u0 = 0, n_iter = 250):
	"""
	Richardson iteration, Jacobi preconditioner
	
		A = FDM coefficient matrix N+1xN+1
		b = h**2*f
		u0 = initial guess (N+1 array)
		n = number of iterations, default = 100
	"""
	if type(u0) == int:
		u0 = np.random.rand(len(b))
#		u0 = np.ones_like(b)
	u_1 = u0	
	n = 0
#	print u_1.shape, A.shape, b.shape
	while n < n_iter:
		n += 1
		D = np.diag(np.diag(A))
		u = u_1 - np.dot(sp.inv(D), np.dot(A, u_1) - b)
		u_1 = u
	u[0] = u[-1] = 0
	return u

def res_op(u):
	"""
		restriction operator
	"""
	U = np.zeros( (len(u) - 1) / 2 + 1 )
	for j in range( 1, len(U) -1):
  		U[j] = 1/4*( u[2*j - 1] + 2*u[2*j] + u[2*j + 1] )
  	U[0] = u[0]; U[-1] = u[-1];#Keeping the boundarys
  	return U
	
def int_op(U):
	"""
		interpolator operator
	"""
	u = np.zeros( (len(U) - 1) * 2 + 1 )
	for j in range( len(U) ):
		u[2*j] = U[j]
	for j in range( 1, len(U) ):
		u[2*j - 1] = 1/2*( U[j] + U[j - 1] )
	return u
	
def run(N, N_start, x0, x1, interpol, U = 0, max_grid = 16):
	"""
		N must be of size 2**k
	"""		
	A = A_matrix( N )
	b = b_vector( N, x0, x1 )
	u = S(A, b, U)
	run_it = False
	if not interpol:
		
		if N > max_grid:
			N /= 2
		  	U = res_op( u )
			run_it = True
			interpol = False
			
		if N == max_grid:
			A_H = A_matrix( N/2 )
			A = A_matrix( N )
			b = b_vector( N, x0, x1 )
#			print U.shape, A.shape, A_H.shape, A.shape, b.shape
			u_solv = U - int_op( np.dot(sp.inv( A_H ), res_op( np.dot( A, U ) - b )) )
			u = u_solv
#			u = U # for kjoring uten solver
			interpol = True
	  		
	if interpol:
		if N < N_start:
			N *= 2
			U = int_op( u )
			run_it = True
			
	if run_it:
		u, A, b = run( N, N_start, x0, x1, interpol, U )
			
	return u, A, b
 	
if __name__ == '__main__':
	n2 = 128*2*2
 	u,A,b = run(n2,n2,0,1,False)
 	ue = exact(0,1,n2)
 	ut = np.dot(sp.inv(A),b)
 	#plt.plot(ut, label = 'test')
 	#plt.plot(ue, label = 'exact')
 	#plt.plot(u, label = 'iteration')
 	#plt.legend()
 	#plt.show()
 	#print A,'\n',b
	h = 1./n2
	error= np.max(np.abs(u-ue))
	print error, error/h, error/(h*h)
	

#	uv = np.array([0,1,2,1,2,1,2,1,0])
#	ures = res_op(uv)
#	
#	uint = int_op(ures)
#	print uint


