import numpy, pylab
from math import pi
from numpy import linalg

def create_stiffness_matrix(N): 
  h = 1.0/(N+1)
  A = numpy.zeros([N,N])
  for i in range(N): 
    A[i,i] = 2.0/(h**2) 
    if i > 0: 
      A[i,i-1] = -1.0/(h**2) 
    if i < N-1: 
      A[i,i+1] = -1.0/(h**2) 
  A = numpy.matrix(A)
  return A 

N = 150 
x = numpy.arange(0, 1, 1.0/(N))
f = numpy.matrix(numpy.sin(pi*x)).transpose()
u_true = (1.0/(pi*pi))*numpy.matrix(numpy.sin(pi*x)).transpose()
u0 = numpy.matrix(numpy.random.random(N)).transpose()
A = create_stiffness_matrix(N)

eigenvalues = numpy.sort(linalg.eigvals(A))
mu_max = eigenvalues[-1]
mu_min = eigenvalues[0]

print "Highes eigenvalue", mu_max
print "Lowest eigenvalue", mu_min

def iterate(tau):
  u_prev = u0; u = u0
  pylab.plot(x,u, 'o')
  for i in range(1001):
    if i == 10 or i == 100 or i == 1000 or i == 1000:   
      pylab.plot(x,u)
    u = u_prev - tau*(A*u_prev - f)    
    u_prev = u 
  pylab.plot(x, u_true)
  pylab.xlabel('x')
  pylab.ylabel('u')
  pylab.legend(["$u^0$", "$u^{10}$", "$u^{100}$", "$u^{1000}$", "$u_{true}$"])
  pylab.show()

tau = 0.9/mu_max
iterate(tau)
tau = 2/(mu_max + mu_min)
iterate(tau)






