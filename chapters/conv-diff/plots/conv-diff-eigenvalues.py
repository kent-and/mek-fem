
from dolfin import *

N = 40 
  
mesh = UnitInterval(N)
V = FunctionSpace(mesh, "CG", 1)

u = TrialFunction(V)
v = TestFunction(V)

alpha_value = 1.0e-5 
alpha = Constant(alpha_value)

beta_value = 0.5 
beta = Constant(beta_value)

f = Expression("x[0]")

f = Constant(0)
h = mesh.hmin()

a1 = -u.dx(0)*v*dx 
a2 = u.dx(0)*v.dx(0)*dx  
a3 = a1 + alpha*a2   
a4 = a3 + beta*h*u.dx(0)*v.dx(0)*dx  
a5 = -f*u.dx(0)*v*dx + alpha*a2  


A1 = assemble(a1) 
A2 = assemble(a2) 
A3 = assemble(a3) 
A4 = assemble(a4) 
A5 = assemble(a5) 

import numpy
from numpy import linalg 
e1 = linalg.eigvals(A1.array())
e2 = linalg.eigvals(A2.array())
e3 = linalg.eigvals(A3.array())
e4 = linalg.eigvals(A4.array())
e5 = linalg.eigvals(A5.array())

import pylab 
pylab.plot(numpy.sort(e1.real))
pylab.plot(numpy.sort(e1.imag))
pylab.legend(["Real part u_x", "Imag part u_x"])
pylab.show()

print  "Smallest eigenvalue ", numpy.sort(numpy.abs(e1.imag))[0] 
print  "Next smallest eigenvalue ", numpy.sort(numpy.abs(e1.imag))[1] 
print  "Largest eigenvalue ", numpy.sort(numpy.abs(e1.imag))[-1]

pylab.plot(numpy.sort(e2.real))
pylab.plot(numpy.sort(e2.imag))
pylab.legend(["Real part u_xx", "Imag part u_xx"])
pylab.show()

print  "Smallest eigenvalue ", numpy.sort(numpy.abs(e2.real))[0] 
print  "Next smallest eigenvalue ", numpy.sort(numpy.abs(e2.real))[1] 
print  "Largest eigenvalue ", numpy.sort(numpy.abs(e2.real))[-1]

print "e3 ", e3
pylab.plot(e3.real,e3.imag, "o")
pylab.legend(["Eigenvalues of -u_x-u_xx"])

print  "Smallest eigenvalue ", numpy.sort(numpy.abs(e3))[0] 
print  "Next smallest eigenvalue ", numpy.sort(numpy.abs(e3))[1] 
print  "Largest eigenvalue ", numpy.sort(numpy.abs(e3))[-1]

pylab.show()

print "e4 ", e4
pylab.plot(e4.real, e4.imag, "o")
pylab.legend(["Real part u_xx", "Imag part u_xx"])
pylab.legend(["Eigenvalues of -u_x-u_xx stab"])
pylab.show()

print  "Smallest eigenvalue ", numpy.sort(numpy.abs(e4))[0] 
print  "Next smallest eigenvalue ", numpy.sort(numpy.abs(e4))[1] 
print  "Largest eigenvalue ", numpy.sort(numpy.abs(e4))[-1]

print "e5 ", e5
pylab.plot(e5.real,e5.imag, "o")
pylab.legend(["Eigenvalues of part -x u_x-u_xx"])

print  "Smallest eigenvalue ", numpy.sort(numpy.abs(e5))[0] 
print  "Next smallest eigenvalue ", numpy.sort(numpy.abs(e5))[1] 
print  "Largest eigenvalue ", numpy.sort(numpy.abs(e5))[-1]

pylab.show()





