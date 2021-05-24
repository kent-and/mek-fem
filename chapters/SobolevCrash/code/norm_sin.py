
from dolfin import *

N = 10000 
mesh = UnitInterval(N)
V = FunctionSpace(mesh, "Lagrange", 1)

for k in [1, 10, 100]: 
  u_ex = Expression("sin(k*pi*x[0])", k=k)
  u = project(u_ex, V)

  L2_norm = sqrt(assemble(u**2*dx))
  print "L2 norm of sin(%d pi x) %e " % (k, L2_norm) 

  L7_norm = pow(assemble(abs(u)**7*dx), 1.0/7)
  print "L7 norm of sin(%d pi x) %e " % (k, L7_norm) 

  H1_norm = sqrt(assemble(u*u*dx + inner(grad(u), grad(u))*dx ))
  print "H1 norm of sin(%d pi x) %e" % (k, H1_norm) 


