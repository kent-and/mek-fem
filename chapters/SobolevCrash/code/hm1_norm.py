
from dolfin import *
from numpy import matrix, diagflat, sqrt
from scipy import linalg 

def boundary(x, on_boundary): return on_boundary

mesh = UnitIntervalMesh(200)
V = FunctionSpace(mesh, "Lagrange", 1)

u = TrialFunction(V)
v = TestFunction(V)

bc = DirichletBC(V, Constant(0), boundary)
A, _ = assemble_system(inner(grad(u), grad(v))*dx, Constant(0)*v*dx, bc)
M, _ = assemble_system(u*v*dx, Constant(0)*v*dx, bc)

AA = matrix(A.array())
MM = matrix(M.array())

l, v = linalg.eigh(AA, MM)
v = matrix(v)
l = matrix(diagflat(l))
for k in [10]: 
#for k in [1, 10, 100]: 
  u_ex = Expression("sin(k*pi*x[0])", k=k)
  u = interpolate(u_ex, V)
  x = matrix(u.vector().array())

  H1_norm = pi*k*sqrt(2)/2  
  print "H1 norm of sin(%d pi x) %e (exact)          " % (k, H1_norm) 
  H1_norm = sqrt(assemble(inner(grad(u), grad(u))*dx)) 
  print "H1 norm of sin(%d pi x) %e (|grad(u)|^2)    " % (k, H1_norm) 
  H1_norm = sqrt(x*AA*x.T)    
  print "H1 norm of sin(%d pi x) %e (x A x' )        " % (k, H1_norm) 
  W = MM.dot(v)
  H1_norm = sqrt(x*W*l*W.T*x.T)   
  print "H1 norm of sin(%d pi x) %e (eig)            " % (k, H1_norm) 

  print "" 

  L2_norm = sqrt(2)/2 
  print "L2 norm of sin(%d pi x) %e (exact)          " % (k, L2_norm) 
  L2_norm = sqrt(assemble(u**2*dx)) 
  print "L2 norm of sin(%d pi x) %e   |u|^2          " % (k, L2_norm) 
  L2_norm = sqrt(x*MM*x.T) 
  print "L1 norm of sin(%d pi x) %e (x M x' )        " % (k, L2_norm) 
  W = MM.dot(v)
  L2_norm = sqrt(x*W*l**0*W.T*x.T)   
  print "L2 norm of sin(%d pi x) %e (eig)            " % (k, L2_norm) 


  print "" 

  Hm1_norm = sqrt(2)/2/k/pi  
  print "H^-1 norm of sin(%d pi x) %e (exact)        " % (k, Hm1_norm) 
  Hm1_norm = sqrt(x*W*l**-1*W.T*x.T)  
  print "H^-1 norm of sin(%d pi x) %e (eig)          " % (k, Hm1_norm) 
  Hm1_norm = sqrt(x*MM*linalg.inv(AA)*MM*x.T)    
  print "H^-1 norm of sin(%d pi x) %e (x inv(A) x') " % (k, Hm1_norm) 






