
from dolfin import *
import numpy
from scipy import linalg, matrix



def boundary(x, on_boundary): return on_boundary

for N in [1000]: 
  mesh = UnitIntervalMesh(N)
  V = FunctionSpace(mesh, "Lagrange", 1)

  u = TrialFunction(V)
  v = TestFunction(V)


  bc = DirichletBC(V, Constant(0), boundary)
  A, _ = assemble_system(inner(grad(u), grad(v))*dx, Constant(0)*v*dx, bc)
  M, _ = assemble_system(u*v*dx, Constant(0)*v*dx, bc)

  AA = matrix(A.array())
  MM = matrix(M.array())

  k = numpy.arange(1, N, 1)
  eig = pi**2*k**2 

  l1, v  = linalg.eigh(AA)
  l2, v  = linalg.eigh(AA, MM)

  print "l1 min, max ", min(l1), max(l1) 
  print "l2 min, max ", min(l2), max(l2) 
  print "eig min, max ", min(eig), max(eig) 

  import pylab 
  pylab.grid(True)
  pylab.gca().xaxis.grid(True, which='minor')
  pylab.gca().yaxis.grid(True, which='minor')
  fsize=35
  lw = 5 


  pylab.loglog(l1[2:], linewidth=lw)
  pylab.loglog(l2[2:], linewidth=lw)
  pylab.loglog(eig, linewidth=lw)
  pylab.legend(["eig(A)", "eig(A,M)", "cont. eig"], loc="upper left")

  pylab.show()




