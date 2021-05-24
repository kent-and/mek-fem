from dolfin import *
from numpy import linalg 

for D in [1, 2]: 
  for N in [4, 8, 16, 32]:
    if   D == 1:  mesh = UnitIntervalMesh(N)
    elif D == 2:  mesh = UnitSquareMesh(N, N)

    V = FunctionSpace(mesh, "Lagrange", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    a = u*v*dx  + inner(grad(u), grad(v))*dx  
    A = assemble(a) 
    e = linalg.eigvals(A.array()) 
    e.sort()
    c = e[-1] / e[0]

    print "D=%d, N=%3d, min eigenvalue=%5.3f, max eigenvalue=%5.3f, cond. number=%5.3f " % (D, N, e[0], e[-1], c) 


