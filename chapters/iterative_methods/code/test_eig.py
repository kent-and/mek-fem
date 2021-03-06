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

    print ("D=%d, N=%2d, min eig=%5.3f, max eig=%3.1f, cond=%3.1f " % (D, N, e[0], e[-1], c)) 


