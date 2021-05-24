
from dolfin import *

for N in [10, 100]:

  mesh = UnitIntervalMesh(N)
  V = FunctionSpace(mesh, "CG", 1)

  u = TrialFunction(V)
  v = TestFunction(V)

  alpha_value = 1.0e-2 
  alpha = Constant(alpha_value)
  f = Constant(0)
  h = mesh.hmin()

  a = (-u.dx(0)*v + alpha*u.dx(0)*v.dx(0))*dx  
  L = f*v*dx  

  u_analytical = Expression("(exp(-x[0]/%e) - 1)/ (exp(-1/%e) - 1)" % (alpha_value, alpha_value), degree=3)  
  def boundary(x):
      return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

  bc = DirichletBC(V, u_analytical, boundary) 

  U = Function(V)
  solve(a == L, U, bc) 

  U_analytical = project(u_analytical, V)

  import pylab 
  pylab.plot(U.vector().array(), linewidth=4)
  pylab.plot(U_analytical.vector().array(), linewidth=4)
  pylab.legend(["Numerical Solution", "Analytical Solution"], fontsize=20)
  pylab.show()



