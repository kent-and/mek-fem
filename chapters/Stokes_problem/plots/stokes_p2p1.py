
from dolfin import *


def u_boundary(x):
  return x[0] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS 

def p_boundary(x):
  return  x[0] > 1.0 - DOLFIN_EPS  



mesh = UnitSquare(40,40) 
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)
W = MixedFunctionSpace([V, Q])

u, p = TrialFunctions(W)
v, q = TestFunctions(W)

f = Constant([0,0])

u_analytical = Expression(["x[1]*(1-x[1])", "0.0"])
p_analytical = Expression("-2+2*x[0]")

bc_u = DirichletBC(W.sub(0), u_analytical, u_boundary)
bc = [bc_u]

a = inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx 
L = inner(f, v)*dx 

UP = Function(W)
solve( a == L, UP, bc) 

print UP.vector().max()
print UP.vector().min()

U, P = UP.split() 




plot(U, title="Numerical velocity")
plot(P, title="Numerical pressure")

U_analytical = project(u_analytical, V)
P_analytical = project(p_analytical, Q)
plot(U_analytical, title="Analytical velocity")
plot(P_analytical, title="Analytical pressure")
interactive()


