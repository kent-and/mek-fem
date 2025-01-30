
from dolfin import * 

N = 10 # number of elements 
P = 2  # order of the polynomial 


def boundary (x): return x[0] < DOLFIN_EPS or x[0] > 1 -DOLFIN_EPS 
f = Expression("M_PI*M_PI*sin(M_PI*x[0])", degree=P+2)
u_analytical = Expression("sin(M_PI*x[0])", degree=P+2)



Ns = [4, 8, 16, 32, 64]
Ps = [1, 2, 3, 4]
L2_errors = {} 

for P in Ps: 
    for N in Ns: 
        mesh = UnitIntervalMesh(N) 
        V = FunctionSpace(mesh, "Lagrange", P) 
        u = TrialFunction(V) 
        v = TestFunction(V) 

        # bc 
        u0 = Constant(0)
        bc = DirichletBC(V, u0, boundary) 

        a = inner(grad(u), grad(v))*dx  
        L = f*v*dx 

        U = Function(V) 
        solve (a ==L, U, bc) 

        L2_error = assemble(pow(U-u_analytical, 2)*dx) 
        L2_error = sqrt(L2_error) 

        L2_errors[(P, N)] = L2_error  
        print (N, P, L2_error) 

import matplotlib.pyplot as plt
for P in Ps: 
    errors = []
    for N in Ns: 
        errors.append(L2_errors[(P, N)])
    plt.loglog(errors)
plt.show()


