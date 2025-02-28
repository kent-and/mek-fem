from dolfin import * 

def boundary (x): return x[0] < DOLFIN_EPS or x[0] > 1 -DOLFIN_EPS 


Ns = [4, 8, 16, 32, 64, 128, 256]
Ps = [1, 2, 3, 4]
L2_errors = {} 

for P in Ps: 
    f = Expression("M_PI*M_PI*sin(M_PI*x[0])", degree=P+2)
    u_analytical = Expression("sin(M_PI*x[0])", degree=P+2)

    for N in Ns: 
        mesh = UnitIntervalMesh(N) 
        h = mesh.hmax()
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

        L2_errors[(P, N)] = (h, L2_error)  

# numbers to estimate the rates 
for P in Ps: 
    for i in range(len(Ns)-1):  
        print ("P ", P, L2_errors[(P, Ns[i])][1] / L2_errors[(P, Ns[i+1])][1])


# log-log plot for assessing whether the convergence rate is 
# polynomial and what the rate is 
from numpy import log
import matplotlib.pyplot as plt
for P in Ps: 
    hs = []
    errors = []
    for N in Ns: 
        h, error = L2_errors[(P, N)]
        hs.append(h)
        errors.append(error)
    plt.loglog(Ns, errors)
plt.legend(["Lagrange %d" % P for P in Ps])
plt.show()


