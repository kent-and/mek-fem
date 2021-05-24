
from dolfin import *

def boundary(x, on_boundary):
    return on_boundary

parameters["krylov_solver"]["relative_tolerance"] = 1.0e-18 
parameters["krylov_solver"]["absolute_tolerance"] = 1.0e-18 
parameters["krylov_solver"]["monitor_convergence"] = True 
parameters["krylov_solver"]["report"] = True 
#parameters["krylov_solver"]["maximum_iterations"] = 50000 
epss = [1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6]
data = {}
Ns= [64, 128, 256, 512, 1024]
#Ns= [8, 16, 32, 64]
for N in Ns:
    for eps in epss:  
	parameters["krylov_solver"]["relative_tolerance"] = eps 

	mesh = UnitSquareMesh(N, N)
	print " N ", N, " dofs ", mesh.num_vertices() 

	V = FunctionSpace(mesh, "P", 1)
	u = TrialFunction(V)
	v = TestFunction(V)

	u_ex = Expression("sin(3.14*x[0])*sin(3.14*x[1])", degree=3) 
	f = Expression("2*3.14*3.14*sin(3.14*x[0])*sin(3.14*x[1])", degree=3) 
	a = inner(grad(u), grad(v))*dx  
	L = f*v*dx 

	U = Function(V)

	A = assemble(a) 
	b = assemble(L)

	bc = DirichletBC(V, u_ex, boundary)
	bc.apply(A)
	bc.apply(b) 

        t0 = time()
	solve(A, U.vector(), b, "gmres", "amg")
        t1 = time()

        cpu_time = t1-t0

	error_L2 = errornorm(u_ex, U, 'L2', degree_rise=3)
	print N, eps, cpu_time, error_L2 
        data[(N, eps)] = (error_L2, cpu_time) 


for N in Ns:
    for eps in epss:  
      D1, D2 = data[(N, eps)]
      print " %3.1e (%3.1e) " % (D1, D2), 
    print ""
	
print "" 

for eps in epss:  
    for N in Ns:
      D1, D2 = data[(N, eps)]
      print " %3.1e (%3.1e) " % (D1, D2), 
    print ""
	
