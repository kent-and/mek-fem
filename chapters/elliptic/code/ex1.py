
from numpy import * 

M = 10**2 
N = 10 
x = arange(0, 1, 1/M) 
f = x*x 
analytical_solution = - 1.0/12*x**4 + 1.0/12*x 

# compute basis functions 
dNs = []
Ns = []
for k in range(1, N): 
    Nk = sin(k*pi*x) 
    Ns.append(Nk) 
    dNk = k*pi*cos(k*pi*x)   # could here use a finite difference, but here analytically calculated  
    dNs.append(dNk) 


# compute matrix (ignore that the matrix is in fact diagonal) 
A = zeros((N-1, N-1)) 
dx = 1.0 / N 
for i, dNi in enumerate(dNs): 
    for j, dNj in enumerate(dNs):  
        A[i,j] = sum(dNi*dNj)*dx  

# compute b 
b = zeros((N-1, 1)) 
dx = 1.0 / N 
for i, Ni in enumerate(Ns): 
    b[i] = sum(f*Ni)*dx  
#print (b) 


sol = linalg.solve(A, b) 

print (sol) 

u = x*0
for i, Ni in enumerate(Ns): 
    u += sol[i] * Ni


import matplotlib.pyplot as plt 
plt.plot(u) 
plt.plot(analytical_solution) 
plt.show()




