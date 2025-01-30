from sympy import * 

x = Symbol("x") 

N = 4  

# compute basis functions 
# basis functions are of the form x^k * (1-x)^(N-k) 
# However, we notice that not all of them will statisfy the boundary condition   
# x^N is zero in 0 but 1 in 1. Hence, the function does not satisfy the bc  
# likewise (1-x)^N is zero in 1 but 1 in zero 
# Hence, these two need to be removed 
dNs = []
Ns = []
for k in range(1, N): 
    Nk = x**k * (1-x)**(N-k)  
    Ns.append(Nk) 

print(Ns) 

A = zeros(N-1, N-1) 
for i, Ni in enumerate(Ns): 
    dNi = diff(Ni, x) # compute the derivative of the basis function  
    for j, Nj in enumerate(Ns): 
        dNj = diff(Nj, x) 
        A[i,j] = integrate(dNi*dNj, (x, 0,1))  # perform the integrals  

# sympy deals represent the integrals as exact rational numbers 
print (A) 
# but can be evaluated to floats for more efficient linear solves 
print (A.evalf()) 
