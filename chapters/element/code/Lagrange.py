
from sympy import * 

def dirac(i,j): 
  if i == j: return 1 
  else: return 0 

k = 5  
x = Symbol("x") 
dx = 1/k 

basis = [x**j for j in range(0, k+1)] 
points = [j*dx for j in range(0, k+1)] 
dofs = symbols("a:%d"%(k+1))

print ("A basis for the polynomial space ", basis) 
print ("A set of points to use for computing the nodal basis ", points) 
print ("The degrees of freedom, represented as symbols ", dofs) 

pol = sum([dofs[i]*basis[i] for i in range(0, (k+1))])   
print ("A generic function is then ", pol)

Lagrange_basis = []
for i in range(0, (k+1)): 
  equations = []
  for j in range(len(points)):  
    p = points[j]
    eq = pol.subs(x,p) - dirac(i,j) 
    equations.append(eq)
  coeff = solve(equations, dofs) 
  Nj = pol.subs(coeff)
  Lagrange_basis.append(Nj) 

print (Lagrange_basis) 

# check that the Lagrange basis is what it is supposed to be 
for i, basis in enumerate(Lagrange_basis): 
    for j, point in enumerate(points): 
        print ("i,j ", i,j, " basis evaluation ",  "%0.3f" % basis.subs({ x: point})) 



