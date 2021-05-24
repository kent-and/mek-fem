
from numpy import * 

def create_stiffness_matrix(N): 
  h = 1.0/(N-1)
  A = zeros([N,N])
  for i in range(N): 
    A[i,i] = 2.0/(h**2) 
    if i > 0: 
      A[i,i-1] = -1.0/(h**2) 
    if i < N-1: 
      A[i,i+1] = -1.0/(h**2) 
  A = matrix(A)
  return A 

eps = 1.0e-6 
Ns = [10, 20, 40, 80, 160, 320] 
for N in Ns: 
  A = create_stiffness_matrix(N)                # creating matrix
  x = arange(0, 1, 1.0/(N))
  f = matrix(sin(3.14*x)).transpose()           # right hand side
  u0 = matrix(random.random(N)).transpose()     # initial guess 
  u_prev = u0 
  r = A*u_prev - f                              # compute the residual 
  norm_of_residual = r.transpose()*r / len(r)   # check for norm of residual 
  tol = eps*norm_of_residual                # relative convergence estimate of 1.0e-6  

  eigenvalues = sort(linalg.eigvals(A))         # compute eigenvalues and tau 
  lambda_max, lambda_min = eigenvalues[-1],  eigenvalues[0]
  tau = 2/(lambda_max + lambda_min)
  K = lambda_max / lambda_min
  rho = (K -1) / (K+1)  
  print "lambda_max ", lambda_max, " lambda_min ", lambda_min, " K ", K, " rho ", rho, " eps ", eps


  estimated_iterations = int(log(eps) / log(rho)) 


  no_iterations= 0                        
  while norm_of_residual > tol: 
    r = A*u_prev - f                          # compute the residual 
    u = u_prev - tau*r                        # the Richardson iteration 
    u_prev = u                                 
    norm_of_residual = r.transpose()*r        # check for norm of residual 
    no_iterations+=1                          # count no iterations  

  print "N ", N, "estimated no. iterations ", estimated_iterations,  " no. of iterations ", no_iterations
    









