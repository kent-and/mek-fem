from scitools.std import cos, pi, linspace, arange
import matplotlib.pyplot as plt

N = 150
a = 0
b = 1
h = float(b-a)/(N+1)
k = arange(1, N+1)

x = k*h
u = cos(pi*x)
plt.plot(x,u, label=r"$\mu_k$ = $\cos$(k$\pi$h)")
plt.xlabel( r"kh",fontsize=20) 
plt.ylabel(r"$\mu_k$",fontsize=30)
plt.axis([0,1,-1,1])
plt.title('Jacobi',fontsize=20)
plt.legend()#fontsize=20)
plt.savefig("jacobi_plot_test.png")
plt.show()
raw_input("Hold!")
