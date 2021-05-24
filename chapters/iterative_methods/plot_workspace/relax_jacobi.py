from scitools.std import cos, pi, linspace, arange
import matplotlib.pyplot as plt


N = 150
a = 0
b = 1
h = float(b-a)/(N+1)
k = arange(1, N+1)

w=2./3

x = k*h
u = 1+w*(cos(pi*x)-1)
plt.plot(x,u, label=r"$\mu_k$=1+$\omega$($\cos$(k$\pi$h)-1)")

plt.xlabel( r"kh",fontsize=20) 
plt.ylabel(r"$\mu_k$($\frac{2}{3}$)",fontsize=30)
plt.axis([0,1,-1,1])
plt.title("Relaxed Jacobi",fontsize=20)
plt.legend()#fontsize=20)
plt.savefig("relax_jacobi_plot_test.png")
plt.show()
raw_input("hold!")
