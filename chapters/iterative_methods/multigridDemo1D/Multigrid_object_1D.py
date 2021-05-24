from __future__ import division
'''
Multigrid 1D demo

Description: Solves a 1D Poisson problem with the multigrid method
'''
__author__ = 'Mikkel Lepperod'
__email__ = 'lep.mik@gmail.com'

import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg


class Int_op(object):
    """
        interpolation operator
    """
    def __mul__(self, var):
        u = np.zeros((len(var) - 1) * 2 + 1)
        for j in range(len(var)):
            u[2 * j] = var[j]
        for j in range(1, len(var)):
            u[2 * j - 1] = 1 / 2 * (var[j] + var[j - 1])
        u = np.matrix(u).transpose()
        return u


class Res_op(object):
    """
        restriction operator
    """
    def __mul__(self, var):
        U = np.zeros((len(var) - 1) / 2 + 1)
        for j in range(1, len(U) - 1):
            U[j] = 1 / 4 * (var[2 * j - 1] + 2 * var[2 * j] + var[2 * j + 1])
        U[0] = var[0]
        U[-1] = var[-1]
        U = np.matrix(U).transpose()
        return U


class Func_init(object):
    """
        function initializer
    """
    def __init__(self, N, x0, x1):
        self.N = N
        self.x0 = x0
        self.x1 = x1
        self.h = (x1 - x0) / (self.N)

    def x(self, N):
        return np.linspace(self.x0, self.x1, N + 1)

    def A(self):
        """
            Coefficient matrix
        """
        c = np.zeros(self.N + 1)
        c[1] = -1
        c[0] = 2
        r = np.zeros(self.N + 1)
        r[1] = -1
        A = linalg.toeplitz(c, r)
        A[0, :] = np.zeros(self.N + 1)
        A[-1, :] = np.zeros(self.N + 1)
        A[0, 0] = A[-1, -1] = 1
        A[1, 0] = A[-2, -1] = 0
        return (1 / (self.h) ** 2) * np.matrix(A)

    def f(self, x):
        """
            source function
        """
        return np.pi ** 2 * np.sin(np.pi * x)

    def b(self):
        """
            rhs vector
        """
        x = np.linspace(self.x0, self.x1, self.N + 1)
        b = self.f(x)
        b[0], b[-1] = 0, 0
        return np.matrix(b).transpose()

    def exact(self):
        """
            calculates the exact solution
        """
        x = np.linspace(self.x0, self.x1, self.N + 1)
        self.u_ex = np.matrix(np.sin(np.pi * x)).transpose()


class Smoother(object):
    """
        smoothing operator, Richardson iteration with and without
        different preconditioners which are initialized with the
        following commands: (relaxed Jacobi is default) iteration =
        'richard_large','richard_small', 'jacobi','gauss_seidel',
        'relaxed_jacobi'
    """
    def __init__(self, n_iter, A, b, iteration=False):
        self.n_iter = n_iter
        self.A = A
        self.b = b
        self.D = np.matrix(np.diag(np.diag(self.A)))
        if iteration == 'richard_large':
            eigenvalues = np.sort(np.linalg.eigvals(self.A))
            mu_max = eigenvalues[-1]
            mu_min = eigenvalues[0]
            self.B = 0.9 / mu_max
        elif iteration == 'richard_small':
            eigenvalues = np.sort(np.linalg.eigvals(self.A))
            mu_max = eigenvalues[-1]
            mu_min = eigenvalues[0]
            self.B = 2 / (mu_max + mu_min)
        elif iteration == 'jacobi':
            self.B = self.D ** (-1)
        elif iteration == 'gauss_seidel':
            L = np.matrix(np.tril(self.A))
            self.B = (self.D + L) ** (-1)
        elif iteration == 'relaxed_jacobi':
            self.B = 2 / 3 * self.D ** (-1)
        else:
            self.B = 2 / 3 * self.D ** (-1)
            print 'iteration command not understood, using default \
            (relaxed jacobi) iteration procedure'

    def __mul__(self, var):
        """
            smoothing function call
        """
        if var == 'random':
            u_prev = np.matrix(np.random.random(len(self.b))).transpose()
        elif type(var) == int:
            u_prev = var * np.matrix(np.ones(len(self.b))).transpose()
        else:
            u_prev = var
        n = 0
        while n < self.n_iter:
            n += 1
            u = u_prev - self.B * (self.A * u_prev - self.b)
            u_prev = u
        u[0] = u[-1] = 0
        return u


class Multigrid(object):
    """
        creates a multigrid instance with a solve method
    """
    def __init__(self, N, x0, x1, n_iter, max_grid=16, iteration=False):
        """
            initializes grid atributes and number of iterations
        """
        self.N = N
        self.n = n_iter
        self.max_grid = max_grid
        self.x0 = x0
        self.x1 = x1
        self.iteration = iteration

    def cycle(self, N, b=0, U='random'):
        """
            solves a problem defined by the function initializer class,
            must be called initially with N = self.N
        """
        F = Func_init(N, self.x0, self.x1)
        A = F.A()
        if type(b) == int:
            b = F.b()
        S = Smoother(self.n, A, b, self.iteration)
        R = Res_op()
        I = Int_op()
        # pre smoothing with initially random guess
        u = S * U
        plt.plot(F.x(N), u, '--', label='{}'.format(N))
        # calculate residual
        b = R * (b - A * u)
        if N == self.max_grid * 2:
            F.N = N / 2
            A_H = F.A()
            #solve the error
            self.u = 4 * A_H ** (-1) * b  # FIXME magic fix (4*A_H)
            plt.plot(F.x(N / 2), self.u, label='{}'.format(N))
        else:
            self.cycle(N / 2, b, U=0)
        # prolong and update
        self.u = (I * self.u) + u
        # post smoothing
        self.u = S * self.u
        plt.plot(F.x(N), self.u, '.', label='{}'.format(N))


def test_run(plot=False):
    """
        solves the poisson problem with a full v cycle multigrid,
        note that N and max_grid must be of size 2**k, where k is a
        natural number.
    """
    M = Multigrid(N=128 * 4, x0=0, x1=1, n_iter=10, max_grid=128,
                  iteration='gauss_seidel')
    F = Func_init(M.N, M.x0, M.x1)
    M.cycle(M.N)
    F.exact()
    if plot:
        plt.plot(F.x(F.N), F.u_ex, '--', label='exact')
        plt.plot(F.x(M.N), M.u, label='iteration')
        plt.legend()
#        plt.show()


def test_smoother(n_iter, iteration=False, show=False):
    """
        A test of the smoothing operator. If one wants several plots;
        set show to false
    """
    N = 128
    x0 = 0
    x1 = 1
    F1 = Func_init(N, x0, x1)
    F2 = Func_init(N / 2, x0, x1)
    S1 = Smoother(n_iter, F1.A(), F1.b(), iteration)
    S2 = Smoother(n_iter, F2.A(), F2.b(), iteration)
    R = Res_op()
    I = Int_op()
    F1.exact()
    F2.exact()
    u_1 = S1 * 0
    u_2 = S2 * (R * u_1)
    u_solv = u_1 - I * (F2.A() ** (-1) * (R * (F1.A() * u_1 - F1.b())))
    error_s = np.mean(np.abs(F1.u_ex - u_solv))
    error_H = np.mean(np.abs(F2.u_ex - u_2))
    error_h = np.mean(np.abs(F1.u_ex - u_1))
    relative_sm_error = error_H / error_h
    relative_so_error = error_s / error_h
    print 'relative solving error = ', relative_so_error
    print 'relative smoothing error = ', relative_sm_error
    plt.plot(F1.x(), u_1, label='smoothing on h')
#    plt.plot(F2.x(), u_2, label='smoothing  on H')
    plt.plot(F1.x(), F1.u_ex, label='exact')
    plt.plot(F1.x(), u_solv, label='solve')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend()
    if show:
        plt.show()


def test_operators(N):
    """
        tests the restriction- and interpolation operator
    """
    uv = np.ones(N)
    uv[1], uv[-1] = 0, 0
    R = Res_op()
    I = Int_op()
    ures = R * uv
    uint = I * ures
    print uint


if __name__ == '__main__':
    np.random.seed(1234)
    test_run(True)
#   test_operators(128)
#    test_smoother(100)
    #run the smoothing operator with several number of iterations and
    # -plot them in one plot
#    for n in range(1,101):
#        if n == 1 or n == 10 or n == 100:
#            test_smoother(n, show = False)
    plt.legend()
    plt.show()
