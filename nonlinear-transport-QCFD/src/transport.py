import numpy as np


def build_grid(start, end, deltax):
    return np.arange(start, end, deltax)

def initial_condition(x):
    return np.exp(-(x - 1)**2 / 0.1)

def convection_part(x, u):
    N = len(x)
    deltax = x[1] - x[0]
    convection_matrix = np.zeros((N, N))
    for i in range(N):
        if i <  N-1:
            convection_matrix[i, i+1] = -u/(2*deltax)
        if i > 0:
            convection_matrix[i, i-1] = u/(2*deltax)
    return convection_matrix

def diffusion_part(x, D):
    N = len(x)
    deltax = x[1] - x[0]
    diffusion_matrix = np.zeros((N, N))
    for i in range(N):
        if i > 0:
            diffusion_matrix[i, i-1] = D/(deltax**2)
        if i < N-1:
            diffusion_matrix[i, i+1] = D/(deltax**2)
        diffusion_matrix[i, i] = -2*D/(deltax**2)
    return diffusion_matrix

def alpha_part(x, alpha):
    N = len(x)
    alpha_matrix = np.zeros((N, N))
    for i in range(N):
        alpha_matrix[i, i] = alpha
    return alpha_matrix

def linear_matrix(x,u,D,alpha):
    N = len(x)
    deltax = x[1] - x[0]
    A = convection_part(x, u) + diffusion_part(x, D) + alpha_part(x, alpha)
    # considering the periodic boundary condition
    A[0, N-1] = D/(deltax**2)+u/(2*deltax)
    A[N-1, 0] = D/(deltax**2)-u/(2*deltax)
    return A

def linear_rhs(x,u,D,alpha,phi):
    return linear_matrix(x,u,D,alpha) @ phi

def nonlinear_part(x,kappa):
    N = len(x) 
    NN = N**2 # the size number comes form the tensor product
    F = np.zeros((N, NN))
    for i in range(N):
        F[i, i*N+i] = kappa
    return F

def nonlinear_rhs(x,u,D,alpha,kappa,phi):
    phiphi = np.kron(phi, phi)
    return nonlinear_part(x,kappa) @ phiphi

def rhs(x,u,D,alpha,kappa,phi):
    return linear_rhs(x,u,D,alpha,phi) + nonlinear_rhs(x,u,D,alpha,kappa,phi)




