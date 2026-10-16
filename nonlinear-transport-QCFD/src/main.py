import schrodingerisation
import transport
import evolution
import carleman
import numpy as np
from qiskit.quantum_info import Statevector

def main():
    # 1. set parameters
    start = 0
    end = 2
    deltax = 0.1

    u = 10 # convection speed
    D = 0.5 # diffusion coefficient
    alpha = -1 # linear source term
    kappa = 0.1 # nonlinear quadratic source term

    T = 1 # time

    K = 2 # carleman order

    r = 8 # trotter layer

    Lp = 16
    Np = 16

    x = transport.build_grid(start, end, deltax)

    phi0 = transport.initial_condition(x)  # physical initial value phi
    A = transport.linear_matrix(x,u,D,alpha)
    F = transport.nonlinear_part(x,kappa)

    Y0 = carleman.build_Y(phi0,K) # extended carleman system initial value Y0
    Gamma = carleman.build_Gamma(A,F,K)

    p, deltap = schrodingerisation.build_p_grid(Lp,Np)
    psi0, eta = schrodingerisation.WPT_method(Y0,Lp,Np)  # need to solev len(eta) equations
    H1,H2 = schrodingerisation.hermitian_split(Gamma)

    # correct the size of carleman system to conform binary
    M = len(psi0)
    dim = 1 << (M - 1).bit_length()  

    psi_t = np.zeros_like(psi0)
    for i, eta_i in enumerate(eta):
        HT = eta_i*H1 - H2
        HT_pad = np.zeros((dim, dim), dtype=complex)
        HT_pad[:M, :M] = HT

        psi0_i = psi0[:,i]
        amplitude = np.linalg.norm(psi0_i)
        if amplitude == 0:
            continue
        
        state0 = np.zeros(dim, dtype=complex)
        state0[:M] = psi0_i / amplitude

    
        qc = evolution.trotter_2nd_evolution(HT_pad,T,r)
        state_t = Statevector(state0).evolve(qc).data
        psi_t[:, i] = amplitude * state_t[:M]

    Wt = schrodingerisation.inverse_fft(psi_t)
    # choose p = 0
    j0 = np.argmin(np.abs(p))   
    if not np.isclose(p[j0], 0.0):  
        raise ValueError("p=0 not exsit in p grid")
    Yt = Wt[:,j0]
    phit = Yt[:len(phi0)]

    Yt_reference = evolution.exact_evolution(Gamma, Y0, T)
    phi_reference = Yt_reference[:len(x)]
    relative_error = (np.linalg.norm(phit - phi_reference)/ np.linalg.norm(phi_reference))
    return x, phit, phi_reference, relative_error




if __name__ == "__main__":
    main()

