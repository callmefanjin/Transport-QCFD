import numpy as np
import evolution

def build_p_grid(Lp,Np):
    deltap = 2 * Lp / Np
    p = np.linspace(-Lp, Lp, Np, endpoint=False)
    return p, deltap


def hermitian_split(Gamma):
    Gamma_dagger = Gamma.conj().T
    H1 = (Gamma+Gamma_dagger) / 2
    H2 = (Gamma-Gamma_dagger) / (2j)
    return H1,H2

def WPT_method(Y0,Lp,Np):
    p,deltap = build_p_grid(Lp,Np)
    W0 = np.outer(Y0, np.exp(-np.abs(p)))
    W0_shifted = np.fft.ifftshift(W0, axes=1)
    # Fourier tranfromation along the p axis
    psi0 = np.fft.fft(W0_shifted, axis=1, norm="ortho")
    eta = 2 * np.pi * np.fft.fftfreq(Np, d=deltap)
    return psi0,eta


def inverse_fft(psi_t):
    """Recover W(t, p) on the original p grid, retaining complex amplitudes."""
    Wt_shifted = np.fft.ifft(psi_t, axis=1, norm="ortho")
    return np.fft.fftshift(Wt_shifted, axes=1)



    
