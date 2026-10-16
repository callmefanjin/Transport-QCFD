import numpy as np

def build_Lk(A,K):
    if K == 1:
        return A.copy()
    N = A.shape[0]
    Lk = np.zeros((N**K,N**K))
    for r in range(K):
        I_left = np.eye(N**r)
        I_right = np.eye(N**(K-r-1))
        mid_term = np.kron(np.kron(I_left, A),I_right)
        Lk = Lk+mid_term
    return Lk

def build_Bk(F,K):
    if K == 1:
        return F.copy()
    N = F.shape[0]
    Bk = np.zeros((N**K,N**(K+1)))
    for r in range(K):
        I_left = np.eye(N**r)
        I_right = np.eye(N**(K-r-1))
        mid_term = np.kron(np.kron(I_left, F),I_right)
        Bk = Bk+mid_term
    return Bk

    
def build_Gamma(A,F,K):
    N = A.shape[0]
    sizes = [N**k for k in range(1, K + 1)] 
    offsets = np.concatenate(([0], np.cumsum(sizes)))
    M = offsets[-1]
    Gamma = np.zeros((M,M))
    for k in range(1, K+1):
        start = offsets[k-1]
        end = offsets[k]
        Gamma[start:end,start:end] = build_Lk(A,k)
        if k<K:
            next_end = offsets[k+1]
            Gamma[start:end, end:next_end] = build_Bk(F,k)
    return Gamma

def build_Y(phi,K):
    N = len(phi)
    sizes = [N**k for k in range(1, K + 1)]  
    offsets = np.concatenate(([0], np.cumsum(sizes))) 

    Y = np.zeros(offsets[-1])
    Phi_last = phi.copy()
    for k in range(1,K+1):
        start = offsets[k-1]
        end = offsets[k]
        Y[start:end] = Phi_last
        if k < K:
            Phi_last = np.kron(Phi_last, phi)
    return Y
