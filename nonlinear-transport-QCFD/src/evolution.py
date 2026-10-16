import numpy as np
from scipy.linalg import expm
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator, SparsePauliOp
from qiskit.circuit.library import PauliEvolutionGate

def exact_evolution(Gamma,Y0,t):
    return expm(t*Gamma) @ Y0

def trotter_2nd_evolution(A, t, r):
    pauli_A = SparsePauliOp.from_operator(Operator(A),rtol=1e-10)
    coeffs = np.real(pauli_A.coeffs) # A must be Hermitian,so coeffs should be real
    labels = pauli_A.paulis.to_labels()
    if len(coeffs) == 0:
        raise ValueError("Pauli decomposition find the coeffs is empty,please check the tolerance of the Pauli decomposition")
    n_qubits = int(np.log2(A.shape[0]))
    qc = QuantumCircuit(n_qubits)
    m = len(labels)
    tau = t / r
    for r_i in range(r):
        for k in range(m-1):
            label_k = labels[k]
            theta_k = coeffs[k]*tau/2
            single = SparsePauliOp([label_k], [theta_k])
            qc.append(PauliEvolutionGate(single, time=1.0), range(n_qubits))
        theta = coeffs[m - 1]*tau
        single = SparsePauliOp([labels[m - 1]], [theta])
        qc.append(PauliEvolutionGate(single, time=1.0), range(n_qubits))
        for k in reversed(range(m - 1)):
            label_k = labels[k]
            theta_k = coeffs[k] * tau / 2
            single = SparsePauliOp([label_k], [theta_k])
            qc.append(PauliEvolutionGate(single, time=1.0), range(n_qubits))
        # remember it is e^(-iAt),there is a minus sign in the exponent
    return qc


