"""Evaluate the original symmetric Pauli product without constructing gate objects.

This is classical statevector evaluation, not a hardware circuit execution.
The decomposition tolerance, term order, and angles match evolution.py.
"""
import numpy as np
from qiskit.quantum_info import Operator, SparsePauliOp


def evolve_pauli_product(H, state0, t, r):
    terms = SparsePauliOp.from_operator(Operator(H), rtol=1e-10)
    coeffs = np.real(terms.coeffs)
    if len(coeffs) == 0:
        raise ValueError("Empty Pauli decomposition")
    if np.any(terms.paulis.phase):
        raise ValueError("Expected canonical Pauli labels without global phases")
    dim = len(state0)
    n = terms.num_qubits
    if dim != 1 << n or r < 1:
        raise ValueError("Invalid state dimension or Trotter repetitions")
    weights = 1 << np.arange(n, dtype=np.int64)
    xm = terms.paulis.x.astype(np.int64) @ weights
    zm = terms.paulis.z.astype(np.int64) @ weights
    ny = np.array([int(x & z).bit_count() for x, z in zip(xm, zm)])
    indices = np.arange(dim, dtype=np.int32)
    permutations = indices[:, None] ^ indices[None, :]
    parity = np.array([i.bit_count() % 2 for i in range(dim)], dtype=np.int8)
    signs = 1 - 2 * parity[indices[:, None] & indices[None, :]]
    angles = coeffs * (t / r) / 2
    angles[-1] *= 2
    cosines = np.cos(angles)
    # (P v)[i] = (-i)^number_of_Y * (-1)^popcount(z & i) * v[i xor x].
    factors = -1j * np.sin(angles) * (-1j)**ny
    state = np.array(state0, dtype=complex, copy=True)
    order = list(range(len(coeffs))) + list(range(len(coeffs)-2, -1, -1))
    for _ in range(r):
        for k in order:
            state = (cosines[k] * state
                     + factors[k] * signs[zm[k]] * state[permutations[xm[k]]])
    return state, r * (2 * len(coeffs) - 1)
