"""Compiled evaluation of exactly the existing ordered Pauli product."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from qiskit.quantum_info import Operator, SparsePauliOp

ROOT = Path(__file__).resolve().parent

def prepare_input(H, state0, t, r, src):
    terms = SparsePauliOp.from_operator(Operator(H), rtol=1e-10)
    coeffs = np.real(terms.coeffs)
    assert len(coeffs) and not np.any(terms.paulis.phase)
    n = terms.num_qubits
    assert len(state0) == 1 << n and r >= 1
    weights = 1 << np.arange(n, dtype=np.int64)
    xm = terms.paulis.x.astype(np.int64) @ weights
    zm = terms.paulis.z.astype(np.int64) @ weights
    ny = np.array([int(x & z).bit_count() for x, z in zip(xm, zm)])
    angles = coeffs * (t/r) / 2
    angles[-1] *= 2
    cosines = np.cos(angles)
    factors = -1j*np.sin(angles)*(-1j)**ny
    with Path(src).open('wb') as f:
        np.array([len(state0), len(coeffs), r], dtype='<i4').tofile(f)
        for arr in (xm, zm): np.asarray(arr, dtype='<i4').tofile(f)
        for arr in (cosines, factors.real, factors.imag, state0.real, state0.imag):
            np.asarray(arr, dtype='<f8').tofile(f)
    return r*(2*len(coeffs)-1)

def evolve_pauli_product(H, state0, t, r):
    with tempfile.TemporaryDirectory(prefix='pauli_', dir=ROOT) as folder:
        src = Path(folder)/'input.bin'
        dst = Path(folder)/'output.bin'
        gates = prepare_input(H, state0, t, r, src)
        subprocess.run([str(ROOT/'PauliProductWorker.exe'), str(src), str(dst)], check=True)
        data = np.fromfile(dst, dtype='<f8')
    dim = len(state0)
    assert len(data) == 2*dim
    return data[:dim] + 1j*data[dim:], gates
