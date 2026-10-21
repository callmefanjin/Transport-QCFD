# Transport-QCFD

Transport-QCFD is a Python computational code for implementing quantum algorithms for transport equations in computational fluid dynamics. The current implementation solves a nonlinear one-dimensional advection–diffusion–reaction equation using Carleman linearisation, Schrödingerisation and Suzuki–Trotter evolution.

## Nonlinear transport

![Quantum-algorithm and finite-difference solutions on two spatial grids](nonlinear-transport-QCFD/figures/selected_grid_results.png)

*Recovered field (`phit`) compared with a nonlinear finite-difference reference. The displayed results use classical statevector simulation.*

## Core capabilities

- Carleman lifting of quadratic nonlinear dynamics.
- Warped-phase transformation and Fourier-mode Hamiltonian evolution.
- Qiskit Pauli circuits and equivalent direct statevector evaluation.
- Field recovery and comparison with a Radau reference solution.

[Model and results](nonlinear-transport-QCFD/README.md) · [Run the code](nonlinear-transport-QCFD/RUNNING.md)
