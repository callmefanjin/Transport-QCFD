# Run the Code

Run these commands from `nonlinear-transport-QCFD/`. The saved results used Python 3.12.

## Setup

```sh
python -m pip install -r requirements.txt
```

Plotting requires Times New Roman.

## Plot saved results

```sh
python examples/plot_results.py
```

## Compute the two cases

```sh
python examples/run_case.py --nx 30 --start -15 --end 15 --mu -10 --np 64 --r 8 --time 1 --pauli-action --name rerun_dx1_np64_r8
python examples/run_case.py --nx 60 --start -15 --end 15 --mu -10 --np 64 --r 32 --time 1 --pauli-action --name rerun_dx0p5_np64_r32
```

These commands use direct statevector evaluation of the ordered Pauli product. Omit `--pauli-action` to use Qiskit circuit simulation. Each run requires a new `--name` and saves its parameters, fields and comparison with the nonlinear reference under `outputs/`. Dense statevector calculations can be expensive.

[Model, parameters and results](README.md)
