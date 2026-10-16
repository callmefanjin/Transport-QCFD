# Running the examples

Run the following commands from this folder. Python 3.12 was used for the recorded results.

```sh
python -m pip install -r requirements.txt
```

## Files

- `README.md`: the extracted results section, with a relative image link.
- `src/`: the five original solver modules, copied unchanged.
- `examples/run_case.py`: a portable copy of the experiment driver; parameter overrides do not edit `src/`.
- `examples/pauli_action.py`: direct evaluation of the original ordered Pauli product.
- `examples/compiled_pauli_action.py` and `PauliProductWorker.cs`: optional equivalent compiled evaluation.
- `figures/`: the combined parameter-labelled PNG and vector PDF.
- `outputs/`: the two selected cases, including parameters, recorded summaries, CSV curves and compact NPZ plotting data.
- `validation/`: previously recorded comparisons with the Qiskit circuit; these are not new hardware runs.

`src/main.py` retains the original development defaults. Use the commands below to reproduce the parameters shown in the figure.

## Regenerate the combined figure from saved data

```sh
python examples/plot_results.py
```

This does not rerun the solver. Times New Roman must be installed to reproduce the figure style; the plotting helper reports an error if it is unavailable.

## Recompute the selected cases

```sh
python examples/run_case.py --nx 30 --start -15 --end 15 --mu -10 --np 64 --r 8 --time 1 --pauli-action --name rerun_dx1_np64_r8
python examples/run_case.py --nx 60 --start -15 --end 15 --mu -10 --np 64 --r 32 --time 1 --pauli-action --name rerun_dx0p5_np64_r32
```

The remaining parameters are `u=10`, `D=0.5`, `alpha=-1`, `kappa=0.1`, `K=2`, and `Lp=16`. The nonlinear reference uses the same spatial grid and Radau time integration. Results are saved under new `outputs/` subfolders; existing cases are not overwritten. These dense statevector calculations can take substantial time.

For the optional compiled backend on Windows, build the included C# source:

```powershell
python examples/build_worker.py
```

Then replace `--pauli-action` in either command with `--compiled-pauli`. The Pauli terms, ordering, rotation angles and repetition count are unchanged. The packaged driver evaluates all frequencies sequentially; the recorded fine-grid run used three workers and conjugate-frequency reuse for the same mathematical product.

The bundled NPZ files contain the curve data needed for plotting, including complex `phit`; the large intermediate Fourier states and temporary checkpoints are omitted. The error is measured against the nonlinear finite-difference solution on the same grid. No spatial-convergence claim is made because the two selected cases also use different `r` values.
