# Full solver results on two spatial grids

We solve $\phi_t+u\phi_x=D\phi_{xx}+\alpha\phi+\kappa\phi^2$ with the periodic boundary condition and Gaussian initial condition shown in the figure. All case parameters are marked on the plots.

**Finite difference** uses second-order central spatial differences and Radau time integration of the full nonlinear system ($\mathrm{rtol}=10^{-10}$, $\mathrm{atol}=10^{-12}$). **phit** uses full-tensor Carleman truncation ($K=2$), WPT and FFT, followed by a second-order symmetric Suzuki--Trotter product over Pauli terms with $r$ steps of size $T/r$. An inverse FFT recovers the field at $p=0$, and the first $N_x$ components give phit. The Pauli rotations are evaluated by equivalent classical statevector operations; these are simulator results, not quantum-hardware measurements.

![Finite-difference and phit solutions on two spatial grids, with all case parameters](./figures/selected_grid_results.png)

The curves show real parts; $E=\|\mathrm{phit}-\phi_{\mathrm{FD}}\|_2/\|\phi_{\mathrm{FD}}\|_2$ uses the full complex solution.

The results show that the solution recovered from the truncated Carleman linear system, through WPT and Trotter evolution, can approximate the spatial distribution of the corresponding nonlinear finite-difference solution.

The parameter tests also show that spatial refinement requires a corresponding increase in the Trotter step count $r$ to control the evolution error. The discrete convection and diffusion coefficients scale as $1/\Delta x$ and $1/\Delta x^2$, respectively, so keeping $T/r$ fixed can give a less accurate product-formula approximation on a finer grid. Spatial refinement should therefore be accompanied by a check of Trotter convergence.
