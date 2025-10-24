# MATLAB Magnetohydrodynamics Simulation

This repository contains a MATLAB project for simulating magnetohydrodynamics (MHD) using a Chebyshev collocation method. It solves a system of partial differential equations to model the behavior of electrically conducting fluids in the presence of magnetic fields.

## How to Run

1.  **Open MATLAB:** Launch the MATLAB environment.
2.  **Set Parameters:** Open the `parameters.m` file and adjust the simulation parameters as needed.
3.  **Run the Simulation:** Execute the `run_parallel.m` script to start the simulation. This script will iterate through the defined parameters and save the results in the `runs` directory.

## File Descriptions

*   `parameters.m`: Configures the simulation's physical and numerical parameters, such as the Ekman number, Prandtl number, and magnetic Prandtl number.
*   `run_parallel.m`: The main script that drives the simulation. It sets up the computational domain, differentiation matrices, and iterates through the parameter space to find solutions.
*   `solver_new.m`: Solves the generalized eigenvalue problem using the `eigs` function to find the eigenvalues and eigenvectors of the system.
*   `build_matrix.m`: Constructs the A and B matrices for the generalized eigenvalue problem Ax = λBx, which represent the discretized differential equations.
*   `BC.m`: Applies the boundary conditions to the system. The file is not present in the repository, but it is called from `run_parallel.m`.
*   `Variable_coeffs_new.m`: Defines the variable coefficients for the differential equations based on the physical parameters.
*   `cheb.m`: Generates the Chebyshev differentiation matrices.
*   `custom_plot.m`: A script for creating custom plots of the simulation results. Not called by the main script, but can be used for post-processing.
*   `call_BC.m`: A helper script for managing boundary conditions.

## Output

The simulation saves its results as `.mat` files in a `runs` directory (which is created if it doesn't exist). Each file is named according to the `chi` and `els` parameters, for example, `chi1_els0.001.mat`. These files contain the computed eigenvalues, eigenvectors, and other relevant data for each simulation run.
