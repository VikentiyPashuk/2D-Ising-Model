# 2D-Ising-Model

## Overview

Fortran code for simulating the 2D Ising model using the Monte Carlo method for integration, a cornerstone of statistical physics, enabling exploration of phase transitions and critical phenomena in materials. In the simulations, negative spins are represented as $-1$ and positive spins as $1$. The coupling constant  $J$ and the Boltzmann constant $k_B$ are both set to $1$. The simulations record magnetization, energy, heat capacity, and magnetic susceptibility for analysis. Visualizations of the results can be found in the `plots.ipynb` notebook.

## Usage

### Running Simulations

To run simulations with different lattice sizes and Monte Carlo steps (MCS), execute the `run_ising.sh` shell script:

```bash
./run_ising.sh
```

This script iterates over various lattice sizes and MCS values, running the `ising-2d` executable and storing results in separate files.

### Analyzing Results

After running simulations, analyze the obtained data using provided analysis tools:

- **Initial and Final Matrices**: Execute `matrix_ising` to generate initial and final spin matrices (`spins_initial_matrix.txt` and `spins_final_matrix.txt`), set the desirable temperature in the `matrix_ising.f90` code to see the different spin arrangements.
- **Critical Temperature**: Compute the critical temperature using `Tc_ising` with 5000 MCS. Results are stored in `tc_ising.dat`.

## Compilation

Compile the Fortran code using the provided Makefile:

```bash
make
```

This will compile the `ising-2d`, `Tc_ising`, and `matrix_ising` executables.

## Visualization

For visualization of magnetization, energy, heat capacity, magnetic susceptibility, intial/final matrices, and critical temperature refer to `plots.ipynb`.

## Cleaning Up

To remove compiled executables, run:

```bash
make clean
```

## Dependencies

- Fortran Compiler (e.g., `gfortran`)
- Jupyter Notebook (for visualization, plotting, and analysis)

## Note

Ensure you have the necessary permissions to execute shell scripts and compile Fortran code.
