![OrbitalTrajectories.jl](assets/logo.svg)

**OrbitalTrajectories.jl** is a modern orbital trajectory design library for Julia,
providing tools and methods for the design, optimisation, and analysis of astrodynamical
models and spacecraft trajectories.

!!! warning "In Development"
    OrbitalTrajectories.jl and its documentation is currently in development
    and is provided as a pre-1.0 release. Please refer to the examples/tutorials
    provided, and feel free to provide feedback to the package maintainers.

## Getting Started

Install the package with ```add OrbitalTrajectories``` in Julia's package
mode (run the Julia console and press ```]```, or alternatively ```using Pkg;
Pkg.add("OrbitalTrajectories")```).

To later update to the newest release, simply do ```update OrbitalTrajectories```.

## Features

- High-performance and highly accurate trajectory propagation in several astrodynamical models:
  - Composable models built directly from symbolic mathematical equations of motion
  - Compiled performance equivalent or exceeding that of MATLAB/Fortran competitors
  - High-order interpolation of resulting trajectories
  - Dynamic conversions between inertial and rotating reference frames
- Implemented astrodynamical models:
  - **Circular Restricted 3-Body Problem** (`ER3BP`)
  - **Elliptic Restricted 3-Body Problem** (`ER3BP`)
  - **Bi-Circular Restricted 4-Body Problem** (`BC4BP`)
  - **Restricted N-Body Problem** (`EphemerisNBP`)
    - Body states computed from NASA NAIF (SPICE) ephemerides downloaded automatically as needed
    - Automatic SPICE dynamic reference frames built for conversions
- Computation of State Transition Matrices (STMs):
  - Using fast, accurate Automatic Differentiation (AD) for every model
  - Using symbolically computed Variational Equations (VE) for some models
  - Hand-tuned variational equations for CR3BP model
- Generic single-shooting differential corrector, supporting:
  - Finding axisymmetric (quasi-)periodic orbits in every model
- Simple parameter continuation for computing families of orbits
- Plotting recipes
  - Plot trajectories in inertial/rotating frames
  - Show primary body, secondary body, and libration points
  - Plot orbit families with custom-valued colourbars
- Support for Julia's highly-composable state-of-the-art ecosystem for machine learning, statistical analysis, and more

## License & References
Distributed under the [Apache License 2.0](LICENSE)

Copyright 2021 Dan Padilha ([dpadilha.com](http://www.dpadilha.com))

If you use OrbitalTrajectories.jl in a scientific project that leads to a publication, we'd appreciate you citing our paper as follows:
```
@inproceedings{OrbitalTrajectories,
  url = {https://www.researchgate.net/publication/348929030_Modern_Numerical_Programming_with_Julia_for_Astrodynamic_Trajectory_Design},
  year = {2021},
  publisher = {AAS/AIAA},
  author = {Padilha, Dan and Dei Tos, Diogene Alessandro and Baresi, Nicola and Kawaguchi, Junichiro},
  title = {Modern Numerical Programming with Julia for Astrodynamic Trajectory Design},
  booktitle = {31st AAS/AIAA Space Flight Mechanics Meeting}
}
```