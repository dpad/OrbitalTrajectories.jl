![OrbitalTrajectories.jl](docs/src/assets/logo.svg?raw=true)

| **Documentation**   |  **Tests**     |  Citation| License
|:--------:|:----------------------:|:-----:|:-----:|
|[![](https://img.shields.io/badge/docs-online-blue.svg)](https://dpad.github.io/OrbitalTrajectories.jl/stable/)| [![CI](https://github.com/dpad/OrbitalTrajectories.jl/workflows/CI/badge.svg)](https://github.com/dpad/OrbitalTrajectories.jl/actions) [![codecov](https://codecov.io/gh/dpad/OrbitalTrajectories.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dpad/OrbitalTrajectories.jl) | [![reference](https://img.shields.io/badge/Thesis%20%28as%20of%20v0.2.0dev%29-Padilha%202021-brightgreen)](https://www.researchgate.net/publication/353906343_Composable_Astrodynamics_Software_via_Multiple_Dispatch_and_Meta-Programming) [![reference](https://img.shields.io/badge/Paper%20%28as%20of%20v0.1.1%29-Padilha%20et%20al%202021-brightgreen)](https://dpadilha.com/Padilha%20-%20AAS%2021-303.pdf) | [![GitHub license](https://img.shields.io/github/license/dpad/OrbitalTrajectories.jl)](LICENSE)

*OrbitalTrajectories.jl is a modern orbital trajectory design, optimisation, and analysis library for Julia, providing methods and tools for designing spacecraft orbits and transfers via high-performance simulations of astrodynamical models.*

[![video](https://img.shields.io/badge/Presentation-AAS%2FAIAA%20Conference-brightgreen)](https://www.youtube.com/watch?v=FMVOUvWNlLE) [![video](https://img.shields.io/badge/Presentation-JuliaCon%202021-brightgreen)](https://www.youtube.com/watch?v=iJr_lU7_7Go)

---

## Getting Started

Install the package with ```add OrbitalTrajectories``` in Julia's package
mode (run the Julia console and press ```]```, or alternatively ```using Pkg;
Pkg.add("OrbitalTrajectories")```).

To later update to the newest release, simply do ```update OrbitalTrajectories```.

---

## License & References
Distributed under the [Apache License 2.0](LICENSE)

Copyright 2021 Dan Padilha ([dpadilha.com](http://www.dpadilha.com))

If you use OrbitalTrajectories.jl in a scientific project that leads to a publication, we'd appreciate you citing our work as follows:
```
@mastersthesis{Padilha2021,
  doi = {10.13140/RG.2.2.14175.79525},
  url = {http://rgdoi.net/10.13140/RG.2.2.14175.79525},
  author = {Padilha, Dan},
  language = {en},
  title = {Composable Astrodynamics Software via Multiple Dispatch and Meta-Programming},
  school = {The University of Tokyo},
  year = {2021},
  month = aug
}

@inproceedings{OrbitalTrajectories,
  url = {https://www.researchgate.net/publication/348929030_Modern_Numerical_Programming_with_Julia_for_Astrodynamic_Trajectory_Design},
  year = {2021},
  publisher = {AAS/AIAA},
  author = {Padilha, Dan and Dei Tos, Diogene Alessandro and Baresi, Nicola and Kawaguchi, Junichiro},
  title = {Modern Numerical Programming with Julia for Astrodynamic Trajectory Design},
  booktitle = {31st AAS/AIAA Space Flight Mechanics Meeting}
}
```
