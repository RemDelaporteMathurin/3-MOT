# README

[![N|Python](https://www.python.org/static/community_logos/python-powered-w-100x40.png)](https://www.python.org)

- [Design goals](#design-goals)
- [Features](#features)
- [Installation](#installation)
- [Todo](#todo)


# Design Goals
3-MOT is a solver using Finite Element Methods (FEM) designed to model multi-material domains based on CAD files. This software provides a user friendly interface to quickly implement complex boundary conditions and simulate multi-physics. 

# Features

- Bulk models (stationary and transient):
    - Tritium Diffusion
    - Heat transfers
    - Computational Fluid Dynamics
- Boundary conditions:
    - Dirichlet
    - Neumann
    - Robin
- Multi-material domains
- Post-processing computations:
    - Fluxes
    - Averages
    - Minimum/Maximum
    - Custom computations

# Installation
- 3-MOT
Install 3-MOT by cloning this git repository and install locally
```sh
git clone https://github.com/RemiTheWarrior/3-MOT
```
3-MOT relies on several softwares:

- FEniCS
To use 3-MOT you will need to install FEniCS first. The instructions can be found [here](https://fenicsproject.org/download/).


- Trelis
To be able to mesh your own geometry, you will need to the software named Trelis. The downloading instructions can be found [here](https://www.csimsoft.com/trelis). If you do not wish to install Trelis, you will be limited to default meshes.

# ToDo










