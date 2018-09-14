# README

[![N|Python](https://www.python.org/static/community_logos/python-powered-w-100x40.png)](https://www.python.org)

- [Design goals](#design-goals)
- [Features](#features)
- [Installation](#installation)
- [My first simulation](#my-first-simulation)
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
- XDMF meshing
    - From CAD files
    - From Trelis/Cubit projects

# Installation
- 3-MOT
Install 3-MOT by cloning this git repository and install locally
```sh
git clone https://github.com/RemiTheWarrior/3-MOT
```
3-MOT relies on several softwares:

- FEniCS
To use 3-MOT you will need to use FEniCS in Docker. The instructions can be found [here](https://fenicsproject.org/download/).
Once the docker container pulled, you should be able to run FEniCS with the command ```fenicsproject run```


- Trelis
To be able to mesh your own geometry, you will need to the software named Trelis. The downloading instructions can be found [here](https://www.csimsoft.com/trelis). If you do not wish to install Trelis, you will be limited to default meshes that can be found in the Examples directory.

# My first simulation

Your first simulation will be a simple diffusion case on a unit cube.

In order to run your first simulation, all the parameters will be stored in ```MOT_parameters_cube_diffusion.json```. JSON files are used in order to store data in a convenient way. 

## Step 1: Specify the mesh

The mesh file is stored in the key named ```"structure_and_materials":{"project_files":}```. The file we will used is called ```"mesh_cube.xml"```. In this file 3-MOT will find the mesh and the markers for the surfaces and the volumes. This will be needed in order to assign different materials or initial values to volumes, or to set boundary conditions on the surfaces.
## Step 2: Define the materials

In this tutorial, the cube is made of only one material say Steel. We then set in the key  ```"structure_and_materials":{"materials":}``` as :
```js
"materials":[
    {
     "volumes":[1],
     "material":"polymer"
    }
]
```
In the key ```"materials"``` each element of the list is a directory which corresponds to a type of material. This directory contains the key ```"material":``` which contains the material as a string (according to the material library ```materials_properties.py```) and the key ```"volumes":``` where are stored all the volumes made of this material as a list. In our case there is only volume we have in our mesh is Volume ```1``` and is made of ```polymer```.

## Step 3: Define the physics to be solved
In this tutorial, we want to model the diffusion of hydrogen and the diffusion of heat. The physics parameters are stored in the key ```"physics":```.
### Step 3.1 : Set the solvers
  In order to tell 3-MOT to solve heat transfer and tritium diffusion we set the keys ```"solve_heat_transfer":1``` and ```"solve_tritium_diffusion":1```. As 3-MOT is also able to simulate the radioactive decay, we set the key ```"solve_with_decay":1```. In order to make the tritium diffusion temperature dependent, set the key ```"couple_tritium_diffusion_heat_transfer":1```.All the other keys are set to ```0```. 
```js
    "solve_with_decay":1,
    "solve_tritium_diffusion":1,
    "solve_laminar_flow":0,
    "solve_heat_transfer":0,
    "couple_tritium_diffusion_heat_transfer":1,
    "couple_tritium_diffusion_laminar_flow":0,
    "couple_heat_transfer_laminar_flow":0,
```
We now need to specify boundary conditions, initial values and source terms for the two physics.

### Step 3.2 : Define the heat diffusion problem

The parameters for the heat diffusion simulation are stored in the key ```"heat_transfers":```. This directory contains 4 keys : ```"boundary_conditions"```, ```"update_properties"```, ```"initial_value"``` and ```"source_terms"```.

- In the key ```"boundary_conditions"```, we can set 3 types of boundary conditions : Dirichlet, Neumann and Robin (convective condition).
We set the boundary conditions as follow:
```js
"boundary_conditions":{
        "neumann":[
          {
            "surface":[3],
            "value":"-1"
          }
            ],
        "dc":[
          {
            "surface":[1,2],
            "value":"300"
          }
        ],
        "robin":[
          {
            "surface":[4],
            "t_amb":"320",
            "h_coeff":"4000"
          }
        ]
      }
```
Surfaces ```5```and ```6``` are by default insulated. 
- We set the key ```"update_properties":0``` to specify that we don't want the thermal properties to be updated at each time step in order to save computing time.
- The initial value is set to ```"initial_value":"303+3*sin(x[0])"```
- The source term is set to zero : ```"source_terms":[]```


### Step 3.3 : Define the tritium diffusion problem
As in Step 3.2, the parameters for the tritium diffusion simulation are stored in the key ```"tritium_diffusion":```. This directory contains 4 keys : ```"boundary_conditions"```, ```"update_properties"```, ```"initial_value"``` and ```"source_terms"```.

- In the key ```"boundary_conditions"```, we can set 3 types of boundary conditions : Dirichlet, Neumann and Robin (convective condition).
We set the boundary conditions as follow:
```js
"boundary_conditions":{
        "neumann":[
            ],
        "dc":[
          {
            "surface":[3,4],
            "value":"0"
          }
        ],
        "robin":[
          {
            "surface":[6],
            "value":"1e-10*c_n"
          }
        ]
      }
```
Surfaces ```1```and ```2``` are by default insulated. 

- We set the key ```"update_properties":1``` to specify that we want the diffusion coefficient to be temperature dependent updated at each time step.
- The initial value is set to ```"initial_value":"0"```
- The source term is set to zero : ```"source_terms":[]```
## Step 4 : Define the solving parameters

The next step is now to define the solving parameters in the key ```"solving_parameters"```. These are the type of study and, if the study is time-dependent, the number of time steps and the final time of the simulation. Our case is transient, the parameters are :
```js
"solving_parameters":{
    "study":"transient",
    "final_time":30000, 
    "number_of_time_steps":100
  }
```

If we wanted to simulate the steady state of the problem, the parameters would be:
```js
"solving_parameters":{
    "study":"steady_state"
  }
```

## Step 5 : Post computing
Our problem is now fully defined. The last step is to tell the solver where to store the solution and eventually define the outputs we need. 

### Step 5.1 : Define the main solution file

The solution file path is stored in the key ```"output_file":```. In our case it will be :
```js
"output_file":"Examples/Cube_tuto/Solution/solution.pvd"
```

Both the temperature and tritium concentration fields will be stored in this file.

### Step 5.2 : Define the post processing computations (optionnal)
Post computing is facilitate in 3-MOT. Surface fluxes, averages, minimum values and maximum values can be computed for each solution. You can even compute your own custom expression. 

```js
"post_processing":{
    "heat_transfers":{
      "surface_flux":[2],
      "volume_average":[],
      "volume_minimum":[1],
      "volume_maximum":[1],
      "custom":[],
      "output_file":"Examples/Cube_tuto/post-processingHT.csv"
    },
    "tritium_diffusion":{
      "surface_flux":[1,2,3],
      "volume_average":[2],
      "volume_minimum":[],
      "volume_maximum":[],
      "custom":["assemble(solution*dx(1))"],
      "output_file":"Examples/Cube_tuto/post-processing-tritium-diffusion_quality5.csv"
      }
```

# ToDo
- Meshing grom GMSH
- Write more tests
- Write more tutorials