from fenics import *
from dolfin import *
import numpy as np
import csv
import sys
import os
import argparse
import json
import ast
from pprint import pprint
from materials_properties import *
import inspect




def get_apreprovars(apreprovars):
    return 'MOT_parameters_RCB.json'
def get_databases(name_database):
    with open(name_database) as f:
        data = json.load(f)
    data=byteify(data)
    return data
def get_solvers(data):
    print('Getting the solvers')
    if data['physics']['solve_heat_transfer']==1:
        solve_temperature=True
    else:
        solve_temperature=False
    if data['physics']['solve_tritium_diffusion']==1:
        solve_diffusion=True
    else:
        solve_diffusion=False
    if data['physics']['diffusion_coeff_temperature_dependent']==1:
        solve_diffusion_coefficient_temperature_dependent=True
    else:
        solve_diffusion_coefficient_temperature_dependent=False
    if data['physics']['solve_with_decay']==1:
        solve_with_decay=True
    else:
        solve_with_decay=False
    calculate_off_gassing=True
    return solve_temperature,solve_diffusion,solve_diffusion_coefficient_temperature_dependent,solve_with_decay

def get_solving_parameters(data):
    Time = float(data["solving_parameters"]['final_time'])  #60000*365.25*24*3600.0# final time 
    num_steps = data['solving_parameters']['number_of_time_steps'] # number of time steps
    dt = Time / num_steps # time step size
    return Time,num_steps,dt

def define_mesh(data):
    print('Defining mesh')
    # Read in Mesh and markers from file
    mesh = Mesh()
    xdmf_in = XDMFFile(mesh.mpi_comm(), str(data['mesh_file']))
    xdmf_in.read(mesh)
    # prepare output file for writing by writing the mesh to the file
    #xdmf_out = XDMFFile(str(data['mesh_file']).split('.')[1]+'_from_fenics.xdmf')

    subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
    print('Number of cell is '+ str(len(subdomains.array())))
    return mesh

def define_functionspaces(data):
    print('Defining Functionspaces')
    V = FunctionSpace(mesh, 'P', 1) #FunctionSpace of the solution c
    V0 = FunctionSpace(mesh, 'DG', 0) #FunctionSpace of the materials properties
    return V,V0

def define_initial_values(solve_temperature,solve_diffusion,data,V):
    print('Defining initial values')
    ##Tritium concentration
    if solve_diffusion==True:
      #print(str(data['physics']['tritium_diffusion']['initial_value']))
      iniC = Expression(str(data['physics']['tritium_diffusion']['initial_value']),degree=2) 
      c_n = interpolate(iniC, V)
    ##Temperature
    if solve_temperature==True:
      #print(str(data['physics']['heat_transfers']['initial_value']))
      iniT = Expression(str(data['physics']['heat_transfers']['initial_value']),degree=2) 
      T_n = interpolate(iniT, V)
    return c_n,T_n







if __name__=="__main__":
    apreprovars=get_apreprovars()
    data=get_databases(apreprovars) #This returns an object data=json.load()
    solve_temperature,solve_diffusion,solve_diffusion_coefficient_temperature_dependent,solve_with_decay=get_solvers(data) #Gets the solvers
    solving_parameters=get_solving_parameters(data) #Gets the parameters (final time, time steps...)
    mesh=define_mesh(data)
    V,V0=define_functionspaces(data)
    initial_values=define_initial_values(solve_temperature,solve_diffusion,solve_diffusion_coefficient_temperature_dependent,solve_with_decay,data,V,V0)
    BC_diffusion=define_BC_diffusion(data)
    BC_heat_transfer=define_BC_heat_transfer(data)
    volume_marker=get_volume_markers(mesh)
    properties=define_materials_properties(V,V0,data,volume_marker)
    F=define_variational_problem_diffusion(solve_diffusion,solve_diffusion_coefficient_temperature_dependent,solve_with_decay,V,V0,data,BC_diffusion)
    FT=define_variational_problem_heat_transfer(solve_temperature,solve_diffusion_coefficient_temperature_dependent,solve_with_decay,V,V0,data,BC_heat_transfer)
    time_stepping(solve_temperature,solve_diffusion,solve_diffusion_coefficient_temperature_dependent,solving_parameters,F,BC_diffusion,FT,BC_heat_transfer)



