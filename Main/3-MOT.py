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
import inspect
#from tqdm import *
#from random import random, randint
#from time import sleep
import math
from scipy import interpolate as scipy_interpolate
from collections import Iterable

from materials_properties import calculate_D, calculate_thermal_conductivity, calculate_specific_heat, calculate_density

def get_apreprovars(apreprovars):
    #return 'MOT_parameters_RCB.json'
    return 'Problems/Breeder_Blanket/Parameters/MOT_parameters_breeder_blankets.json'
    #return 'MOT_parameters_breeder_blankets_connected.json'
    #return 'MOT_parameters_CFD.json'


def byteify(input):
    if isinstance(input, dict):
        return {byteify(key): byteify(value)
                for key, value in input.iteritems()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input


def get_databases(name_database):
    with open(name_database) as f:
        data = json.load(f)
    #data = byteify(data)
    return data


def get_solvers(data):
    print('Getting the solvers')
    if data['solving_parameters']['study'] == "steady_state":
        solve_transient = False
    else:
        if data['solving_parameters']['study'] == 'transient':
            solve_transient = True
    if data['physics']['solve_heat_transfer'] == 1:
        solve_heat_transfer = True
    else:
        solve_heat_transfer = False
    if data['physics']['solve_tritium_diffusion'] == 1:
        solve_diffusion = True
    else:
        solve_diffusion = False
    if data['physics']['solve_laminar_flow'] == 1:
        solve_laminar_flow = True
    else:
        solve_laminar_flow = False
    if data['physics']['diffusion_coeff_temperature_dependent'] == 1:
        solve_diffusion_coefficient_temperature_dependent = True
    else:
        solve_diffusion_coefficient_temperature_dependent = False
    if data['physics']['solve_with_decay'] == 1:
        solve_with_decay = True
    else:
        solve_with_decay = False
    calculate_off_gassing = True
    return solve_transient, solve_heat_transfer, solve_laminar_flow, \
        solve_diffusion, solve_diffusion_coefficient_temperature_dependent, \
        solve_with_decay


def get_solving_parameters(data):
    print('Getting the solving parameters')
    Time = 0
    num_steps = 0
    dt = 0
    if solve_transient is True:
        Time = float(data["solving_parameters"]['final_time'])  #60000*365.25*24*3600.0# final time
        num_steps = data['solving_parameters']['number_of_time_steps']  #number of time steps
        dt = Time / num_steps  # time step size
    return Time, num_steps, dt


def define_mesh(data):
    print('Defining mesh')
    # Read in Mesh and markers from file
    mesh = Mesh()
    xdmf_in = XDMFFile(mesh.mpi_comm(), str(data['mesh_file']))
    xdmf_in.read(mesh)
    # prepare output file for writing by writing the mesh to the file
    # xdmf_out = XDMFFile(str(data['mesh_file']).split('.')[1]+'_from_fenics.xdmf')

    subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
    print('Number of cell is ' + str(len(subdomains.array())))
    n0 = FacetNormal(mesh)
    return mesh, xdmf_in, n0


def define_functionspaces(data):
    print('Defining Functionspaces')
    V = FunctionSpace(mesh, 'P', 1)  # FunctionSpace of the solution c
    V0 = FunctionSpace(mesh, 'DG', 0)  # FunctionSpace of the materials properties
    U = VectorFunctionSpace(mesh, 'P', 2) # FunctionSpace of velocity
    return V, V0, U


def get_surface_marker(mesh, xdmf_in):
    print('Marking the surfaces')
    # read in Surface markers for boundary conditions
    surface_marker_mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
    xdmf_in.read(surface_marker_mvc, "surface_marker")
    surface_marker_mvc.rename("surface_marker", "surface_marker")
    # xdmf_out.write(surface_marker_mvc, xdmf_encoding)
    surface_marker = MeshFunction("size_t", mesh, surface_marker_mvc)

    ds = Measure('ds', domain=mesh, subdomain_data=surface_marker)
    return surface_marker, ds


def define_BC_diffusion(data, solve_diffusion, V, surface_marker, ds):
    # #Tritium Diffusion
    bcs_c = list()
    Neumann_BC_c_diffusion = []
    Robin_BC_c_diffusion = []
    if solve_diffusion is True:
        print('Defining BC tritium diffusion')
        # DC
        for DC in data['physics']['tritium_diffusion']['boundary_conditions']['dc']:
            value_DC = Expression(str(DC['value']), t=0, degree=2)
            if type(DC['surface']) == list:
                for surface in DC['surface']:
                    bci_c = DirichletBC(V, value_DC, surface_marker, surface)
                    bcs_c.append(bci_c)
            else:
                bci_c = DirichletBC(V, value_DC, surface_marker, DC['surface'])
                bcs_c.append(bci_c)
        # Neumann
        for Neumann in data['physics']['tritium_diffusion']['boundary_conditions']['neumann']:
            value = Expression(str(Neumann['value']), t=0, degree=2)
            if type(Neumann['surface']) == list:
                for surface in Neumann['surface']:
                    Neumann_BC_c_diffusion.append([ds(surface), value])
            else:
                Neumann_BC_c_diffusion.append([ds(Neumann['surface']), value])
        # Robins
        for Robin in data['physics']['tritium_diffusion']['boundary_conditions']['robin']:
            value = Function(V)
            value = eval(Robin['value'])
            if type(Robin['surface']) == list:
                for surface in Robin['surface']:
                    Robin_BC_c_diffusion.append([ds(surface), value])
            else:
                Robin_BC_c_diffusion.append([ds(Robin['surface']), value])
    return bcs_c, Neumann_BC_c_diffusion, Robin_BC_c_diffusion


def define_BC_heat_transfer(data, solve_heat_transfer, V, surface_marker, ds):
    # #Temperature
    bcs_T = list()
    Neumann_BC_T_diffusion = []
    Robin_BC_T_diffusion = []

    if solve_heat_transfer is True:
        print('Defining BC heat transfer')
        # DC
        for DC in data['physics']['heat_transfers']['boundary_conditions']['dc']:
            value_DC = Expression(str(DC['value']), t=0, degree=2)
            if type(DC['surface'])==list:
                for surface in DC['surface']:
                    bci_T = DirichletBC(V, value_DC, surface_marker, surface)
                    bcs_T.append(bci_T)
                    # print(bci_T)
            else:
                # print(DC)
                bci_T = DirichletBC(V, value_DC, surface_marker, DC['surface'])
                bcs_T.append(bci_T)
                # print(bci_T)

        # Neumann
        for Neumann in data['physics']['heat_transfers']['boundary_conditions']['neumann']:
            value = Expression(str(Neumann['value']), t=0, degree=2)

            if type(Neumann['surface']) == list:
                for surface in Neumann['surface']:
                    Neumann_BC_T_diffusion.append([ds(surface), value])
            else:
                Neumann_BC_T_diffusion.append([ds(Neumann['surface']), value])

        # Robins
        for Robin in data['physics']['heat_transfers']['boundary_conditions']['robin']:

            if type(Robin['surface']) == list:
                for surface in Robin['surface']:
                    Robin_BC_T_diffusion.append([ds(surface), Robin['hc_coeff'], Robin['t_amb']])
            else:
                Robin_BC_T_diffusion.append([ds(Robin['surface']), Robin['hc_coeff'], Robin['t_amb']])

    return bcs_T, Neumann_BC_T_diffusion, Robin_BC_T_diffusion



def define_BC_laminar_flow(data, solve_laminar_flow, U, V, surface_marker, ds):
    print("Defining BC laminar flow")
    bcu = []
    for DC in data['physics']['laminar_flow']['boundary_conditions_velocity']['dc']:
        value = Expression((DC['valuex'],DC['valuey'],DC['valuez']),t=0,degree=2)
        for surface in DC['surface']:
            bci = DirichletBC(U, value, surface_marker, surface)
            bcu.append(bci)

    bcp = []
    for DC in data['physics']['laminar_flow']['boundary_conditions_pressure']['dc']:
        value = Expression(str(DC['value']), t=0, degree=2)
        for surface in DC['surface']:
            bci = DirichletBC(V, value, surface_marker, surface)
            bcp.append(bci)
    return bcu, bcp


def get_volume_markers(mesh):

    print('Marking the volumes')
    # read in the volume markers
    volume_marker_mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim())
    xdmf_in.read(volume_marker_mvc, "volume_marker_volume_id")

    volume_marker = MeshFunction("size_t", mesh, volume_marker_mvc)
    dx = Measure('dx', domain=mesh, subdomain_data=volume_marker)
    return volume_marker, dx


def define_initial_values(solve_heat_transfer, solve_diffusion, data, V):
    # #Tritium concentration
    c_n = Function(V)
    T_n = Function(V)
    if solve_diffusion is True:
        print('Defining initial values tritium diffusion')
        iniC = Expression(str(data['physics']['tritium_diffusion']['initial_value']), degree=2)
        c_n = interpolate(iniC, V)
    # #Temperature
    if solve_heat_transfer is True:
        print('Defining initial values heat transfer')
        iniT = Expression(str(data['physics']['heat_transfers']['initial_value']), degree=2)
        T_n = interpolate(iniT, V)
    return c_n, T_n


def define_source_terms(solve_heat_transfer, solve_diffusion, dx, data, sV):
    print('Defining source terms')
    Source_c_diffusion = list()
    Source_T_diffusion = list()

    if solve_diffusion is True:
        for source in data["physics"]["tritium_diffusion"]["source_terms"]:
            value = Expression(str(source["value"]), t=0, degree=2)
            for volume in source["volumes"]:
                Source_c_diffusion.append([dx(volume), value])

    if solve_heat_transfer is True:
        for source in data["physics"]["heat_transfers"]["source_terms"]:
            if str(source["value"]).endswith('.xml'):
                mesh_file = Mesh(source['mesh'])
                V_file = FunctionSpace(mesh_file, 'DG', 0)
                value_file = Function(V_file, source['value'])
                value = interpolate(value_file, V)
                File('vol_source.pvd') << value
            else:
                value = Expression(str(source["value"]), t=0, degree=2)
            for volume in source["volumes"]:
                Source_T_diffusion.append([dx(volume), value])

    return Source_c_diffusion, Source_T_diffusion


def calculate_D(T, material_id):
    R = 8.314  # Perfect gas constant
    k_B = 8.6e-5
    if material_id == "concrete":  # Concrete
        return 2e-6
    elif material_id == "polymer":  # Polymer
        return 2.0e-7*np.exp(-29000.0/R/T)
    elif material_id == "steel":  # steel
        return 7.3e-7*np.exp(-6.3e3/T)
    elif material_id == "tungsten":
        return 4.1e-7*np.exp(-0.39/k_B/T)
    elif material_id == "eurofer":
        return 8.1e-8*np.exp(-14470/R/T)
    elif material_id == "lithium_lead":
        return 2.5e-7*np.exp(-27000/R/T)
    else:
        raise ValueError("!!ERROR!! Unable to find "+str(material_id)+" as material ID in the database "+str(inspect.stack()[0][3]))



def calculate_thermal_conductivity(T, material_id):

    if material_id == "concrete":
        return 0.5

    elif material_id == "tungsten":
       # temperature_c =       [20, 100, 200, 300, 400, 500, 600, 700]
        temperature_k =        [293.15, 373.15, 473.15, 573.15, 673.15, 773.15, 873.15, 973.15]
        thermal_conductivity = [172.8,  164.8,  155.5,  147.2,  139.8,  133.1,  127.2,  122.1]

    elif material_id == "lithium_lead":
        #temperature_c =       [20,     300,    350,    400,    450,    500,    550,    600,    650,    700]
        temperature_k =        [293.15, 573.15, 623.15, 673.15, 723.15, 773.15, 823.15, 873.15, 923.15, 973.15]
        thermal_conductivity = [7.69,   13.18,  14.16,  15.14,  16.12,  17.10,  18.08,  19.06,  20.04,  21.02]

    elif material_id == "eurofer":
        #temperature_c =       [20,     50,     100,    150,    200, 250, 300, 350, 400, 450, 500, 550,600]
        temperature_k =        [293.15, 323.15, 373.15, 423.15, 473.15, 523.15, 573.15, 623.15, 673.15, 723.15, 773.15, 823.15, 873.15]
        thermal_conductivity = [27.63,  28.73,  29.87,  30.32,  30.28,  29.95,  29.51,  29.10,  28.84,  28.82,  29.08,  29.62,  30.38]

    else:

        raise ValueError("!!ERROR!! Unable to find "+str(material_id)+" as material ID in the database "+str(inspect.stack()[0][3]))

    interpolated_object = scipy_interpolate.interp1d(temperature_k, thermal_conductivity) # this object could be created once on inititation to speed up the code
    return float(interpolated_object.__call__(T))


def calculate_specific_heat(T, material_id):

    if material_id == "concrete":
        return 880

    elif material_id == "tungsten":
       # temperature_c = [20, 100, 200, 300, 400, 500, 600, 700]
        temperature_k =[293.15, 373.15, 473.15, 573.15, 673.15, 773.15, 873.15, 973.15]
        specific_heat = [129, 131.6, 134.7, 137.8, 140.9, 133.1, 127.2, 122.1]

    elif material_id == "lithium_lead":
        #temperature_c = [20, 300, 350, 400, 450, 500, 550, 600, 650, 700]
        temperature_k = [293.15, 573.15, 623.15, 673.15, 723.15, 773.15, 823.15, 873.15, 923.15, 973.15]
        specific_heat = [192, 190, 189, 189, 188, 188, 187, 187, 187, 186]

    elif material_id == "eurofer":
        #temperature_c = [20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600]
        temperature_k =  [293.15, 323.15, 373.15, 423.15, 473.15, 523.15, 573.15, 623.15, 673.15, 723.15, 773.15, 823.15, 873.15]
        specific_heat =[439, 462, 490, 509, 523, 534, 546, 562, 584, 616, 660, 721, 800]
    else:

        raise ValueError("!!ERROR!! Unable to find "+str(material_id)+" as material ID in the database "+str(inspect.stack()[0][3]))

    interpolated_object = scipy_interpolate.interp1d(temperature_k, specific_heat) # this object could be created once on inititation to speed up the code
    return float(interpolated_object.__call__(T))


def calculate_density(T, material_id):
    if material_id == "concrete":
        return 2400
    elif material_id == "tungsten":
       # temperature_c = [20, 100, 200, 300, 400, 500, 600, 700]
        temperature_k =[293.15, 373.15, 473.15, 573.15, 673.15, 773.15, 873.15, 973.15]
        density = [19298, 19279, 19254, 19229, 19205, 19178, 19152, 19125 ]

    elif material_id == "lithium_lead":
        #temperature_c = [20, 300, 350, 400, 450, 500, 550, 600, 650, 700]
        temperature_k = [293.15, 573.15, 623.15, 673.15, 723.15, 773.15, 823.15, 873.15, 923.15, 973.15]
        density = [10172, 9839, 9779, 9720, 9661, 9601, 9542, 9482, 9423, 9363]

    elif material_id == "eurofer":
        #temperature_c = [20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600]
        temperature_k = [293.15, 323.15, 373.15, 423.15, 473.15, 523.15, 573.15, 623.15, 673.15, 723.15, 773.15, 823.15, 873.15]
        density =       [7760,   7753,   7740,   7727,   7713,   7699,   7685,   7670,   7655,   7640,   7625,   7610, 7594]
    else:
        raise ValueError("!!ERROR!! Unable to find "+str(material_id)+" as material ID in the database "+str(inspect.stack()[0][3]))

    interpolated_object = scipy_interpolate.interp1d(temperature_k, density) # this object could be created once on inititation to speed up the code
    return float(interpolated_object.__call__(T))


def which_material_is_it(volume_id, data):
    for material in data["structure_and_materials"]["materials"]:
        for volumes in material["volumes"]:
            if volume_id in [volumes]:
                material_id = material["material"]
                break
    return material_id


def define_materials_properties(V0, data, volume_marker):
    # ##Defining materials properties
    print('Defining the materials properties')

    D = Function(V0)  # Diffusion coefficient
    thermal_conductivity = Function(V0)
    specific_heat = Function(V0)
    density = Function(V0)
    # #Assigning each to each cell its properties
    for cell_no in range(len(volume_marker.array())):
        volume_id = volume_marker.array()[cell_no]  # This is the volume id (Trelis)
        material_id = which_material_is_it(volume_id, data)
        if solve_heat_transfer is True:
            thermal_conductivity.vector()[cell_no] = calculate_thermal_conductivity(data['physics']['heat_transfers']['initial_value'], material_id)
            density.vector()[cell_no] = calculate_density(data['physics']['heat_transfers']['initial_value'], material_id)
            specific_heat.vector()[cell_no]=calculate_specific_heat(data['physics']['heat_transfers']['initial_value'], material_id)
        if solve_diffusion is True:
            D.vector()[cell_no] = calculate_D(data['physics']['heat_transfers']['initial_value'], material_id)

        # print(D.vector()[cell_no])
    return D, thermal_conductivity, specific_heat, density


def define_variational_problem_diffusion(solve_diffusion, solve_transient, solve_with_decay, V, Neumann_BC_c_diffusion, Robin_BC_c_diffusion, Source_c_diffusion, data):
    if solve_diffusion is True:
        print('Defining variation problem tritium diffusion')
        c = TrialFunction(V)  # c is the tritium concentration
        vc = TestFunction(V)
        if solve_with_decay is True:
            decay = np.log(2)/(12.33*365.25*24*3600)  # Tritium Decay constant [s-1]
        else:
            decay = 0

        if solve_transient is True:
            F = ((c-c_n)/dt)*vc*dx
        else:
            F = 0
        for source in Source_c_diffusion:
            F += -vc*source[1]*source[0]

        F += D*dot(grad(c), grad(vc))*dx + decay*c*vc*dx
        for Neumann in Neumann_BC_c_diffusion:
            F += vc * Neumann[1]*Neumann[0]
        for Robin in Robin_BC_c_diffusion:
            F += D*vc*Robin[1]*Robin[0]
        return F
    return False, False


def define_variational_problem_heat_transfer(solve_heat_transfer, solve_transient, V, Neumann_BC_T_diffusion, Robin_BC_T_diffusion, Source_T_diffusion, data):

    if solve_heat_transfer is True:
        print('Defining variation problem heat transfer')
        T = TrialFunction(V)  # T is the temperature
        vT = TestFunction(V)
        if solve_transient is True:
            FT = specific_heat*density*((T-T_n)/dt)*vT*dx
        else:
            FT = 0
        for source in Source_T_diffusion:
            FT += -vT*source[1]*source[0]
        FT += thermal_conductivity*dot(grad(T), grad(vT))*dx  # This is the heat transfer equation
        for Neumann in Neumann_BC_T_diffusion:
            FT += - vT * Neumann[1]*Neumann[0]
        for Robin in Robin_BC_T_diffusion:
            FT += vT * Robin[1] * (T-Robin[2])*Robin[0]
        return FT
    return False, False


def define_variational_problem_laminar_flow(solve_laminar_flow, solve_transient, U, V, data):

    if solve_laminar_flow is True:
        print('Defining variation problem laminar flow')
        mu = 1.875e-4
        density = 11600
        # Define trial and test functions
        u = TrialFunction(U)
        p = TrialFunction(V)
        v = TestFunction(U)
        q = TestFunction(V)

        u1 = Function(U)
        u0 = Function(U)
        p1 = Function(V)
        # Define coefficients
        k = Constant(dt)
        f = Constant((0, 0, 0))

        # Tentative velocity step
        F1 = density*(1/k)*inner(u - u0, v)*dx + density*inner(grad(u0)*u0, v)*dx + \
             mu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
        a1 = lhs(F1)
        L1 = rhs(F1)

        # Pressure update
        a2 = inner(grad(p), grad(q))*dx
        L2 = -(1/k)*div(u1)*q*dx

        # Velocity update
        a3 = inner(u, v)*dx
        L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx


        # Assemble matrices
        A1 = assemble(a1)
        A2 = assemble(a2)
        A3 = assemble(a3)
        return A1, L1, A2, L2, A3, L3, u1, u0, p1

    return False, False, False, False, False, False, False, False, False


def update_D(mesh, volume_marker, D, T):
    for cell in cells(mesh):
        cell_no = cell.index()
        material_id = volume_marker.array()[cell_no]
        Ta = 0
        for i in range(0, 12, 3):
            Ta += T(cell.get_vertex_coordinates()[i], cell.get_vertex_coordinates()[i+1], cell.get_vertex_coordinates()[i+2])
        Ta = Ta/4
        D_value = calculate_D(Ta, material_id)
        D.vector()[cell_no] = D_value  # Assigning for each cell the corresponding diffusion coeff


def update_bc(t, physic):
    if physic == "tritium_diffusion":
        for Neumann in Neumann_BC_c_diffusion:
            Neumann[1].t = t
    if physic == "heat_transfers":
        for Neumann in Neumann_BC_T_diffusion:
            Neumann[1].t = t

    bcs = list()
    for DC in data['physics'][physic]['boundary_conditions']['dc']:
        value_DC = Expression(str(DC['value']), t=t, degree=2)
        if type(DC['surface']) == list:
            for surface in DC['surface']:
                bci = DirichletBC(V, value_DC, surface_marker, surface)
                bcs.append(bci)
        else:
            print(DC)
            bci = DirichletBC(V, float(value_DC), surface_marker, DC['surface'])
            bcs.append(bci)
    return bcs


def update_source_term(t, physic):
    if physic == "tritium_diffusion":
        for source in Source_c_diffusion:
            source[1].t = t
    if physic == "heat_transfers":
        for source in Source_T_diffusion:
            source[1].t = t
    return


def calculate_flux(surface, ds, solution, property, n0):
    flux = assemble(-property*dot(grad(solution), n0)*ds(surface))
    return flux


def calculate_average_volume(volume, dx, solution, property):
    average = assemble(solution*dx(volume))/assemble(1*dx(volume))
    return average


def write_in_file(header, values, output_file):
    with open(output_file, "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerows([header])

        for val in values:
            writer.writerows([val])
    return


def calculate_minimum_volume(u, cell_function, subdomain_id, V):
    dofmap = V.dofmap()
    mesh = V.mesh()
    mini = u.vector().max()
    for cell in cells(mesh):
        if cell_function[cell.index()] == subdomain_id:
            dofs = dofmap.cell_dofs(cell.index())
            for dof in dofs:
                try:
                    [dof][0]
                    if u.vector()[dof][0] < mini:
                        mini = u.vector()[dof][0]
                except:
                    if u.vector()[dof] < mini:
                        mini = u.vector()[dof]
    return mini


def calculate_maximum_volume(u, cell_function, subdomain_id, V):
    dofmap = V.dofmap()
    mesh = V.mesh()
    maxi = u.vector().min()
    for cell in cells(mesh):
        if cell_function[cell.index()] == subdomain_id:
            dofs = dofmap.cell_dofs(cell.index())
            for dof in dofs:
                try:
                    [dof][0]
                    if u.vector()[dof][0] > maxi:
                        maxi = u.vector()[dof][0]
                except:
                    if u.vector()[dof] > maxi:
                        maxi = u.vector()[dof]
    return maxi


def post_processing(data, solution, physic, header, values, t, ds, dx, volume_marker, property, n0):

    output_file = data["post_processing"][physic]["output_file"]
    tab = []
    tab.append(t)
    for surface in data["post_processing"][physic]["surface_flux"]:
        tab.append(calculate_flux(surface, ds, solution, property, n0))
    for volume in data["post_processing"][physic]["volume_average"]:
        tab.append(calculate_average_volume(volume, dx, solution, property))
    for volume in data["post_processing"][physic]["volume_minimum"]:
        tab.append(calculate_minimum_volume(solution, volume_marker, volume, V))
    for volume in data["post_processing"][physic]["volume_maximum"]:
        tab.append(calculate_maximum_volume(solution, volume_marker, volume, V))
    for expression in data["post_processing"][physic]["custom"]:
        tab.append(eval(expression))
    values.append(tab)

    write_in_file(header, values, output_file)
    return


def initialise_post_processing(data, physic):
    output_file = data["post_processing"][physic]["output_file"]
    header = ['t(s)']
    i = 0
    for surface in data["post_processing"][physic]["surface_flux"]:
        header.append('Flux through surface '+str(surface))
    for volume in data["post_processing"][physic]["volume_average"]:
        header.append('Average volume '+str(volume))
    for volume in data["post_processing"][physic]["volume_minimum"]:
        header.append('Minimum volume '+str(volume))
    for volume in data["post_processing"][physic]["volume_maximum"]:
        header.append('Maximum volume '+str(volume))
    for expression in data["post_processing"][physic]["custom"]:
        i+=1
        header.append('Custom ' + str(i))
    return header


def time_stepping(data,
                  solve_heat_transfer,
                  solve_diffusion,
                  solve_laminar_flow,
                  solve_diffusion_coefficient_temperature_dependent,
                  Time,
                  num_steps,
                  dt,
                  V,
                  D,
                  thermal_conductivity,
                  F,
                  bcs_c,
                  FT,
                  q,
                  bcs_T,
                  ds,
                  dx,
                  volume_marker,
                  header_heat_transfers,
                  header_tritium_diffusion,
                  values_heat_transfers,
                  values_tritium_diffusion):
    # pbar = tqdm(total=100,bar_format='{desc}: {percentage:3.0f}%|{bar}|{n:.0f}/{total_fmt} [{elapsed}<{remaining}, ' '{rate_fmt}{postfix}]')
    ## Time-stepping
    try:
        set_log_active(False)
    except:
        print('active log not available')
    print('Time stepping')
    T = Function(V)
    c = Function(V)

    off_gassing = list()
    output_file = File('solution.pvd')#File(data["output_file"])
    t = 0
    # Use amg preconditioner if available
    prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"
    # Use nonzero guesses - essential for CG with non-symmetric BC
    parameters['krylov_solver']['nonzero_initial_guess'] = True

    for n in range(num_steps):
        t += dt

        print('{0}% {1} \r'.format(100*t/Time, '|'*int(50*t/Time)),)


        # Compute solution concentration
        if solve_diffusion is True:
            update_source_term(t, 'tritium_diffusion')
            solve(lhs(F) == rhs(F), c, bcs_c)
            output_file << (c, t)
            c_n.assign(c)
            bcs_c = update_bc(t, "tritium_diffusion")
            post_processing(data, c, "tritium_diffusion", header_tritium_diffusion, values_tritium_diffusion, t, ds, dx, volume_marker, D, n0)

        # Compute solution temperature
        if solve_heat_transfer is True:
            update_source_term(t, 'heat_transfers')
            solve(lhs(FT) == rhs(FT), T, bcs_T)
            output_file << (T, t)
            T_n.assign(T)
            bcs_T = update_bc(t, "heat_transfers")
            post_processing(data, T, "heat_transfers", header_heat_transfers, values_heat_transfers, t, ds, dx, volume_marker, thermal_conductivity, n0)

        # Compute solution velocity and pressure
        if solve_laminar_flow is True:

            # Compute tentative velocity step
            begin("Computing tentative velocity")
            b1 = assemble(L1)
            [bc.apply(A1, b1) for bc in bcu]
            solve(A1, u1.vector(), b1, "bicgstab", "default")
            end()

            # Pressure correction
            begin("Computing pressure correction")
            b2 = assemble(L2)
            [bc.apply(A2, b2) for bc in bcp]
            [bc.apply(p1.vector()) for bc in bcp]
            solve(A2, p1.vector(), b2, "bicgstab", prec)
            end()

            # Velocity correction
            begin("Computing velocity correction")
            b3 = assemble(L3)
            [bc.apply(A3, b3) for bc in bcu]
            solve(A3, u1.vector(), b3, "bicgstab", "default")
            end()
            output_file << u1
            #output_file << (p1,t)
            u0.assign(u1)

        # Update the materials properties
        if solve_diffusion_coefficient_temperature_dependent is True and solve_heat_transfer is True and solve_diffusion is True:
            D = update_D(mesh, volume_marker, D, T)
    return


def solving(data,
            solve_heat_transfer, solve_diffusion, solve_laminar_flow, solve_diffusion_coefficient_temperature_dependent, Time, num_steps, dt, V, D, thermal_conductivity, F, Source_c_diffusion, bcs_c, FT, Source_T_diffusion, bcs_T, ds, dx, volume_marker, n0):
    values_heat_transfers = []
    values_tritium_diffusion = []
    header_heat_transfers = ''
    header_tritium_diffusion = ''
    if solve_heat_transfer is True:
        header_heat_transfers = initialise_post_processing(data, "heat_transfers")
    if solve_diffusion is True:
        header_tritium_diffusion = initialise_post_processing(data, "tritium_diffusion")

    if solve_transient is True:
        time_stepping(data=data,
                      solve_heat_transfer=solve_heat_transfer,
                      solve_diffusion=solve_diffusion,
                      solve_laminar_flow=solve_laminar_flow,
                      solve_diffusion_coefficient_temperature_dependent=solve_diffusion_coefficient_temperature_dependent,
                      Time=Time,
                      num_steps=num_steps,
                      dt=dt,
                      V=V,
                      D=D,
                      thermal_conductivity=thermal_conductivity,
                      F=F,
                      bcs_c=bcs_c,
                      FT=FT,
                      q=Source_T_diffusion,
                      bcs_T=bcs_T,
                      ds=ds,
                      dx=dx,
                      volume_marker=volume_marker,
                      header_heat_transfers=header_heat_transfers,
                      header_tritium_diffusion=header_tritium_diffusion,
                      values_heat_transfers=values_heat_transfers,
                      values_tritium_diffusion=values_tritium_diffusion)


    else:
        output_file = File(data["output_file"])
        if solve_heat_transfer is True:
            T = Function(V)
            solve(lhs(FT) == rhs(FT), T, bcs_T)
            output_file << (T, 0.0)
            post_processing(data, T, "heat_transfers", header_heat_transfers, values_heat_transfers, 0, ds, dx, volume_marker, thermal_conductivity, n0)
        if solve_diffusion is True:
          c = Function(V)
          solve(lhs(F) == rhs(F), c, bcs_c)
          output_file << (c, 0.0)
          post_processing(data, c, "tritium_diffusion", header_tritium_diffusion, values_tritium_diffusion, 0, ds, dx, volume_marker, D, n0)



if __name__ == "__main__":

    apreprovars = get_apreprovars(2)
    data = get_databases(apreprovars)  # This returns an object data=json.load()
    solve_transient, solve_heat_transfer, solve_laminar_flow, solve_diffusion, solve_diffusion_coefficient_temperature_dependent, solve_with_decay = get_solvers(data)  # Gets the solvers

    Time, num_steps, dt = get_solving_parameters(data)  # Gets the parameters (final time, time steps...)

    mesh, xdmf_in, n0 = define_mesh(data)

    V, V0, U = define_functionspaces(data)

    volume_marker, dx = get_volume_markers(mesh)#, xdmf_in)

    c_n, T_n = define_initial_values(solve_heat_transfer, solve_diffusion, data, V)

    Source_c_diffusion, Source_T_diffusion = define_source_terms(solve_heat_transfer, solve_diffusion, dx, data, V)

    surface_marker, ds = get_surface_marker(mesh, xdmf_in)

    bcs_c, Neumann_BC_c_diffusion, Robin_BC_c_diffusion = define_BC_diffusion(data, solve_diffusion, V, surface_marker, ds)

    bcs_T, Neumann_BC_T_diffusion, Robin_BC_T_diffusion = define_BC_heat_transfer(data, solve_heat_transfer, V, surface_marker, ds)

    bcu, bcp = define_BC_laminar_flow(data,solve_laminar_flow, U, V, surface_marker, ds)

    D, thermal_conductivity, specific_heat, density = define_materials_properties(V0, data, volume_marker)

    F = define_variational_problem_diffusion(solve_diffusion, solve_transient, solve_with_decay, V, Neumann_BC_c_diffusion, Robin_BC_c_diffusion, Source_c_diffusion, data)

    FT = define_variational_problem_heat_transfer(solve_heat_transfer, solve_transient, V, Neumann_BC_T_diffusion, Robin_BC_T_diffusion, Source_T_diffusion, data)

    A1, L1, A2, L2, A3, L3, u1, u0, p1 = define_variational_problem_laminar_flow(solve_laminar_flow, solve_transient, U, V, data)

    solving(data, solve_heat_transfer, solve_diffusion, solve_laminar_flow, solve_diffusion_coefficient_temperature_dependent, Time, num_steps, dt, V, D, thermal_conductivity, F, Source_c_diffusion, bcs_c, FT, Source_T_diffusion, bcs_T, ds, dx, volume_marker, n0)
