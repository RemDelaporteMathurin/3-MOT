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

from materials_properties import calculate_D, calculate_thermal_conductivity, calculate_specific_heat, calculate_density, calculate_mu


def get_apreprovars(apreprovars):
    #return 'Problems/RCB/Parameters/MOT_parameters_RCB.json'
    return 'Problems/Breeder_Blanket/Parameters/MOT_parameters_breeder_blankets.json'
    #return 'Problems/Square_Pipe/Parameters/MOT_parameters_square_pipe.json'
    #return 'MOT_parameters_breeder_blankets_connected.json'


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

    update_thermal_properties = False
    update_tritium_diffusion_properties = False
    update_fluid_properties = False
    if data['solving_parameters']['study'] == "steady_state":
        solve_transient = False
    else:
        if data['solving_parameters']['study'] == 'transient':
            solve_transient = True
    if data['physics']['solve_heat_transfer'] == 1:
        solve_heat_transfer = True
        if data['physics']['heat_transfers']['update_properties'] == 1:
            update_thermal_properties = True
        else:
            update_thermal_properties = False
    else:
        solve_heat_transfer = False
    if data['physics']['solve_tritium_diffusion'] == 1:
        solve_diffusion = True
        if data['physics']['tritium_diffusion']['update_properties'] == 1:
            update_tritium_diffusion_properties = True
        else:
            update_tritium_diffusion_properties = False
    else:
        solve_diffusion = False
    if data['physics']['solve_laminar_flow'] == 1:
        solve_laminar_flow = True
        if data['physics']['laminar_flow']['update_properties'] == 1:
            update_fluid_properties = True
        else:
            update_fluid_properties = False
    else:
        solve_laminar_flow = False
    if data['physics']['couple_tritium_diffusion_heat_transfer'] == 1:
        couple_tritium_diffusion_heat_transfer = True
    else:
        couple_tritium_diffusion_heat_transfer = False
    if data['physics']['solve_with_decay'] == 1:
        solve_with_decay = True
    else:
        solve_with_decay = False
    if data['physics']['couple_tritium_diffusion_laminar_flow'] == 1:
        couple_tritium_diffusion_laminar_flow = True
    else:
        couple_tritium_diffusion_laminar_flow = False
    if data['physics']['couple_heat_transfer_laminar_flow'] == 1:
        couple_heat_transfer_laminar_flow = True
    else:
        couple_heat_transfer_laminar_flow = False
    return solve_transient, solve_heat_transfer, solve_laminar_flow, \
        solve_diffusion, solve_with_decay, \
        couple_tritium_diffusion_heat_transfer, \
        couple_heat_transfer_laminar_flow, \
        couple_tritium_diffusion_laminar_flow, update_thermal_properties, \
        update_tritium_diffusion_properties, update_fluid_properties


def get_solving_parameters(data, solve_transient):
    print('Getting the solving parameters')
    Time = 0
    num_steps = 0
    dt = 0
    if solve_transient is True:
        Time = float(data["solving_parameters"]['final_time'])  #60000*365.25*24*3600.0# final time
        num_steps = data['solving_parameters']['number_of_time_steps']  #number of time steps
        dt = Time / num_steps  # time step size
    return Time, num_steps, dt


def get_volume_markers(mesh, xdmf_in):

    print('Marking the volumes')
    # read in the volume markers
    volume_marker_mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim())
    xdmf_in.read(volume_marker_mvc, "volume_marker_volume_id")

    volume_marker = MeshFunction("size_t", mesh, volume_marker_mvc)
    dx = Measure('dx', domain=mesh, subdomain_data=volume_marker)
    return volume_marker, dx


def get_surface_marker(mesh, xdmf_in):
    print('Marking the surfaces')
    # read in Surface markers for boundary conditions
    surface_marker_mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
    xdmf_in.read(surface_marker_mvc, "surface_marker")
    surface_marker_mvc.rename("surface_marker", "surface_marker")
    #xdmf_out = XDMFFile(str(data['mesh_file']).split('.')[1]+'_from_fenics.xdmf')
    #xdmf_out.write(surface_marker_mvc)
    surface_marker = MeshFunction("size_t", mesh, surface_marker_mvc)

    ds = Measure('ds', domain=mesh, subdomain_data=surface_marker)
    return surface_marker, ds


def define_mesh(data, solve_laminar_flow):
    print('Defining mesh')
    # Read in Mesh and markers from file
    mesh = Mesh()
    xdmf_in = XDMFFile(mesh.mpi_comm(), str(data['mesh_file']))
    xdmf_in.read(mesh)
    # prepare output file for writing by writing the mesh to the file
    #xdmf_out = XDMFFile(str(data['mesh_file']).split('.')[1]+'_from_fenics.xdmf')
    #xdmf_out.write()
    ncells = MeshFunction("size_t", mesh, mesh.topology().dim())
    print('Number of cell is ' + str(len(ncells.array())))
    n0 = FacetNormal(mesh)

    volume_marker, dx = get_volume_markers(mesh, xdmf_in)
    surface_marker, ds = get_surface_marker(mesh, xdmf_in)

    if solve_laminar_flow is True:
        domains_fluid = MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
        
        for volume in data["physics"]["laminar_flow"]["volumes"]:
            domains_fluid.array()[volume_marker.array() == volume] = 1
        mesh_fluid = SubMesh(mesh, domains_fluid, 1)
        surface_marker_fluid = MeshFunction("size_t", mesh_fluid, mesh_fluid.topology().dim() - 1, 0)
        ncells_fluid = MeshFunction("size_t", mesh_fluid, mesh_fluid.topology().dim())
        volume_marker_fluid = MeshFunction("size_t", mesh_fluid, mesh_fluid.topology().dim())
        print('Creating surface markers for SubMesh')
        vmap = mesh_fluid.data().array('parent_vertex_indices', 0)
        cmap = mesh_fluid.data().array('parent_cell_indices', mesh_fluid.topology().dim())

        n = 0
        for c in cells(mesh_fluid):
            print(str(100*n/len(ncells_fluid))+' %', end='\r')
            parent_cell = Cell(mesh, cmap[c.index()])
            volume_marker_fluid.array()[c.index()] = volume_marker.array()[parent_cell.index()]

            for f in facets(parent_cell):
                for g in facets(c):
                    g_vertices = vmap[g.entities(0)]
                    if set(f.entities(0)) == set(g_vertices):
                        surface_marker_fluid.array()[g.index()] = surface_marker.array()[f.index()]
            n += 1

        dx_fluid = Measure('dx', domain=mesh_fluid)
        ds_fluid = Measure('ds', domain=mesh_fluid, subdomain_data=surface_marker_fluid)
        n_fluid = FacetNormal(mesh_fluid)
    else:
        mesh_fluid = False
        dx_fluid = False
        ds_fluid = False
        n_fluid = False
        surface_marker_fluid = False
        volume_marker_fluid = False

    return mesh, n0, volume_marker, dx, surface_marker, ds, mesh_fluid, volume_marker_fluid, dx_fluid, surface_marker_fluid, ds_fluid, n_fluid


def define_functionspaces(data, mesh, mesh_fluid):
    print('Defining Functionspaces')
    V = FunctionSpace(mesh, 'P', 1)  # FunctionSpace of the solution c and T
    V0 = FunctionSpace(mesh, 'DG', 0)  # FunctionSpace of the materials properties
    if mesh_fluid is not False:
        Q = FunctionSpace(mesh_fluid, 'P', 1)  # Functionspace of pressure
        V2 = VectorFunctionSpace(mesh_fluid, 'P', 2)  # FunctionSpace of velocity
        V3  = VectorElement('P', mesh_fluid.ufl_cell(), 2)
        P  = FiniteElement('P', mesh_fluid.ufl_cell(), 1)
        TH = MixedElement([V3, P])
        W  = FunctionSpace(mesh_fluid, TH)
    else:
        Q = False#FunctionSpace(mesh, 'P', 1)  # Functionspace of pressure
        V2 = False#VectorFunctionSpace(mesh, 'P', 2)  # FunctionSpace of velocity
        V3  = False#VectorElement('P', mesh.ufl_cell(), 2)
        P  = False#FiniteElement('P', mesh.ufl_cell(), 1)
        TH = False#MixedElement([V3, P])
        W  = False#FunctionSpace(mesh, TH)
        
    return V, V0, V2, Q, W


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


def define_BC_laminar_flow(data, solve_laminar_flow, U, Q, W, surface_marker_fluid):
    if solve_laminar_flow is True:
        print("Defining BC laminar flow")
        bcu_transient = []
        bcs_stationary = []
        for DC in data['physics']['laminar_flow']['boundary_conditions_velocity']:
            value = Expression((DC['valuex'], DC['valuey'], DC['valuez']), t=0, degree=2)
            for surface in DC['surface']:
                bci = DirichletBC(U, value, surface_marker_fluid, surface)
                bci_  = DirichletBC(W.sub(0), value, surface_marker_fluid, surface)
                bcu_transient.append(bci)
                bcs_stationary.append(bci_)

        bcp_transient = []
        for DC in data['physics']['laminar_flow']['boundary_conditions_pressure']:
            value = Expression(str(DC['value']), t=0, degree=2)
            for surface in DC['surface']:
                bci = DirichletBC(Q, value, surface_marker_fluid, surface)
                bci_  = DirichletBC(W.sub(1), value, surface_marker_fluid, surface)
                bcp_transient.append(bci)
                bcs_stationary.append(bci_)

        
        return bcs_stationary, bcu_transient, bcp_transient
    else:
        return False, False, False


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


def which_material_is_it(volume_id, data):
    for material in data["structure_and_materials"]["materials"]:
        for volumes in material["volumes"]:
            if volume_id in [volumes]:
                material_id = material["material"]
                break
    return material_id


def define_materials_properties(V0, data, volume_marker, solve_heat_transfer, solve_diffusion, solve_laminar_flow):
    # ##Defining materials properties
    print('Defining the materials properties')

    D = Function(V0)  # Diffusion coefficient
    thermal_conductivity = Function(V0)
    specific_heat = Function(V0)
    density = Function(V0)
    mu = Function(V0)
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
        if solve_laminar_flow is True:
            mu.vector()[cell_no] = calculate_mu(data['physics']['heat_transfers']['initial_value'], material_id)
        # print(D.vector()[cell_no])
    return D, thermal_conductivity, specific_heat, density, mu


def define_variational_problem_diffusion(solve_diffusion, solve_transient, dt, solve_with_decay, V, Neumann_BC_c_diffusion, Robin_BC_c_diffusion, Source_c_diffusion):
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


def define_variational_problem_heat_transfer(solve_heat_transfer, solve_transient, dt, V, specific_heat, density, thermal_conductivity, Neumann_BC_T_diffusion, Robin_BC_T_diffusion, Source_T_diffusion):

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
        
        if couple_heat_transfer_laminar_flow is True:

            for volume in data["physics"]["laminar_flow"]["volumes"]:
                FT += specific_heat*density*(dot(u_, grad(T)))*vT*dx(2)

        return FT
    return False, False


def define_variational_problem_laminar_flow(solve_laminar_flow, solve_transient, V2, Q, W, dx_fluid, ds_fluid, n_fluid, bcu_transient, bcp_transient, bcs_stationary, dt, mu, density):

    if solve_laminar_flow is True:
        print('Defining variation problem laminar flow')
        if solve_transient is True:
            # Define strain-rate tensor
            def epsilon(u):
                return sym(nabla_grad(u))
            # Define stress tensor
            def sigma(u, p):
                return 2*mu*epsilon(u) - p*Identity(len(u))

            #### Define variational problem for step 1
            F1 = density*dot((u - u_n) / k, v)*dx_fluid \
                + density*dot(dot(u_n, nabla_grad(u_n)), v)*dx_fluid \
                + inner(sigma(U, p_n), epsilon(v))*dx_fluid \
                + dot(p_n*n_fluid, v)*ds_fluid - dot(mu*nabla_grad(U)*n_fluid, v)*ds_fluid \
                - dot(f, v)*dx_fluid
            a1 = lhs(F1)
            L1 = rhs(F1)
            # Define variational problem for step 2
            a2 = dot(nabla_grad(p), nabla_grad(q))*dx_fluid
            L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx_fluid - (1/k)*div(u_)*q*dx_fluid
            # Define variational problem for step 3
            a3 = dot(u, v)*dx_fluid
            L3 = dot(u_, v)*dx_fluid - k*dot(nabla_grad(p_ - p_n), v)*dx_fluid
            # Assemble matrices
            A1 = assemble(a1)
            A2 = assemble(a2)
            A3 = assemble(a3)

            # Apply boundary conditions to matrices
            [bc.apply(A1) for bc in bcu_transient]
            [bc.apply(A2) for bc in bcp_transient]

            return False, A1, L1, A2, L2, A3, L3
        else:
            (v_, q_) = TestFunctions(W)
            w = Function(W)
            (u, p) = split(w)
            ### Steady Part of the Momentum Equation
            def steady(u, dx_fluid):
                T = -p*I + 2*mu*sym(grad(u))
                return (inner(grad(u)*u, v_) + inner(T, grad(v_)) - inner(f, v_)) * dx_fluid


            I = Identity(u.geometric_dimension())

            F = steady(u, dx_fluid) + q_*div(u)*dx_fluid

            J = derivative(F, w)
            problem = NonlinearVariationalProblem(F, w, bcs_stationary, J)
            solver = NonlinearVariationalSolver(problem)
            return solver, False, False, False, False, False, False
    return False, False, False, False, False, False, False


def update_properties(mesh, mesh_fluid, volume_marker, volume_marker_fluid, T, solve_diffusion, solve_heat_transfer, solve_laminar_flow, update_tritium_diffusion_properties, update_thermal_properties, update_fluid_properties):
    for cell in cells(mesh):
        cell_no = cell.index()
        material_id = which_material_is_it(volume_marker.array()[cell_no], data)
        Ta = 0
        for i in range(0, 12, 3):
            Ta += T(cell.get_vertex_coordinates()[i], cell.get_vertex_coordinates()[i+1], cell.get_vertex_coordinates()[i+2])
        Ta = Ta/4
        if solve_diffusion is True and update_tritium_diffusion_properties is True:
            D_value = calculate_D(Ta, material_id)
            D.vector()[cell_no] = D_value  # Assigning for each cell the corresponding diffusion coeff
        if solve_heat_transfer is True and update_thermal_properties is True:
                thermal_conductivity_value = calculate_thermal_conductivity(Ta, material_id)
                thermal_conductivity.vector()[cell_no] = thermal_conductivity_value            
                specific_heat_value = calculate_specific_heat(Ta, material_id)
                specific_heat.vector()[cell_no] = specific_heat_value 
        if solve_laminar_flow is True and update_fluid_properties is True:
            mu_value = calculate_mu(Ta, material_id)
            mu.vector()[cell_no] = mu_value     
        if (solve_heat_transfer is True and update_thermal_properties is True):
                density_value = calculate_density(Ta, material_id)
                density.vector()[cell_no] = density_value         
    
    if solve_laminar_flow is True:
        for cell_fluid in cells(mesh_fluid):
            cell_no = cell_fluid.index()
            cmap = mesh_fluid.data().array('parent_cell_indices', mesh_fluid.topology().dim())
            cell_parent = Cell(mesh, cmap[cell_no])
            fluid_id = which_material_is_it(volume_marker_fluid.array()[cell_no], data)

            Ta = 0
            for i in range(0, 12, 3):
                Ta += T(cell_parent.get_vertex_coordinates()[i], cell_parent.get_vertex_coordinates()[i+1], cell_parent.get_vertex_coordinates()[i+2])
            Ta = Ta/4

            if update_fluid_properties is True :
                density_value = calculate_density(Ta, material_id)
                density.vector()[cell_no] = density_value
                mu_value = calculate_mu(Ta, material_id)
                mu.vector()[cell_no] = mu_value     


    return D, thermal_conductivity, density, specific_heat, mu


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
            bci = DirichletBC(V, float(value_DC), surface_marker, DC['surface'])
            bcs.append(bci)
    return bcs


def update_source_term(t, physic, Source_c_diffusion, Source_T_diffusion):
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
                  couple_tritium_diffusion_heat_transfer,
                  Time,
                  num_steps,
                  dt,
                  V,
                  U,
                  Q,
                  D,
                  thermal_conductivity,
                  mu,
                  density,
                  F,
                  bcs_c,
                  FT,
                  bcs_T,
                  A1, L1, A2, L2, A3, L3,
                  bcu_transient, bcp_transient,
                  ds,
                  dx,
                  Source_c_diffusion, Source_T_diffusion,
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
    output_file = File(data["output_file"])
    output_file1 = File('Problems/Square_pipe/solution_flow.pvd')
    t = 0
    # Use amg preconditioner if available
    #prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"
    # Use nonzero guesses - essential for CG with non-symmetric BC
    #parameters['krylov_solver']['nonzero_initial_guess'] = True

    for n in range(num_steps):
        t += dt
        print(100*t/Time, end=' % \r')
        # Compute solution velocity and pressure

        
        # Compute solution concentration
        if solve_diffusion is True:
            update_source_term(t, 'tritium_diffusion', Source_c_diffusion, Source_T_diffusion)
            solve(lhs(F) == rhs(F), c, bcs_c)
            output_file << (c, t)
            c_n.assign(c)
            bcs_c = update_bc(t, "tritium_diffusion")
            post_processing(data, c, "tritium_diffusion", header_tritium_diffusion, values_tritium_diffusion, t, ds, dx, volume_marker, D, n0)
        # Compute solution temperature

        if solve_laminar_flow is True:
            # Step 1: Tentative velocity step
            b1 = assemble(L1)
            [bc.apply(b1) for bc in bcu_transient]
            solve(A1, u_.vector(), b1)
            # Step 2: Pressure correction step
            b2 = assemble(L2)
            [bc.apply(b2) for bc in bcp]
            solve(A2, p_.vector(), b2)
            # Step 3: Velocity correction step
            b3 = assemble(L3)
            solve(A3, u_.vector(), b3)
            output_file << (u_, t)
            output_file << (p_, t)
            # Update previous solution
            u_n.assign(u_)
            p_n.assign(p_)
            print(u_(0,0,0))
        if solve_heat_transfer is True:
            update_source_term(t, 'heat_transfers', Source_c_diffusion, Source_T_diffusion)
            solve(lhs(FT) == rhs(FT), T, bcs_T)
            output_file << (T, t)
            T_n.assign(T)
            bcs_T = update_bc(t, "heat_transfers")
            post_processing(data, T, "heat_transfers", header_heat_transfers, values_heat_transfers, t, ds, dx, volume_marker, thermal_conductivity, n0)
        # Update the materials properties
        D, thermal_conductivity, density, specific_heat, mu = \
        update_properties(mesh, mesh_fluid, volume_marker, volume_marker_fluid, T, solve_diffusion, solve_heat_transfer, solve_laminar_flow, update_tritium_diffusion_properties, update_thermal_properties, update_fluid_properties)
    return


def solving(data, 
            solve_transient, 
            solve_heat_transfer, 
            solve_diffusion, 
            solve_laminar_flow, 
            couple_tritium_diffusion_heat_transfer, 
            Time, num_steps, dt, 
            V, V2, Q, D, 
            thermal_conductivity, mu, density, 
            F, Source_c_diffusion, bcs_c, 
            FT, Source_T_diffusion, bcs_T, 
            A1, L1, A2, L2, A3, L3, bcu_transient, bcp_transient, bcs_stationary, 
            ds, dx, volume_marker, n0):
    values_heat_transfers = []
    values_tritium_diffusion = []
    header_heat_transfers = ''
    header_tritium_diffusion = ''
    print('Solving')
    if solve_heat_transfer is True:
        header_heat_transfers = initialise_post_processing(data, "heat_transfers")
    if solve_diffusion is True:
        header_tritium_diffusion = initialise_post_processing(data, "tritium_diffusion")

    if solve_transient is True:
        time_stepping(data=data,
                      solve_heat_transfer=solve_heat_transfer,
                      solve_diffusion=solve_diffusion,
                      solve_laminar_flow=solve_laminar_flow,
                      couple_tritium_diffusion_heat_transfer=couple_tritium_diffusion_heat_transfer,
                      Time=Time,
                      num_steps=num_steps,
                      dt=dt,
                      V=V,
                      U=U,
                      Q=Q,
                      D=D,
                      thermal_conductivity=thermal_conductivity,
                      mu = mu,
                      density = density,
                      F=F,
                      bcs_c=bcs_c,
                      FT=FT,
                      Source_T_diffusion=Source_T_diffusion,
                      Source_c_diffusion=Source_c_diffusion,
                      bcs_T=bcs_T,
                      A1=A1, 
                      L1=L1, 
                      A2=A2, 
                      L2=L2, 
                      A3=A3, 
                      L3=L3, 
                      bcu_transient=bcu_transient, bcp_transient=bcp_transient,
                      ds=ds,
                      dx=dx,
                      volume_marker=volume_marker,
                      header_heat_transfers=header_heat_transfers,
                      header_tritium_diffusion=header_tritium_diffusion,
                      values_heat_transfers=values_heat_transfers,
                      values_tritium_diffusion=values_tritium_diffusion)


    else:
        (v_, q_) = TestFunctions(W)
        w = Function(W)
        (u, p) = split(w)
        ### Steady Part of the Momentum Equation
        def steady(u, dx_fluid):
            T = -p*I + 2*mu*sym(grad(u))
            return (inner(grad(u)*u, v_) + inner(T, grad(v_)) - inner(f, v_)) * dx_fluid
        I = Identity(u.geometric_dimension())
        F = steady(u, dx_fluid) + q_*div(u)*dx_fluid
        J = derivative(F, w)
        problem = NonlinearVariationalProblem(F, w, bcs_stationary, J)
        solver = NonlinearVariationalSolver(problem)
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
        if solve_laminar_flow is True:
            solver.solve()
            (u, p) = w.split()
            output_file << (u, 0.0)
            output_file << (p, 0.0)





if __name__ == "__main__":

    apreprovars = get_apreprovars(2)

    data = get_databases(apreprovars)  # This returns an object data=json.load()
    
    solve_transient, solve_heat_transfer, solve_laminar_flow, \
        solve_diffusion, solve_with_decay, \
        couple_tritium_diffusion_heat_transfer, \
        couple_heat_transfer_laminar_flow, \
        couple_tritium_diffusion_laminar_flow, update_thermal_properties, \
        update_tritium_diffusion_properties, update_fluid_properties = get_solvers(data)  # Gets the solvers

    Time, num_steps, dt = get_solving_parameters(data, solve_transient)  # Gets the parameters (final time, time steps...)

    mesh, n0, volume_marker, dx, surface_marker, ds, mesh_fluid, volume_marker_fluid, dx_fluid, surface_marker_fluid, ds_fluid, n_fluid = define_mesh(data, solve_laminar_flow)

    V, V0, V2, Q, W = define_functionspaces(data, mesh, mesh_fluid)

    c_n, T_n = define_initial_values(solve_heat_transfer, solve_diffusion, data, V)

    Source_c_diffusion, Source_T_diffusion = define_source_terms(solve_heat_transfer, solve_diffusion, dx, data, V)

    bcs_c, Neumann_BC_c_diffusion, Robin_BC_c_diffusion = define_BC_diffusion(data, solve_diffusion, V, surface_marker, ds)

    bcs_T, Neumann_BC_T_diffusion, Robin_BC_T_diffusion = define_BC_heat_transfer(data, solve_heat_transfer, V, surface_marker, ds)

    bcs_stationary, bcu_transient, bcp_transient = define_BC_laminar_flow(data, solve_laminar_flow, V2, Q, W, surface_marker_fluid)

    D, thermal_conductivity, specific_heat, density, mu = define_materials_properties(V0, data, volume_marker, solve_heat_transfer, solve_diffusion, solve_laminar_flow)
    
    if solve_laminar_flow is True:
        ### Unknown and test functions
        (v_, q_) = TestFunctions(W)
        w = Function(W)
        (u_stationary, p_stationary) = split(w)
        u_n = Function(V2)
        u_  = Function(V2)
        p_n = Function(Q)
        p_  = Function(Q)
        u = TrialFunction(V2)
        v = TestFunction(V2)
        p = TrialFunction(Q)
        q = TestFunction(Q)
        U   = 0.5*(u_n + u)
        f   = Constant((0, 0, 0))
        k   = Constant(dt)

    F = define_variational_problem_diffusion(solve_diffusion, solve_transient, dt, solve_with_decay, V, Neumann_BC_c_diffusion, Robin_BC_c_diffusion, Source_c_diffusion)

    FT = define_variational_problem_heat_transfer(solve_heat_transfer, solve_transient, dt, V, specific_heat, density, thermal_conductivity, Neumann_BC_T_diffusion, Robin_BC_T_diffusion, Source_T_diffusion)

    solver, A1, L1, A2, L2, A3, L3 = define_variational_problem_laminar_flow(solve_laminar_flow, solve_transient, V2, Q, W, dx_fluid, ds_fluid, n_fluid, bcu_transient, bcp_transient, bcs_stationary, dt, mu, density)

    solving(data, solve_transient, solve_heat_transfer, solve_diffusion, solve_laminar_flow, couple_tritium_diffusion_heat_transfer, Time, num_steps, dt, V, V2, Q, D, thermal_conductivity, mu, density, F, Source_c_diffusion, bcs_c, FT, Source_T_diffusion, bcs_T, A1, L1, A2, L2, A3, L3, bcu_transient, bcp_transient, bcs_stationary, ds, dx, volume_marker, n0)
