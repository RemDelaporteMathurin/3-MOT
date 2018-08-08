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
    #return 'MOT_parameters_breeder_blankets.json'

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
    data=byteify(data)
    return data

def get_solvers(data):
    print('Getting the solvers')
    if data['physics']['solve_heat_transfer']==1:
        solve_heat_transfer=True
    else:
        solve_heat_transfer=False
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
    return solve_heat_transfer,solve_diffusion,solve_diffusion_coefficient_temperature_dependent,solve_with_decay

def get_solving_parameters(data):
    print('Getting the solving parameters')
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
    return mesh,xdmf_in

def define_functionspaces(data):
    print('Defining Functionspaces')
    V = FunctionSpace(mesh, 'P', 1) #FunctionSpace of the solution c
    V0 = FunctionSpace(mesh, 'DG', 0) #FunctionSpace of the materials properties
    return V,V0

def define_initial_values(solve_heat_transfer,solve_diffusion,data,V):
    ##Tritium concentration
    c_n=Function(V)
    T_n=Function(V)
    if solve_diffusion==True:
      print('Defining initial values tritium diffusion')
      #print(str(data['physics']['tritium_diffusion']['initial_value']))
      iniC = Expression(str(data['physics']['tritium_diffusion']['initial_value']),degree=2) 
      c_n = interpolate(iniC, V)
    ##Temperature
    if solve_heat_transfer==True:
      print('Defining initial values heat transfer')
      #print(str(data['physics']['heat_transfers']['initial_value']))
      iniT = Expression(str(data['physics']['heat_transfers']['initial_value']),degree=2) 
      T_n = interpolate(iniT, V)
    return c_n,T_n

def get_surface_marker(mesh,xdmf_in):
    print('Marking the surfaces')
    # read in Surface markers for boundary conditions
    surface_marker_mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
    xdmf_in.read(surface_marker_mvc, "surface_marker")
    surface_marker_mvc.rename("surface_marker", "surface_marker")
    #xdmf_out.write(surface_marker_mvc, xdmf_encoding)
    surface_marker = MeshFunction("size_t", mesh, surface_marker_mvc)

    ds = Measure('ds', domain=mesh, subdomain_data = surface_marker)
    return surface_marker,ds

def define_BC_diffusion(data,solve_diffusion,V,surface_marker,ds):
    ##Tritium Diffusion
    bcs_c=list()
    Neumann_BC_c_diffusion=[]
    Robin_BC_c_diffusion=[]
    if solve_diffusion==True:
      print('Defining BC tritium diffusion')
      #DC
      for DC in data['physics']['tritium_diffusion']['boundary_conditions']['dc']:
        value_DC=Expression(DC['value'],t=0,degree=2)
        if type(DC['surface'])==list:
          for surface in DC['surface']:
            #print(surface)
            bci_c=DirichletBC(V,value_DC,surface_marker,surface)
            bcs_c.append(bci_c)
            #print(bci_T)
        else:
          #print(DC)
          bci_c=DirichletBC(V,value_DC,surface_marker,DC['surface'])
          bcs_c.append(bci_c)
          #print(bci_T)

      #Neumann
      for Neumann in data['physics']['tritium_diffusion']['boundary_conditions']['neumann']:
        value=Expression(Neumann['value'],t=0,degree=2)
        if type(Neumann['surface'])==list:
          for surface in Neumann['surface']:
            Neumann_BC_c_diffusion.append([ds(surface),value])
        else:
          Neumann_BC_c_diffusion.append([ds(Neumann['surface']),value])
      
      #Robins
      for Robin in data['physics']['tritium_diffusion']['boundary_conditions']['robin']:
        value=Function(V)
        value=eval(Robin['value'])
        #value=conditional(gt(c_n, 0), k*(c_n)**0.74, Constant(0.0))
        if type(Robin['surface'])==list:
          for surface in Robin['surface']:
            Robin_BC_c_diffusion.append([ds(surface),value])
        else:
          Robin_BC_c_diffusion.append([ds(Robin['surface']),value])
    return bcs_c,Neumann_BC_c_diffusion,Robin_BC_c_diffusion

def define_BC_heat_transfer(data,solve_heat_transfer,V,surface_marker,ds):
    ##Temperature
    bcs_T=list()
    Neumann_BC_T_diffusion=[]
    Robin_BC_T_diffusion=[]

    if solve_heat_transfer==True:
      print('Defining BC heat transfer')
      #DC
      for DC in data['physics']['heat_transfers']['boundary_conditions']['dc']:
        #value_DC=DC['value'] #todo make this value able to be an Expression (time or space dependent)
        value_DC=Expression(DC['value'],t=0,degree=2)
        if type(DC['surface'])==list:
          for surface in DC['surface']:
            #print(surface)
            bci_T=DirichletBC(V,value_DC,surface_marker,surface)
            bcs_T.append(bci_T)
            #print(bci_T)
        else:
          #print(DC)
          bci_T=DirichletBC(V,value_DC,surface_marker,DC['surface'])
          bcs_T.append(bci_T)
          #print(bci_T)

      #Neumann
      for Neumann in data['physics']['heat_transfers']['boundary_conditions']['neumann']:
        value=Expression(Neumann['value'],t=0,degree=2)

        if type(Neumann['surface'])==list:
          for surface in Neumann['surface']:
            Neumann_BC_T_diffusion.append([ds(surface),value])
        else:
          Neumann_BC_T_diffusion.append([ds(Neumann['surface']),value])

      #Robins
      for Robin in data['physics']['heat_transfers']['boundary_conditions']['robin']:
      
        if type(Robin['surface'])==list:
          for surface in Robin['surface']:
            Robin_BC_T_diffusion.append([ds(surface),Robin['hc_coeff'],Robin['t_amb']])
        else:
          Robin_BC_T_diffusion.append([ds(Robin['surface']),Robin['hc_coeff'],Robin['t_amb']])
    
      #print(Neumann_BC_T_diffusion)
    return bcs_T,Neumann_BC_T_diffusion,Robin_BC_T_diffusion

def get_volume_markers(mesh):
    print('Marking the volumes')
    #read in the volume markers
    volume_marker_mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim())
    xdmf_in.read(volume_marker_mvc, "volume_marker_material")

    #volume_marker_mvc.rename("volume_marker_material", "volume_marker_material")
    #xdmf_out.write(volume_marker_mvc, xdmf_encoding)

    volume_marker = MeshFunction("size_t", mesh, volume_marker_mvc)
    dx = Measure('dx', domain=mesh, subdomain_data=volume_marker)
    return volume_marker,dx

def calculate_D(T,material_id):
  R=8.314 #Perfect gas constant
  if material_id=="concrete": #Concrete
    return 2e-6#7.3e-7*np.exp(-6.3e3/T)
  elif material_id=="polymer": #Polymer
    return 2.0e-7*np.exp(-29000.0/R/T)
  elif material_id=="steel": #steel
    return 7.3e-7*np.exp(-6.3e3/T)#1e-16#2e-6
  else:
    raise ValueError("!!ERROR!! Unable to find "+str(material_id)+" as material ID in the database "+str(inspect.stack()[0][3]))

def calculate_thermal_conductivity(T,material_id):
  R=8.314 #Perfect gas constant
  if material_id=="concrete":
    return 0.5
  elif material_id=="tungsten":
    return 150
  elif material_id=="lithium_lead":
    return 50
  elif material_id=="eurofer":
    return 29
  else:
    raise ValueError("!!ERROR!! Unable to find "+str(material_id)+" as material ID in the database "+str(inspect.stack()[0][3]))

def calculate_specific_heat(T,material_id):
  R=8.314 #Perfect gas constant
  if material_id=="concrete":
    return 880
  elif material_id=="tungsten":
    return 130
  elif material_id=="lithium_lead":
    return 500
  elif material_id=="eurofer":
    return 675
  else:
    raise ValueError("!!ERROR!! Unable to find "+str(material_id)+" as material ID in the database "+str(inspect.stack()[0][3]))

def calculate_density(T,material_id):

  R=8.314 #Perfect gas constant
  if material_id=="concrete":
    return 2400
  elif material_id=="tungsten":
    return 19600
  elif material_id=="lithium_lead":
    return 11600
  elif material_id=="eurofer":
    return 7625
  else:
    raise ValueError("!!ERROR!! Unable to find "+str(material_id)+" as material ID in the database "+str(inspect.stack()[0][3]))

def which_material_is_it(volume_id,data):
  for material in data["structure_and_materials"]["materials"]:
    for volumes in material["volumes"]:
      if volume_id in [volumes]:
        #print('Coucou toi')
        material_id=material["material"]
        break
  return material_id

def define_materials_properties(V0,data,volume_marker):
    ###Defining materials properties
    print('Defining the materials properties')
    
    D  = Function(V0) #Diffusion coefficient
    thermal_conductivity=Function(V0)
    specific_heat=Function(V0)
    density=Function(V0)
    ##Assigning each to each cell its properties
    for cell_no in range(len(volume_marker.array())):
      volume_id=volume_marker.array()[cell_no] #This is the volume id (Trellis)
      material_id=which_material_is_it(volume_id,data)
      if solve_heat_transfer==True:
        thermal_conductivity.vector()[cell_no]=calculate_thermal_conductivity(data['physics']['heat_transfers']['initial_value'],material_id)
        density.vector()[cell_no]=calculate_density(data['physics']['heat_transfers']['initial_value'],material_id)
        specific_heat.vector()[cell_no]=calculate_specific_heat(data['physics']['heat_transfers']['initial_value'],material_id)
      if solve_diffusion==True:
        D.vector()[cell_no]=calculate_D(data['physics']['heat_transfers']['initial_value'],material_id)
    
      #print(D.vector()[cell_no])
    return D,thermal_conductivity,specific_heat,density

def define_variational_problem_diffusion(solve_diffusion,solve_with_decay,V,data):
    if solve_diffusion==True:
        print('Defining variation problem tritium diffusion')
        c = TrialFunction(V)#c is the tritium concentration
        vc = TestFunction(V)
        if solve_with_decay==True:
            decay=np.log(2)/(12.33*365.25*24*3600) #Tritium Decay constant [s-1]
        else:
            decay=0
    
        f = Expression(str(data['physics']['tritium_diffusion']['source_term']),t=0,degree=2)#This is the tritium volumetric source term 
        F=((c-c_n)/dt)*vc*dx + D*dot(grad(c), grad(vc))*dx + (-f+decay*c)*vc*dx 
        for Neumann in Neumann_BC_c_diffusion:
            F += vc * Neumann[1]*Neumann[0] 
        for Robin in Robin_BC_c_diffusion:
            F += D*vc*Robin[1]*Robin[0]
        return F,f
    return False,False
    
def define_variational_problem_heat_transfer(solve_heat_transfer,V,data):
    
    if solve_heat_transfer==True:
        print('Defining variation problem heat transfer')
        T = TrialFunction(V) #T is the temperature
        vT = TestFunction(V)
        q = Expression(str(data['physics']['heat_transfers']['source_term']),t=0,degree=2) #q is the volumetric heat source term

        FT = specific_heat*density*((T-T_n)/dt)*vT*dx +thermal_conductivity*dot(grad(T), grad(vT))*dx - q*vT*dx #This is the heat transfer equation         
        for Neumann in Neumann_BC_T_diffusion:
            #print(Neumann)
            FT += - vT * Neumann[1]*Neumann[0]    
        for Robin in Robin_BC_T_diffusion:
            FT += vT* Robin[1] * (T-Robin[2])*Robin[0]
        return FT,q
    return False,False

def update_D(mesh,volume_marker,D,T):
    for cell in cells(mesh):
      cell_no=cell.index()
      material_id=volume_marker.array()[cell_no]
      Ta=0
      for i in range(0,12,3):
        Ta+=T(cell.get_vertex_coordinates()[i],cell.get_vertex_coordinates()[i+1],cell.get_vertex_coordinates()[i+2])
      Ta=Ta/4
      D_value = calculate_D(Ta,material_id)
      D.vector()[cell_no] = D_value #Assigning for each cell the corresponding diffusion coeff

def update_bc(t,physic):
  if physic=="tritium_diffusion":
    for Neumann in Neumann_BC_c_diffusion:
      Neumann[1].t=t
  if physic=="heat_transfers":
    for Neumann in Neumann_BC_T_diffusion:
      Neumann[1].t=t


  bcs=list()
  for DC in data['physics'][physic]['boundary_conditions']['dc']:
    #value_DC=DC['value'] #todo make this value able to be an Expression (time or space dependent)
    value_DC=Expression(DC['value'],t=t,degree=2)
    if type(DC['surface'])==list:
      for surface in DC['surface']:
        #print(surface)
        bci=DirichletBC(V,value_DC,surface_marker,surface)
        bcs.append(bci)
        #print(bci_T)
    else:
      print(DC)
      bci=DirichletBC(V,value_DC,surface_marker,DC['surface'])
      bcs.append(bci)
      #print(bci_T)
  return bcs

def time_stepping(solve_heat_transfer,solve_diffusion,solve_diffusion_coefficient_temperature_dependent,Time,num_steps,dt,F,f,bcs_c,FT,q,bcs_T,ds):
    ### Time-stepping
    print('Time stepping')
    T = Function(V)
    c = Function(V)
    off_gassing=list()
    output_file  = File(data["output_file"])
    t=0

    n0 = FacetNormal(mesh)
    for n in range(num_steps):

    
      # Update current time
      print("t= "+str(t)+" s")
      print("t= "+str(t/3600/24/365.25)+" years")
      print(str(100*t/Time)+" %")
      t += dt
    

      # Compute solution concentration
      if solve_diffusion==True:
        f.t += dt
        solve(lhs(F)==rhs(F),c,bcs_c)
        output_file << (c,t)
        c_n.assign(c)
        bcs_c=update_bc(t,"tritium_diffusion")
      # Compute solution temperature
      if solve_heat_transfer==True:
        q.t += dt
        solve(lhs(FT)==rhs(FT), T, bcs_T)
        output_file << (T,t)
        T_n.assign(T)
        bcs_T=update_bc(t,"heat_transfers")
      #Update the materials properties
      if solve_diffusion_coefficient_temperature_dependent==True and solve_heat_transfer==True and solve_diffusion==True:
        D=update_D(mesh,volume_marker,D,T)
      
      flux_1 =-assemble(dot(2e-6*grad(c), n0)*ds(1))
      
      #assemble(conditional(gt(c_n, 0), 5.08e-6*(c_n)**0.74, Constant(0.0))*ds(1))#assemble(c*ds(1))#+assemble(dot(grad(c), n0)*ds(2))+assemble(dot(grad(c), n0)*ds(3))+assemble(dot(grad(c), n0)*ds(4))+assemble(dot(grad(c), n0)*ds(5))+assemble(dot(grad(c), n0)*ds(6))
      off_gassing.append([flux_1,t/365/24/3600])
      with open("OG20degC.csv", "w") as output:
        writer = csv.writer(output, lineterminator='\n')
        writer.writerow('ct')
        for val in off_gassing:
          writer.writerows([val])
    
    return

if __name__=="__main__":
    apreprovars=get_apreprovars(2)
    data=get_databases(apreprovars) #This returns an object data=json.load()
    solve_heat_transfer,solve_diffusion,solve_diffusion_coefficient_temperature_dependent,solve_with_decay=get_solvers(data) #Gets the solvers

    Time, num_steps,dt=get_solving_parameters(data) #Gets the parameters (final time, time steps...)

    mesh, xdmf_in=define_mesh(data)

    V, V0=define_functionspaces(data)

    c_n, T_n=define_initial_values(solve_heat_transfer,solve_diffusion,data,V)

    surface_marker, ds=get_surface_marker(mesh,xdmf_in)

    bcs_c, Neumann_BC_c_diffusion,Robin_BC_c_diffusion=define_BC_diffusion(data,solve_diffusion,V,surface_marker,ds)

    bcs_T, Neumann_BC_T_diffusion,Robin_BC_T_diffusion=define_BC_heat_transfer(data,solve_heat_transfer,V,surface_marker,ds)

    volume_marker, dx=get_volume_markers(mesh)

    D,thermal_conductivity,specific_heat,density=define_materials_properties(V0,data,volume_marker)

    F,f=define_variational_problem_diffusion(solve_diffusion,solve_with_decay,V,data)

    FT,q=define_variational_problem_heat_transfer(solve_heat_transfer,V,data)

    time_stepping(solve_heat_transfer,solve_diffusion,solve_diffusion_coefficient_temperature_dependent,Time,num_steps,dt,F,f,bcs_c,FT,q,bcs_T,ds)
