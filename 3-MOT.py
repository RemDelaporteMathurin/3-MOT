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
#os.system('dolfin-convert geo/mesh.inp geo/coucou.xml')

#parser = argparse.ArgumentParser(description='make a mesh')
#parser.add_argument("thickness")
#args = parser.parse_args()

#print(args)

#steel_thickness = float(args.thickness)

#parser = argparse.ArgumentParser()

# parser.add_argument("-m","--mesh", help="XDMF Mesh file input ", required=False, default='mesh_and_markers.xdmf')
# parser.add_argument("-om","--output_mesh", help="XDMF Mesh file output ", required=False, default='mesh_and_markers_from_fenics.xdmf')
# parser.add_argument("-ot","--output_temperature", help="XDMF Mesh file output with temperature", required=False, default='temperature_output.xdmf')
#parser.add_argument("-j","--json_input", help="XDMF Mesh file output with temperature", required=True)

#args = parser.parse_args()


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

def update_bc(t,physic):
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



print('Getting the databases')


#xdmf_encoding = XDMFFile.Encoding.ASCII

#xdmf_encoding = XDMFFile.Encoding.HDF5


materialDB='3-MOT_materials.json'
MOT_parameters='MOT_parameters_RCB.json'
with open(MOT_parameters) as f:
    data = json.load(f)

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

data=byteify(data)

print('Defining the solving parameters')
Time = float(data["solving_parameters"]['final_time'])  #60000*365.25*24*3600.0# final time 
num_steps = data['solving_parameters']['number_of_time_steps'] # number of time steps
dt = Time / num_steps # time step size
t=0 #Initialising time to 0s



print('Defining mesh')
# Read in Mesh and markers from file
mesh = Mesh()
xdmf_in = XDMFFile(mesh.mpi_comm(), str(data['mesh_file']))
xdmf_in.read(mesh)
# prepare output file for writing by writing the mesh to the file
#xdmf_out = XDMFFile(str(data['mesh_file']).split('.')[1]+'_from_fenics.xdmf')

subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
print('Number of cell is '+ str(len(subdomains.array())))


print('Defining Functionspaces')
V = FunctionSpace(mesh, 'P', 1) #FunctionSpace of the solution c
V0 = FunctionSpace(mesh, 'DG', 0) #FunctionSpace of the materials properties

### Define initial values
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


### Boundary Conditions
print('Defining boundary conditions')

# read in Surface markers for boundary conditions
surface_marker_mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
xdmf_in.read(surface_marker_mvc, "surface_marker")

surface_marker_mvc.rename("surface_marker", "surface_marker")
#xdmf_out.write(surface_marker_mvc, xdmf_encoding)

surface_marker = MeshFunction("size_t", mesh, surface_marker_mvc)

ds = Measure('ds', domain=mesh, subdomain_data = surface_marker)

##Tritium Diffusion
if solve_diffusion==True:
  #DC
  bcs_c=list()
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
  Neumann_BC_c_diffusion=[]
  for Neumann in data['physics']['tritium_diffusion']['boundary_conditions']['neumann']:
    value=Neumann['value']
    Neumann_BC_c_diffusion.append([Neumann['value'],Neumann['surface']])
    
    if type(Neumann['surface'])==list:
      for surface in Neumann['surface']:
        Neumann_BC_c_diffusion.append([ds(surface),value])
    else:
      Neumann_BC_c_diffusion.append([ds(Neumann['surface']),value])


#Robins
Robin_BC_c_diffusion=[]
for Robin in data['physics']['tritium_diffusion']['boundary_conditions']['robin']:
  value=Function(V)
  k=3.56e-8
  value=eval(Robin['value'])
  #value=conditional(gt(c_n, 0), k*(c_n)**0.74, Constant(0.0))
  if type(Robin['surface'])==list:
    for surface in Robin['surface']:
      Robin_BC_c_diffusion.append([ds(surface),value])
  else:
    Robin_BC_c_diffusion.append([ds(Robin['surface']),value])



##Temperature
if solve_temperature==True:
  #DC
  bcs_T=list()
  #bci_T=DirichletBC(V,23,surface_marker,22)
  #bcs_T.append(bci_T)
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
  Neumann_BC_T_diffusion=[]
  for Neumann in data['physics']['heat_transfers']['boundary_conditions']['neumann']:
    value=Neumann['value']
    
    if type(Neumann['surface'])==list:
      for surface in Neumann['surface']:
        Neumann_BC_T_diffusion.append([ds(surface),value])
    else:
      Neumann_BC_T_diffusion.append([ds(Neumann['surface']),value])


  #Robins
  Robin_BC_T_diffusion=[]
  for Robin in data['physics']['heat_transfers']['boundary_conditions']['robin']:
    
    if type(Robin['surface'])==list:
      for surface in Robin['surface']:
        Robin_BC_T_diffusion.append([ds(surface),Robin['hc_coeff'],Robin['t_amb']])
    else:
      Robin_BC_T_diffusion.append([ds(Robin['surface']),Robin['hc_coeff'],Robin['t_amb']])
  
  #print(Neumann_BC_T_diffusion)


#read in the volume markers
volume_marker_mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
xdmf_in.read(volume_marker_mvc, "volume_marker_material")

#volume_marker_mvc.rename("volume_marker_material", "volume_marker_material")
#xdmf_out.write(volume_marker_mvc, xdmf_encoding)

volume_marker = MeshFunction("size_t", mesh, volume_marker_mvc)


###Defining materials properties
print('Defining the materials properties')

D  = Function(V0) #Diffusion coefficient

def calculate_D(T,material_id):
  R=8.314 #Perfect gas constant
  if material_id==1: #Concrete
    return 2e-6#7.3e-7*np.exp(-6.3e3/T)
  elif material_id==2: #Polymer
    return 2.0e-7*np.exp(-29000.0/R/T)
  elif material_id==3: #steel
    return 7.3e-7*np.exp(-6.3e3/T)#1e-16#2e-6

if solve_with_decay==True:
  decay=np.log(2)/(12.33*365.25*24*3600) #Tritium Decay constant [s-1]
else:
  decay=0

thermal_conductivity=Function(V0)
thermal_conductivity_values=[150.0,122.0,16.12,29.0]


##Assigning each to each cell its properties
for cell_no in range(len(volume_marker.array())):
  volume_id=volume_marker.array()[cell_no] #This is the volume id
  #print(str(volume_id))
  material_id=volume_id
  thermal_conductivity.vector()[cell_no]=thermal_conductivity_values[material_id]

  D.vector()[cell_no]=calculate_D(data['physics']['heat_transfers']['initial_value'],material_id)
  
  #print(D.vector()[cell_no])


### Define variational problem
print('Defining the variational problem')


if solve_diffusion==True:
  c = TrialFunction(V)#c is the tritium concentration
  vc = TestFunction(V)
  f = Expression(str(data['physics']['tritium_diffusion']['source_term']),t=0,degree=2)#This is the tritium volumetric source term 
  F=((c-c_n)/dt)*vc*dx + D*dot(grad(c), grad(vc))*dx + (-f+decay*c)*vc*dx 
  for Neumann in Neumann_BC_c_diffusion:
    F += vT * Neumann[1]*Neumann[0] 
  for Robin in Robin_BC_c_diffusion:
    F += D*vc*Robin[1]*Robin[0]
  ac,Lc= lhs(F),rhs(F)


if solve_temperature==True:
  T = TrialFunction(V) #T is the temperature
  vT = TestFunction(V)
  q = Expression(str(data['physics']['heat_transfers']['source_term']),t=0,degree=2) #q is the volumetric heat source term
  
  FT = 1e6*((T-T_n)/dt)*vT*dx +thermal_conductivity*dot(grad(T), grad(vT))*dx - q*vT*dx #This is the heat transfer equation     

  for Neumann in Neumann_BC_T_diffusion:
    #print(Neumann)
    FT += -1/thermal_conductivity* vT * Neumann[1]*Neumann[0] 

  for Robin in Robin_BC_T_diffusion:
    FT += 1/thermal_conductivity *vT* Robin[1] * (T-Robin[2])*Robin[0]


  aT, LT = lhs(FT), rhs(FT) #Rearranging the equation

### Time-stepping
T = Function(V)
c = Function(V)
off_gassing=list()
output_file  = File("Solution.pvd")

if calculate_off_gassing==True:
  file_off_gassing = "off_gassing.csv"

for n in range(num_steps):

  
  # Update current time
  print("t= "+str(t)+" s")
  print("t= "+str(t/3600/24/365.25)+" years")
  print(str(100*t/Time)+" %")
  t += dt
  

  # Compute solution concentration
  if solve_diffusion==True:
    f.t += dt
    solve(ac==Lc,c,bcs_c)
    output_file << (c,t)
    c_n.assign(c)
    bcs_c=update_bc(t,"tritium_diffusion")
  # Compute solution temperature
  if solve_temperature==True:
    q.t += dt
    solve(aT == LT, T, bcs_T)
    output_file << (T,t)
    T_n.assign(T)
    bcs_T=update_bc(t,"heat_transfers")
  #Update the materials properties
  if solve_diffusion_coefficient_temperature_dependent==True and solve_temperature==True and solve_diffusion==True:
    for cell in cells(mesh):
      cell_no=cell.index()
      material_id=volume_marker.array()[cell_no]
      Ta=0
      for i in range(0,12,3):
        Ta+=T(cell.get_vertex_coordinates()[i],cell.get_vertex_coordinates()[i+1],cell.get_vertex_coordinates()[i+2])
      Ta=Ta/4
      D_value = calculate_D(Ta,material_id)
      D.vector()[cell_no] = D_value #Assigning for each cell the corresponding diffusion coeff
  
  #Calculate off-gassing
  if calculate_off_gassing==True:
    #g=conditional(gt(c_n, 0), k*(c_n)**0.74, Constant(0.0))
    g=Robin_BC_c_diffusion[0][1]
    #off_gassing_per_day=3600*24*(assemble(g*ds(1))+assemble(g*ds(2))+assemble(g*ds(3))+assemble(g*ds(4))+assemble(g*ds(5))+assemble(g*ds(6)))/dt #off-gassing in mol/day
    off_gassing_per_day=assemble(g*ds)/dt*24*3600
    #off_gassing+=assemble(g*ds(2))*24*3600 +assemble(g*ds(3))*24*3600 +assemble(g*ds(4)*24*3600)  +assemble(g*ds(5)*24*3600)  +assemble(g*ds(6)*24*3600 )
    i=0
    print(off_gassing_per_day)
    off_gassing.append([off_gassing_per_day,t/3600/24/365])

    with open(file_off_gassing, "w") as output:
      writer = csv.writer(output, lineterminator='\n')
      writer.writerow('ct')
      for val in off_gassing:
        writer.writerows([val])