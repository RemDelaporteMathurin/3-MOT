#from fenics import *
from dolfin import *
import numpy as np
import csv
import sys
import os
import argparse
import json
import ast
from pprint import pprint

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



print('Getting the databases')

xdmf_encoding = XDMFFile.Encoding.ASCII
#xdmf_encoding = XDMFFile.Encoding.HDF5


materialDB='3-MOT_materials.json'
MOT_parameters='MOT_parameters_RCB.json'
with open(MOT_parameters) as f:
    data = json.load(f)

print('Getting the solvers')
solve_temperature=True
solve_diffusion=True
solve_diffusion_coefficient_temperature_dependent=True
solve_with_decay=False
calculate_off_gassing=False

#data=byteify(data)
pprint(data)

print('Defining the solving parameters')
Time = data["solving_parameters"]['final_time']  #60000*365.25*24*3600.0# final time 
num_steps = data['solving_parameters']['number_of_time_steps'] # number of time steps
dt = Time / num_steps # time step size
t=0 #Initialising time to 0s

# prepare output file for writing by writing the mesh to the file
<<<<<<< HEAD
xdmf_out = XDMFFile(MPI.comm_world, data['mesh_file'].split('.')[1]+'_from_fenics.xdmf')
#xdmf_out.write(mesh, xdmf_encoding)
=======

xdmf_out = XDMFFile(MPI.comm_world, data['mesh_file'].split('.')[1]+'_from_fenics.xdmf')
xdmf_out.write(mesh, xdmf_encoding)
>>>>>>> 048ef430e0e8f24a3cdd91de2773d4fd70c7dbce


cells=2 

print('Defining mesh')
#Create mesh and define function space

# Read in Mesh and markers from file
mesh = Mesh()
xdmf_in = XDMFFile(MPI.comm_world, data['mesh_file'])
xdmf_in.read(mesh)

subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
print('Number of cell is '+ str(len(subdomains.array())))


print('Defining Functionspaces')
V = FunctionSpace(mesh, 'P', 1) #FunctionSpace of the solution c
V0 = FunctionSpace(mesh, 'DG', 0) #FunctionSpace of the materials properties

### Define initial values
print('Defining initial values')
##Tritium concentration
iniC = Expression(str(data['physics']['tritium_diffusion']['initial_value']),degree=1) 
c_n = interpolate(iniC, V)
##Temperature
iniT = Expression(str(data['physics']['heat_Transfers']['initial_value']),degree=1)
T_n = interpolate(iniT, V)


### Boundary Conditions
print('Defining boundary conditions')

# read in Surface markers for boundary conditions
surface_marker_mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
xdmf_in.read(surface_marker_mvc, "surface_marker")

surface_marker_mvc.rename("surface_marker", "surface_marker")
xdmf_out.write(surface_marker_mvc, xdmf_encoding)

surface_marker = MeshFunction("size_t", mesh, surface_marker_mvc)

ds = Measure('ds', domain=mesh, subdomain_data = surface_marker)

##Tritium Diffusion
if solve_diffusion==True:
  #DC
  bcs_c=list()
  for BC in data['physics']['tritium_diffusion']['boundary_conditions']['dc']:
    value_BC=BC['value'] #todo make this value able to be an Expression (time or space dependent)

    #bci_c=DirichletBC(V,value_BC,surface_marker,BC['surface'])
    #bcs_c.append(bci_c)


  #Neumann
  Neumann_BC_c_diffusion=[]
  for Neumann in data['physics']['tritium_diffusion']['boundary_conditions']['neumann']:
    value=Neumann['value']
    Neumann_BC_c_diffusion.append(value*ds(Neumann['surface']))


  #Robins
  Robin_BC_c_diffusion=[]
  for Robin in data['physics']['tritium_diffusion']['boundary_conditions']['robin']:
    value=Function(V)
    value=ast.literal_eval(Robin['expression'])
    Robin_BC_c_diffusion.append(value*ds(Robin))



##Temperature


if solve_temperature==True:
  #DC
  bcs_T=list()
  for BC in data['physics']['heat_Transfers']['boundary_conditions']['dc']:
    value_BC=BC['value'] #todo make this value able to be an Expression (time or space dependent)
    
    bci_T=DirichletBC(V,value_BC,surface_marker,BC['surface'])
    bcs_T.append(bci_T)


  #Neumann
  Neumann_BC_T_diffusion=[]
  for Neumann in data['physics']['heat_Transfers']['boundary_conditions']['neumann']:
    value=Neumann['value']
    Neumann_BC_T_diffusion.append(value*ds(Neumann['surface']))


  #Robins
  Robin_BC_T_diffusion=[]
  for Robin in data['physics']['heat_Transfers']['boundary_conditions']['robin']:
    Robin_BC_T_diffusion.append([ds(Robin['surface']),Robin['hc_coeff'],Robin['t_amb']])


#read in the volume markers
volume_marker_mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
xdmf_in.read(volume_marker_mvc, "volume_marker_material")

volume_marker_mvc.rename("volume_marker_material", "volume_marker_material")
xdmf_out.write(volume_marker_mvc, xdmf_encoding)

volume_marker = MeshFunction("size_t", mesh, surface_marker_mvc)


###Defining materials properties
print('Defining the materials properties')

D  = Function(V0) #Diffusion coefficient
def calculate_D(T,material_id):
  R=8.314 #Perfect gas constant
  if material_id==0: #Steel
    return 7.3e-7*np.exp(-6.3e3/T)
  elif material_id==1: #Polymer
    return 2.0e-7*np.exp(-29000.0/R/T)
  elif material_id==2: #Concrete
    return 2e-6

if solve_with_decay==True:
  decay=np.log(2)/(12.33*365.25*24*3600) #Tritium Decay constant [s-1]
else:
  decay=0

thermal_diffusivity=Function(V0)
thermal_diffusivity_values=[4e-6,0.15e-6,0.54e-6]


##Assigning each to each cell its properties
for cell_no in range(len(volume_marker.array())):
  material_id=volume_marker.array()[cell_no]
  
  thermal_diffusivity.vector()[cell_no]=thermal_diffusivity_values[subdomain_no]

  D.vector()[cell_no]=calculate_D(data['physics']['heat_Transfers']['initial_value'],material_id)
  #print(D.vector()[cell_no])
  

### Define variational problem
print('Defining the variational problem')


if solve_diffusion==True:
  c = TrialFunction(V)#c is the tritium concentration
  vc = TestFunction(V)
  f = Expression(str(data['physics']['tritium_diffusion']['source_term']),t=0,degree=2)#This is the tritium volumetric source term 
  F=((c-c_n)/dt)*vc*dx + D*dot(grad(c), grad(vc))*dx + (f+decay*c)*vc*dx 
  for Neumann in Neumann_BC_c_diffusion:
    F += D*Neumann*vc
  for Robin in Robin_BC_c_diffusion:
    F += D*Robin*vc
  ac,Lc= lhs(F),rhs(F)


if solve_temperature==True:
  T = TrialFunction(V) #T is the temperature
  vT = TestFunction(V)
  q = Expression(str(data['physics']['heat_transfers']['source_term']),t=0,degree=2) #q is the volumetric heat source term
  
  FT = ((T-T_n)/dt)*vT*dx + thermal_diffusivity*dot(grad(T), grad(vT))*dx + q*vT*dx #This is the heat transfer equation

  for Neumann in Neumann_BC_T_diffusion:
    FT += thermal_diffusivity * Neumann * vT

  for Robin in Robin_BC_T_diffusion:
    FT += thermal_diffusivity * Robin(1) * (T_n-Robin(2))*Robin(0)


  aT, LT = lhs(FT), rhs(FT) #Rearranging the equation



### Time-stepping
T = Function(V)
c=Function(V)
off_gassing=list()
output_file = 'T.pvd'

if calculate_off_gassing==True:
  file_off_gassing = "Solutions/Off-gassing/off_gassing.csv"

for n in range(num_steps):

  
  # Update current time
  print("t= "+str(t)+" s")
  print("t= "+str(t/3600/24/365.25)+" years")
  print(str(100*t/Time)+" %")
  t += dt
  f.t += dt


  # Compute solution concentration
  if solve_diffusion==True:
    solve(ac==Lc,c,bcs_c)
    output_file << (c,t)
  # Compute solution temperature
  if solve_temperature==True:
    solve(aT == LT, T, bcs_T)
    output_file << (T,t)

  # Update previous solution
  c_n.assign(c)
  T_n.assign(T)

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
    off_gassing_per_day=assemble(g*ds(1,domain=mesh))*24*3600 #off-gassing in mol/day
    i=0
    print(off_gassing_per_day)
    while i<int(dt/24/3600):
      off_gassing.append(off_gassing_per_day)
      i+=1
    with open(file_off_gassing, "w") as output:
      writer = csv.writer(output, lineterminator='\n')
      writer.writerow('c')
      for val in off_gassing:
        writer.writerow([val])