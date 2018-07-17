from fenics import *
from dolfin import *
import numpy as np
import csv
import sys
import os
import argparse
import json


#os.system('dolfin-convert geo/mesh.inp geo/coucou.xml')


#parser = argparse.ArgumentParser(description='make a mesh')
#parser.add_argument("thickness")
#args = parser.parse_args()

#print(args)

#steel_thickness = float(args.thickness)

print('Getting the databases')
materialDB='3-MOT_materials.json'
_3-MOT_parameters='3-MOT_parameters.json'

print('Getting the solvers')
solve_temperature=False
solve_diffusion=True
solve_diffusion_coefficient_temperature_dependent=True
solve_with_decay=False
calculate_off_gassing=False

print('Defining the solving parameters')
Time =5e10#60000*365.25*24*3600.0# final time
num_steps = 10 # number of time steps
dt = Time / num_steps # time step size
t=0 #Initialising time to 0s


cells=2

print('Defining mesh')
#Create mesh and define function space
#mesh=Mesh('geo/mesh_rcb_12_nodes_in_steel.xml')
mesh=Mesh('geo/mesh_rcb_'+str(cells)+'_nodes_in_steel.xml')
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
print('Number of cell is '+ str(len(subdomains.array())))
#WARNING : MAKE SURE THESE DIMENSIONS ARE THE SAME AS IN THE CUBIT SCRIPT
concrete_thickness=240e-3
steel_thickness=2e-3
polymer_thickness=20e-3
internal_cavity_dimension=1738e-3

print('Defining Functionspaces')
V = FunctionSpace(mesh, 'P', 1) #FunctionSpace of the solution c
V0 = FunctionSpace(mesh, 'DG', 0) #FunctionSpace of the materials properties

### Define initial values
print('Defining initial values')
##Tritium concentration
iniC = Expression('0',degree=1)
c_n = interpolate(iniC, V)
##Temperature
initial_temperature=273.15+14
iniT = Expression(str(initial_temperature),degree=1)
T_n = interpolate(iniT, V)


### Boundary Conditions
print('Defining boundary conditions')
boundary_parts=MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
inside =  CompiledSubDomain("(fabs(x[0])<=side && fabs(x[1])<side && fabs(x[2])<side) && on_boundary", side = (internal_cavity_dimension+steel_thickness)/2)
outside = CompiledSubDomain("(fabs(x[0])>=side || fabs(x[1])>=side || fabs(x[2])>=side) && on_boundary", side = (internal_cavity_dimension+concrete_thickness+polymer_thickness+steel_thickness)/2)
#Marking the Boundaries 
markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1,0) #initialization to 0
ds=ds(subdomain_data=markers) #Defining the measurements ds
inside.mark(markers, 0)
outside.mark(markers,1)

##Tritium concentration
inside_bc_c=Expression('1', t=0, degree=1) #0.4*exp(-t*log(2)/(12.3*365.25*24*3600))
outside_bc_c=Expression('0', t=0, degree=2)
bci_c=DirichletBC(V,inside_bc_c,inside)
bco_c=DirichletBC(V,outside_bc_c,outside)
#g = Function(V)
k=3.56e-8
g = Constant(0.0)
#g=conditional(gt(c_n, 0), k*(c_n)**0.74, Constant(0.0))#
bcs_c=list()
bcs_c.append(bci_c)
#bcs_c.append(bco_c)

##Temperature
inside_bc_T=0
outside_bc_T=Expression('14+273.15+7*cos(2*3.14*t/365.25/24/3600)+16*cos(2*3.14*t/24/3600)', t=0, degree=2)
bci_T=DirichletBC(V,inside_bc_T,inside)
bco_T=DirichletBC(V,outside_bc_T,outside)
bcs_T=list()
bcs_T.append(bco_T)


###Defining the subdomains
print('Defining the subdomains')
class steel(SubDomain):
  def inside(self, x, on_boundary):
      return True if abs(x[0])<=(internal_cavity_dimension+steel_thickness)/2 and abs(x[1])<=(internal_cavity_dimension+steel_thickness)/2 and abs(x[2])<=(internal_cavity_dimension+steel_thickness)/2 else False

class polymer(SubDomain):
  def inside(self, x, on_boundary):
      return True if abs(x[0]) <= (steel_thickness+internal_cavity_dimension+polymer_thickness)/2 and abs(x[1]) <= (steel_thickness+internal_cavity_dimension+polymer_thickness)/2 and abs(x[2]) <= (steel_thickness+internal_cavity_dimension+polymer_thickness)/2 else False

class concrete(SubDomain):
  def inside(self, x, on_boundary):
      return True if abs(x[0]) <= (steel_thickness+internal_cavity_dimension+polymer_thickness+concrete_thickness)/2 and abs(x[1]) <= (steel_thickness+internal_cavity_dimension+polymer_thickness+concrete_thickness)/2 and abs(x[2]) <= (steel_thickness+internal_cavity_dimension+polymer_thickness+concrete_thickness)/2 else False

##marking subdomains
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
subdomains.set_all(4)
concrete = concrete()
concrete.mark(subdomains, 2)
polymer = polymer()
polymer.mark(subdomains, 1)
steel = steel()
steel.mark(subdomains, 0)



###Defining materials properties
print('Defining the materials properties')
density=Function(V0) #density [kg/m3]
density_values=[7850,1000,2000] #steel,polymer,concrete

D  = Function(V0) #Diffusion coefficient
def calculate_D(T,subdomain_no):
  R=8.314 #Perfect gas constant
  if subdomain_no==0: #Steel
    return 1e-17#7.3e-7*np.exp(-6.3e3/T)
  elif subdomain_no==1: #Polymer
    return 1e-7#2.0e-7*np.exp(-29000.0/R/T)
  elif subdomain_no==2: #Concrete
    return 1e-4#2e-6

if solve_with_decay==True:
  decay=np.log(2)/(12.33*365.25*24*3600) #Tritium Decay constant [s-1]
else:
  decay=0

thermal_diffusivity=Function(V0)
thermal_diffusivity_values=[4e-6,0.15e-6,0.54e-6]


##Assigning each to each cell its properties
for cell_no in range(len(subdomains.array())):
  subdomain_no = subdomains.array()[cell_no]
  if subdomain_no==4:
    print("WARNING : SOME CELLS IN SUBDOMAINS AREN'T DETECTED")
  thermal_diffusivity.vector()[cell_no]=thermal_diffusivity_values[subdomain_no]
  D.vector()[cell_no]=calculate_D(initial_temperature,subdomain_no)
  #print(D.vector()[cell_no])
  
print(str(cell_no)+' Cells found')

### Define variational problem
print('Defining the variational problem')
c = TrialFunction(V)#c is the tritium concentration
#c=Function(V)
vc = TestFunction(V)
f = Expression('0',t=0,degree=2)#This is the tritium volumetric source term 
F=((c-c_n)/dt)*vc*dx + D*dot(grad(c), grad(vc))*dx + (f+decay*c)*vc*dx +D*g*vc*ds
ac,Lc= lhs(F),rhs(F)



T = TrialFunction(V) #T is the temperature
vT = TestFunction(V)
q = Constant(0) #q is the volumetric heat source term
FT = T*vT*dx + dt*thermal_diffusivity*dot(grad(T), grad(vT))*dx - (T_n + dt*q)*vT*dx #This is the heat transfer equation
aT, LT = lhs(FT), rhs(FT) #Rearranging the equation





### Time-stepping
T = Function(V)
off_gassing=list()
c=Function(V)
if solve_diffusion==True:
  #fileC = File("Validation/Comparison between Cells in steel/"+str(cells)+"Cells/"+str(cells)+"Cells.pvd")
  fileC = File("Solutions/Compare/"+str(cells)+"Cells/"+str(cells)+"Cells.pvd")

if solve_temperature==True:
  fileT = File("solutions/Test/solutionT.pvd")

if calculate_off_gassing==True:
  file_off_gassing = "Solutions/Off-gassing/"+str(cells)+"Cells/"+str(cells)+"off_gassing.csv"

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
    fileC << (c,t)
  # Compute solution temperature
  if solve_temperature==True:
    solve(aT == LT, T, bcs_T)
    fileT << (T,t)

  # Update previous solution
  c_n.assign(c)
  T_n.assign(T)

  #Updating boundary conditions
  outside_bc_T.t += dt
  bco_T=DirichletBC(V,outside_bc_T,outside)
  bcs_T=list()
  bcs_T.append(bco_T)


  inside_bc_c.t += dt
  bci_c=DirichletBC(V,inside_bc_c,inside)
  bcs_c=list()
  bcs_c.append(bci_c)
  #bcs_c.append(bco_c)

  #Update the materials properties
  if solve_diffusion_coefficient_temperature_dependent==True and solve_temperature==True and solve_diffusion==True:
    for cell in cells(mesh):
      cell_no=cell.index()
      subdomain_no = subdomains.array()[cell_no]
      Ta=0
      for i in range(0,12,3):
        Ta+=T(cell.get_vertex_coordinates()[i],cell.get_vertex_coordinates()[i+1],cell.get_vertex_coordinates()[i+2])
      Ta=Ta/4
      D_value = calculate_D(Ta,subdomain_no)
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