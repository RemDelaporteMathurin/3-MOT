from fenics import *
from dolfin import *
import numpy as np
import csv


print('Getting the solvers')
solve_temperature=False
solve_diffusion=True
solve_diffusion_coefficient_temperature_dependent=False
solve_with_decay=False
calculate_off_gassing=False

print('Defining the solving parameters')
Time =8.64e12      # final time
num_steps = 5 # number of time steps
dt = Time / num_steps # time step size
t=0 #Initialising time to 0s

print('Defining mesh')
#Create mesh and define function space
#mesh=Mesh('geo/mesh_rcb_12_nodes_in_steel.xml')
#mesh=Mesh('geo/mesh_3D_test_validation.xml')
mesh=BoxMesh(Point(0,0,0), Point(10,10,0.1), 20, 20, 30)

#WARNING : MAKE SURE THESE DIMENSIONS ARE THE SAME AS IN THE CUBIT SCRIPT
concrete_thickness=0.1
steel_thickness=0
polymer_thickness=0
internal_cavity_dimension=1738e-2

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

inside =  CompiledSubDomain("x[2]<=0.001 && on_boundary")
outside = CompiledSubDomain("x[2]>=-0.001+0.1 && on_boundary")

##Tritium concentration
inside_bc_c=Expression('0', t=0, degree=1) #t<12.3*365*24*3600.0 ? 1.0 :1.0
outside_bc_c=Expression('1', t=0, degree=2)
bci_c=DirichletBC(V,inside_bc_c,inside)
bco_c=DirichletBC(V,outside_bc_c,outside)
#g = Function(V)
g = Constant(0.0)
#g=conditional(gt(c_n, 0), 1e-10*c_n**0.74, Constant(0.0))#
bcs_c=list()
bcs_c.append(bci_c)
bcs_c.append(bco_c)
##Temperature
inside_bc_T=0
outside_bc_T=Expression('14+273.15+7*cos(2*3.14*t/365.25/24/3600)+16*cos(2*3.14*t/24/3600)', t=0, degree=2)
bci_T=DirichletBC(V,inside_bc_T,inside)
bco_T=DirichletBC(V,outside_bc_T,outside)
bcs_T=list()
bcs_T.append(bco_T)




D=Constant(5e-16)

### Define variational problem
print('Defining the variational problem')
c = Function(V)#c is the tritium concentration
vc = TestFunction(V)
f = Expression('0',t=0,degree=2)#This is the tritium volumetric source term 
F=((c-c_n)/dt)*vc*dx + D*dot(grad(c), grad(vc))*dx + f*vc*dx +g*vc*ds
#F=((c-c_n)/dt)*vc*dx + (f+decay*c_n)*vc*dx +g*vc*ds 
ac,Lc= lhs(F),rhs(F)

#c= Function(V0)

### Time-stepping

if solve_diffusion==True:
  fileC = File("solutions/Validation/solutionC.pvd")

if solve_temperature==True:
  fileT = File("solutions/Validation/solutionT.pvd")

if calculate_off_gassing==True:
  file_off_gassing = "solutions/Test/off_gassing.csv"

for n in range(num_steps):

  
  # Update current time
  print("t= "+str(t)+" s")
  print(str(100*t/Time)+" %")
  t += dt
  f.t += dt

  # Compute solution concentration
  if solve_diffusion==True:
    problem_diffusion = NonlinearVariationalProblem(F, c, bcs_c, derivative(F,c))
    solver_diffusion  = NonlinearVariationalSolver(problem_diffusion)
    prm_diffusion = solver_diffusion.parameters
    prm_diffusion['newton_solver']['absolute_tolerance'] = 1E-20
    prm_diffusion['newton_solver']['relative_tolerance'] = 1E-3
    solver_diffusion.solve()
    #solve(F==0,c,J=derivative(F,c),bcs=bcs_c)
    #solve(ac==Lc,c,bcs_c)
    fileC << (c,t)

  # Update previous solution
  c_n.assign(c)
