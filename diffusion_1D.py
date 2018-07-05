from fenics import *
from dolfin import *
import numpy as np
Time = 50.0  # final time
num_steps = 50 # number of time steps
dt = Time / num_steps # time step size 6912000000000

# Create mesh and define function space
nodes=100
size=3e-6
#Dx=size/nodes
mesh = IntervalMesh(nodes,0,size)
V = FunctionSpace(mesh, 'P', 1) #FunctionSpace of the solution c
V0 = FunctionSpace(mesh, 'DG', 0) #FunctionSpace of the materials properties

#Defining the subdomains

concrete_thickness=240e-3
steel_thickness=2e-3
polymer_thickness=20e-3
air_thickness=15e-3 #This is the external layer of air that will be simulated
internal_cavity_dimension=1738e-3




# Boundary Conditions
tol = 1e-13

def boundary_L(x, on_boundary):
    return on_boundary and (near(x[0], 0, tol))

def boundary_R(x, on_boundary):
    return on_boundary and (near(x[0], size, tol))

inside = 0.0
outside = 0.0
bci=DirichletBC(V,inside,boundary_L)
bco=DirichletBC(V,outside,boundary_R)
bcs=[bci,bco]

##Defining materials properties
D  = 1.16e-13

# Define initial value
iniC = Expression('2.26e24',degree=1)
c_n = interpolate(iniC, V)


# Define variational problem
c = TrialFunction(V) #c is the tritium concentration
v = TestFunction(V)
f = Constant(0) #This is the tritium volumetric source term

F = c*v*dx + dt*D*dot(grad(c), grad(v))*dx - (c_n + dt*f)*v*dx #This is the tritium diffusion equation
a, L = lhs(F), rhs(F) #Rearranging the equation


# Time-stepping
c = Function(V)
t=0
fileC = File("Solutions/Test/"+str(nodes)+"nodes.pvd")
#File("Validation/"+str(nodes)+ " nodes "+ str(num_steps) + " time steps 1e8dy/solutions/D="+ str(D)+"/D="+ str(D)+".pvd")

for n in range(num_steps):
    # Update current time
    print("t= "+str(t)+" s")
    print(str(100*t/Time)+" %")
    t += dt
    # Compute solution
    solve(a == L, c, bcs)
    fileC << (c,t)
    # Update previous solution
    c_n.assign(c)