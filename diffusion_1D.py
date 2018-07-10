from fenics import *
from dolfin import *
import numpy as np
Time = 400.0  # final time
num_steps = 400 # number of time steps
dt = Time / num_steps # time step size 6912000000000

# Create mesh and define function space
nodes=500
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

inside = 0
outside = 0.0
bci=DirichletBC(V,inside,boundary_L)
bco=DirichletBC(V,outside,boundary_R)
bcs=[bci,bco]

##Defining materials properties
D  = 8.13e-14#4.1e-7*np.exp(-0.39/(8.6e-5*300))/(2**0.5)/1400

# Define initial value
iniC = Expression('0',degree=1)
c_n = interpolate(iniC, V)


# Define variational problem
c = TrialFunction(V) #c is the tritium concentration
v = TestFunction(V)
f = Constant((1-0.56)*2.5e19/6e28) #This is the tritium volumetric source term
fx = Expression('1/(width*pow(2*3.14,0.5))*exp(-0.5*(pow((x[0]-center)/width,2)))',center=5e-9,width=2e-9,degree=2)#5e-9
F = c*v*dx + dt*D*dot(grad(c), grad(v))*dx - (c_n + dt*f*fx)*v*dx #This is the tritium diffusion equation
a, L = lhs(F), rhs(F) #Rearranging the equation


# Time-stepping
c = Function(V)
t=0
fileC = File("Solutions/Test/test.pvd")
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