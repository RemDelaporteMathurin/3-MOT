from fenics import *
import numpy as np
import csv


N=10

mesh = UnitCubeMesh(N,N,N)

Time = 3*12.33*365.25*24*3600.0     # final time
num_steps = 100   # number of time steps
dt = Time / num_steps # time step size
t=0 #Initialising time to 0s


decay=np.log(2)/(12.33*365.25*24*3600)
diff=0.00


# Create mesh and define function space
mesh = UnitCubeMesh(N,N,N)
V = FunctionSpace(mesh, 'P', 1)

iniC = Expression('1',degree=6)
c_n = interpolate(iniC, V)

# Define boundary condition
c_D = Expression('1',t=0, degree=6)
c_d=Function(V)

def boundary(x, on_boundary):
    return on_boundary
bc = DirichletBC(V, c_D, boundary)
# Define variational problem
c = TrialFunction(V)
v = TestFunction(V)

FC=c*v*dx+dt*(decay*c*v*dx)-c_n*v*dx
ac,Lc=lhs(FC),rhs(FC)
c = Function(V)
num=File('Validation/3D decay/Solution.pvd')


file_decay = "Validation/3D decay/decay.csv"


values=[]

for n in range(num_steps):
  t += dt
  print("t= "+str(t)+" s")
  print(str(100*t/Time)+" %")
  c_D.t += dt
  c_d.interpolate(c_D)
  bc = DirichletBC(V, c_D, boundary)
  # Compute solution
  
  solve(ac == Lc, c, bc)
  num << (c,t)
  c_n.assign(c)

  values.append([t/24/3600/365,c(0.5,0.5,0.5)])


with open(file_decay, "w") as output:
    writer = csv.writer(output, lineterminator='\n')
    writer.writerows(['tc'])
    writer.writerows(values)







