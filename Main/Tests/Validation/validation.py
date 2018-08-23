from fenics import *
import numpy as np
 
N=25

mesh = UnitCubeMesh(N,N,N)
subdomains = CellFunction("size_t", mesh)

print('Number of cell is '+ str(len(subdomains.array())))
compute_error=True

time_start=time()
Time = 60.0       # final time
num_steps = 50   # number of time steps
dt = Time / num_steps # time step size
t=0 #Initialising time to 0s


decay=np.log(2)/(12.33*365*24*3600)
diff=0.01


# Create mesh and define function space
mesh = UnitCubeMesh(N,N,N)
V = FunctionSpace(mesh, 'P', 1)

iniC = Expression('1 + x[0]*x[0] +x[1]*x[1]+ 2*x[2]*x[2]',degree=6)
c_n = interpolate(iniC, V)

# Define boundary condition
c_D = Expression('1 + x[0]*x[0] +x[1]*x[1]+ 2*x[2]*x[2]+t',t=0, degree=6)
c_d=Function(V)

def boundary(x, on_boundary):
    return on_boundary
bc = DirichletBC(V, c_D, boundary)
# Define variational problem
c = TrialFunction(V)
v = TestFunction(V)
f = Expression('1-8*D+ d*(1+x[0]*x[0] +x[1]*x[1]+ 2*x[2]*x[2])',D=diff,d=decay,t=0, degree=6)
a = dot(grad(c), grad(v))*dx
L = f*v*dx

FC=c*v*dx+diff*dt*dot(grad(c), grad(v))*dx+dt*(decay*c_n*v*dx-f*v*dx)-c_n*v*dx
ac,Lc=lhs(FC),rhs(FC)
c = Function(V)
difference=Function(V)
num=File('Validation/3D/Numerical Solution.pvd')
real=File('Validation/3D/Manufactured Solution.pvd')
error=File('Validation/3D/Error.pvd')

error_L2=0

for n in range(num_steps):
  t += dt
  f.t = t
  print("t= "+str(t)+" s")
  print(str(100*t/Time)+" %")
  c_D.t += dt
  c_d.interpolate(c_D)
  bc = DirichletBC(V, c_D, boundary)
  # Compute solution
  
  solve(ac == Lc, c, bc)
  
  
  real << (c_d,t)
  num << (c,t)
  if compute_error==True:
    # Compute error in L2 norm
    error_L2 += errornorm(c_D, c, 'L2')
    # Compute maximum error at vertices
    #vertex_values_c_D = c_D.compute_vertex_values(mesh)
    #vertex_values_c = c.compute_vertex_values(mesh)

    #error_max = np.max(np.abs(vertex_values_c_D - vertex_values_c))
    #error_average = np.average(np.abs(vertex_values_c_D - vertex_values_c))
    # Print errors
    #print('error_L2  =', error_L2)
    #print('error_max =', error_max)

  c_n.assign(c)

time_end=time()
print('Computation time is '+ str(time_end-time_start)+' s')
if compute_error==True:
  error_L2=error_L2/num_steps
  #print('Average error is '+str(error_average))
  print('L2 error is '+str(error_L2))

print('Number of cell is '+ str(len(subdomains.array())))