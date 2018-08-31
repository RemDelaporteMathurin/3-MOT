"""
FEniCS tutorial demo program: Incompressible Navier-Stokes equations
for channel flow (Poisseuille) on the unit square using the
Incremental Pressure Correction Scheme (IPCS).

  u' + u . nabla(u)) - div(sigma(u, p)) = f
                                 div(u) = 0
"""

from __future__ import print_function
from fenics import *
import numpy as np

T = 100.0            # final time
num_steps = 500     # number of time steps
dt = T / num_steps  # time step size
mu = 1              # kinematic viscosity
rho = 1             # density

# Create mesh and define function spaces
mesh = UnitCubeMesh(10, 10, 10)

tol = 1e-14


class Fluid(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] >= 0.25 and x[1] <= 0.75 and x[2] >= 0.25  and x[2] <= 0.75


class SolidTop(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] >= 0.75


class SolidBottom(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] <= 0.25


class SolidLeft(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] <= 0.25


class SolidRight(SubDomain):
    def inside(self, x, on_boundary):
        return x[2] >= 0.75

subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())

fluid = Fluid()
solidTop = SolidTop()
solidBottom = SolidBottom()
solidLeft = SolidLeft()
solidRight = SolidRight()

fluid.mark(subdomains, 0)
solidTop.mark(subdomains, 1)
solidBottom.mark(subdomains, 1)
solidLeft.mark(subdomains, 1)
solidRight.mark(subdomains, 1)


fluid = SubMesh(mesh, subdomains, 0)

W = FunctionSpace(mesh, 'P', 1)  # FS of temperature
Q = FunctionSpace(fluid, 'P', 1)  # FS of pressure
V = VectorFunctionSpace(fluid, 'P', 2)  # FS of velocity

dx = Measure('dx', domain=mesh, subdomain_data=subdomains)

dx_fluid = Measure('dx', domain=fluid)


## Define boundaries demi-square

class Inflow(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0) and x[1] >= 0.2 and x[1] <= 0.8 and x[2] >= 0.2  and x[2] <= 0.8


class Outflow(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 1) and x[1] >= 0.2 and x[1] <= 0.8 and x[2] >= 0.2  and x[2] <= 0.8 - DOLFIN_EPS


class Walls(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[1], 0.2) or near(x[1], 0.8) or near(x[2], 0.2) or near(x[2], 0.8))

inflow = Inflow()
outflow = Outflow()
walls = Walls()

boundaries_fluid = MeshFunction("size_t", fluid, fluid.topology().dim() - 1)
boundaries_fluid.set_all(0)
inflow.mark(boundaries_fluid, 0)
outflow.mark(boundaries_fluid, 1)
walls.mark(boundaries_fluid, 2)

ds = Measure('ds', domain=mesh)

ds_fluid = Measure('ds', domain=fluid)

# Define boundary conditions
bcu_noslip = DirichletBC(V, Constant((0, 0, 0)), walls)

bcp_outflow = DirichletBC(Q, Constant(0), outflow)
bcp_inflow = DirichletBC(Q, Constant(8), inflow)
bcu = [bcu_noslip]
bcp = [bcp_inflow,bcp_outflow]

# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# Define functions for solutions at previous and current time steps
u_n = Function(V)
u_  = Function(V)
p_n = Function(Q)
p_  = Function(Q)

# Define expressions used in variational forms
U   = 0.5*(u_n + u)
n   = FacetNormal(fluid)
f   = Constant((0, 0, 0))
k   = Constant(dt)
mu  = Constant(mu)
rho = Constant(rho)

# Define strain-rate tensor
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))


print('Area fluid', assemble(1*ds_fluid))
print('Area total', assemble(1*ds))


## Define variational problem for step 1
F1 = rho*dot((u - u_n) / k, v)*dx_fluid \
   + rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx_fluid \
   + inner(sigma(U, p_n), epsilon(v))*dx_fluid \
   + dot(p_n*n, v)*ds_fluid - dot(mu*nabla_grad(U)*n, v)*ds_fluid \
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
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]


# Define variational problem heat diffusion
T = TrialFunction(W)
T_n = Function(W)
vT = TestFunction(W)
FT = ((T-T_n)/k)*vT*dx   # Transient term
FT += 0.01*dot(grad(T), grad(vT))*dx(0)  # Diffusion (conduction term) in fluid
FT += 0.001*dot(grad(T), grad(vT))*dx(1)  # Diffusion (conduction term) in solid
FT += (dot(u_, grad(T)))*vT*dx(0)  # Advection term in fluid
T_ = Function(W)



bcsT = [DirichletBC(W, Constant(8), inflow)]
output_file = File('Problems/Square_Pipe_Non_Automatic/Solution/solution.pvd')



# Time-stepping
t = 0
for n in range(num_steps):

    # Update current time
    t += dt

    # Step 1: Tentative velocity step
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(), b1)

    # Step 2: Pressure correction step
    b2 = assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p_.vector(), b2)

    # Step 3: Velocity correction step
    b3 = assemble(L3)
    solve(A3, u_.vector(), b3)

    solve(lhs(FT)==rhs(FT), T_, bcsT)
    T_n.assign(T_)

    output_file << (T_, t)
    output_file << (u_, t)
    output_file << (p_, t)
    # Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)