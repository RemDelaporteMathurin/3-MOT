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

T = 10.0           # final time
num_steps = 500    # number of time steps
dt = T / num_steps # time step size
mu = 1             # kinematic viscosity
rho = 1            # density

# Create mesh and define function spaces
mesh = UnitSquareMesh(40, 40)
V = VectorFunctionSpace(mesh, 'P', 2)
Q = FunctionSpace(mesh, 'P', 1)

tol = 1e-14
class Omega_0(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] <= 0.5 + tol


class Omega_1(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] >= 0.5 - tol

subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())

subdomain_0 = Omega_0()
subdomain_1 = Omega_1()
subdomain_0.mark(subdomains, 0)
subdomain_1.mark(subdomains, 1)


dx = Measure('dx', domain=mesh, subdomain_data=subdomains)





# Define boundaries
inflow  = 'near(x[0], 0)'
outflow = 'near(x[0], 1)'
walls   = 'near(x[1], 0) || near(x[1], 1)'

# Define boundaries whole demi-square
inflow  = 'near(x[0], 0) && x[1]<0.5'
outflow = 'near(x[0], 1) && x[1]<0.5'
walls   = 'near(x[1], 0) || near(x[1], 0.5) || near(x[1], 1) || (near(x[0], 0) && x[1]>0.5) || (near(x[0], 1) && x[1]>0.5) '

# Define boundary conditions
bcu_noslip  = DirichletBC(V, Constant((0, 0)), walls)
bcp_inflow  = DirichletBC(Q, Constant(20), inflow)
bcp_outflow = DirichletBC(Q, Constant(0), outflow)
bcu = [bcu_noslip]
bcp = [bcp_inflow, bcp_outflow]

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
n   = FacetNormal(mesh)
f   = Constant((0, 0))
k   = Constant(dt)
mu  = Constant(mu)
rho = Constant(rho)

# Define strain-rate tensor
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))

# Define variational problem for step 1
F1 = rho*dot((u - u_n) / k, v)*dx(0) + \
     rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx(0) \
   + inner(sigma(U, p_n), epsilon(v))*dx(0) \
   + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
   - dot(f, v)*dx(0)
a1 = lhs(F1)
L1 = rhs(F1)

# Define variational problem for step 2
a2 = dot(nabla_grad(p), nabla_grad(q))*dx(0)
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx(0) - (1/k)*div(u_)*q*dx(0)

# Define variational problem for step 3
a3 = dot(u, v)*dx(0)
L3 = dot(u_, v)*dx(0) - k*dot(nabla_grad(p_ - p_n), v)*dx(0)

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Apply boundary conditions to matrices
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

T = TrialFunction(Q)
T_n = Function(Q)
vT = TestFunction(Q)
FT = ((T-T_n)/k)*vT*dx  # Transient term
FT += 0.1*dot(grad(T), grad(vT))*dx # Diffusion (conduction term) 
FT += (dot(u_, grad(T)))*vT*dx(0) # Advection term
T_ = Function(Q)


bcsT = [DirichletBC(Q, Constant(8), inflow)]
output_file = File('solution.pvd')



# Time-stepping
t = 0
for n in range(num_steps):

    # Update current time
    t += dt

    # Step 1: Tentative velocity step
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(), b1)
#
    # Step 2: Pressure correction step
    b2 = assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p_.vector(), b2)
#
    # Step 3: Velocity correction step
    b3 = assemble(L3)
    solve(A3, u_.vector(), b3)

    solve(lhs(FT)==rhs(FT), T_, bcsT)
    T_n.assign(T_)

    output_file << (T_, t)
    output_file << (u_n, t)
    # Update previous solution
    u_n.assign(u_)
    p_n.assign(p_)
    #print(u_n(0.5,0.75))
