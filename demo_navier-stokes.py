"""This demo program solves the incompressible Navier-Stokes equations
on an L-shaped domain using Chorin's splitting method."""

# Copyright (C) 2010-2011 Anders Logg
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Mikael Mortensen 2011
#
# First added:  2010-08-30
# Last changed: 2011-06-30

# Begin demo

from __future__ import print_function
from dolfin import *

# Print log messages only from the root process in parallel
parameters["std_out_all_processes"] = False;

# Load mesh from file
mesh = Mesh()
xdmf_in = XDMFFile(mesh.mpi_comm(), 'mesh_and_markers_tube.xdmf')
xdmf_in.read(mesh)
# prepare output file for writing by writing the mesh to the file
# xdmf_out = XDMFFile(str(data['mesh_file']).split('.')[1]+'_from_fenics.xdmf')
# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)

# Define trial and test functions
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

# Set parameter values
dt = 0.01
T = 3.0
nu = 0.01

# Define time-dependent pressure boundary condition
p_in = Expression("sin(3.0*t)", t=0.0, degree=2)


print('Marking the surfaces')
# read in Surface markers for boundary conditions
surface_marker_mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
xdmf_in.read(surface_marker_mvc, "surface_marker")
surface_marker_mvc.rename("surface_marker", "surface_marker")
# xdmf_out.write(surface_marker_mvc, xdmf_encoding)
surface_marker = MeshFunction("size_t", mesh, surface_marker_mvc)
ds = Measure('ds', domain=mesh, subdomain_data=surface_marker)


# Define boundary conditions
noslip3  = DirichletBC(V, (0, 0, 0), surface_marker, 3)
noslip4  = DirichletBC(V, (0, 0, 0), surface_marker, 4)
noslip5  = DirichletBC(V, (0, 0, 0), surface_marker, 5)
noslip6  = DirichletBC(V, (0, 0, 0), surface_marker, 6)
noslip11  = DirichletBC(V, (0, 0, 0), surface_marker, 11)
noslip12  = DirichletBC(V, (0, 0, 0), surface_marker, 12)
noslip13  = DirichletBC(V, (0, 0, 0), surface_marker, 13)
noslip14  = DirichletBC(V, (0, 0, 0), surface_marker, 14)
noslip15  = DirichletBC(V, (0, 0, 0), surface_marker, 15)
noslip16  = DirichletBC(V, (0, 0, 0), surface_marker, 16)



inflow  = DirichletBC(Q, 1, surface_marker, 1)
outflow = DirichletBC(Q, 0, surface_marker, 2)
bcu = [noslip3,noslip4,noslip5,noslip6,noslip11,noslip12,noslip13,noslip14,noslip15,noslip16]
bcp = [inflow, outflow]

# Create functions
u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

# Use nonzero guesses - essential for CG with non-symmetric BC
parameters['krylov_solver']['nonzero_initial_guess'] = True

# Create files for storing solution
ufile = File("results/velocity.pvd")
pfile = File("results/pressure.pvd")


# Time-stepping
t = dt
while t < T + DOLFIN_EPS:

    # Update pressure boundary condition
    p_in.t = t

    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "bicgstab", "default")
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    [bc.apply(p1.vector()) for bc in bcp]
    solve(A2, p1.vector(), b2, "bicgstab", prec)
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "bicgstab", "default")
    end()

    # Save to file
    ufile << u1
    pfile << p1
    print(u1(0.0015/2,0.0015/2,0.0015/2))
    print(u1(0,0,0))
    print(u1(0.06,0.06,0.1))


    # Move to next time step
    u0.assign(u1)
    t += dt
    print("t =", t)

# Hold plot
#interactive()
