from __future__ import print_function, division
from fenics import *

parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['optimize'] = True

Reynolds = 400
density = 1
mu = density/Reynolds
print(mu)

N = 100
#mesh = RectangleMesh(Point(0, -1), Point(10, 1), 5*N, N, 'crossed')
#inflw_bndr   = 'near(x[0], 0.0)'
#outflw_bndr  = 'near(x[0], 10.0)'
#no_slip_bndr = 'near(x[1], -1) || near(x[1], 1)'


no_slip_bndr = ' near(x[1], 0) || near(x[0], 0) || near(x[0], 1)'
driven_bndr = 'near(x[1],1)'

mesh = UnitSquareMesh(N, N, 'crossed')


### Space definition
V  = VectorElement('P', mesh.ufl_cell(), 2)
P  = FiniteElement('P', mesh.ufl_cell(), 1)
TH = MixedElement([V, P])
W  = FunctionSpace(mesh, TH)

### Define Boundary Conditions

no_slip_bc  = DirichletBC(W.sub(0), Constant((0.0, 0.0)), no_slip_bndr)

#p_outflw_bc = DirichletBC(W.sub(1), Constant(0.0), outflw_bndr)


driven_bc = DirichletBC(W.sub(0), Constant((1.0, 0.0)), driven_bndr)

bcs = [no_slip_bc, driven_bc]
#bcs = [no_slip_bc, u_inflw_bc, p_outflw_bc]

### Steady Part of the Momentum Equation
def steady(u):
    T = -p*I + 2*mu*sym(grad(u))
    return density*(inner(grad(u)*u, v_) + inner(T, grad(v_)) - inner(f, v_)) * dx

### Unknown and test functions
(v_, p_) = TestFunctions(W)
w = Function(W)
(v, p) = split(w)

f = Constant((0.0, 0.0))

I = Identity(v.geometric_dimension())

F = steady(v) + p_*div(v)*dx

J = derivative(F, w)
problem = NonlinearVariationalProblem(F, w, bcs, J)
solver  = NonlinearVariationalSolver(problem)

### Create Files for Storing the Solution
vfile = File('results/2d-stationary-channel-ns-velocity.pvd')
    
### Compute Solution
solver.solve()

### Extract Solutions
(v, p) = w.split()
    

vfile << v