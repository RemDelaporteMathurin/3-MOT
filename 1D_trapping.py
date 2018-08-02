from fenics import *
from dolfin import *
import numpy as np
import csv
import sys
import os
import argparse
import matplotlib.pyplot as plt

#os.system('dolfin-convert geo/mesh.inp geo/coucou.xml')
#parser = argparse.ArgumentParser(description='make a mesh')
#parser.add_argument("thickness")
#args = parser.parse_args()

#print(args)

#steel_thickness = float(args.thickness)



print('Getting the solvers')

print('Defining the solving parameters')
fluence=1e22 #D/m2
flux=2.5e19 #D/m2/s
implantation_time=fluence/flux
resting_time=50.0
TDS_time=50.0
Time =implantation_time+resting_time+TDS_time #60000*365.25*24*3600.0# final time
num_steps = int((implantation_time+resting_time+TDS_time))# number of time steps
dt = Time / num_steps # time step size
t=0 #Initialising time to 0s

print(Time)
print(num_steps)


print('Defining mesh')
nodes=500
size=3e-6
#Dx=size/nodes
mesh = IntervalMesh(nodes,0,size)

print('Defining Functionspaces')
V = FunctionSpace(mesh, 'P', 1) #FunctionSpace of the solution c
V0 = FunctionSpace(mesh, 'DG', 0) #FunctionSpace of the materials properties

### Define initial values
print('Defining initial values')

##Tritium concentration
iniC = Expression('0',degree=1)
#c_sol_n = interpolate(iniC, V)
c_sol_n = interpolate(Constant(0),V)
c_trap_n = interpolate(Constant(0),V)
c_trap2_n = interpolate(Constant(0.0),V)
c_trap3_n = interpolate(Constant(0.0),V)
n3_n=interpolate(Constant(0.0),V) #trap 3 density at t=0s
##Temperature
initial_temperature=273.15+14
iniT = Expression(str(initial_temperature),degree=1)
T_n = interpolate(iniT, V)


### Boundary Conditions
tol=1e-13
print('Defining boundary conditions')
def boundary_L(x, on_boundary):
    return on_boundary and (near(x[0], 0, tol))

def boundary_R(x, on_boundary):
    return on_boundary and (near(x[0], size, tol))
##Tritium concentration
inside_bc_c=Expression('0', t=0, degree=1) #0.4*exp(-t*log(2)/(12.3*365.25*24*3600))
outside_bc_c=Expression('0', t=0, degree=2)
bci_c=DirichletBC(V,inside_bc_c,boundary_L)
bco_c=DirichletBC(V,outside_bc_c,boundary_R)

g = Constant(0.0)

bcs_c=list()
bcs_c.append(bci_c)
bcs_c.append(bco_c)
##Temperature

def T_var(t):
  if t<implantation_time: 
    return 300
  elif t<implantation_time+resting_time: 
    return 300
  else:
     return 300+8*(t-(implantation_time+resting_time))
temp=Expression('t <= (implantation_time+resting_time) ? 300: 300+8*(t-(implantation_time+resting_time))',implantation_time=implantation_time,resting_time=resting_time,t=0,degree=2)

###Defining materials properties
print('Defining the materials properties')


k_B=8.6e-5 #Boltzann constant eV/K
density_W=6e28 #density of tungsten in m^-3
def calculate_D(T,subdomain_no):
  R=8.314 #Perfect gas constant
  if subdomain_no==0: #Steel
    return 4.1e-7*np.exp(-0.39/(k_B*T))/(2**0.5)

D=Function(V0)
for cell in cells(mesh):
  cell_no=cell.index()
  D.vector()[cell_no] = calculate_D(T_var(0),0)
  #print('coucou')

#print(D)



alpha=110e-12 #lattice constant ()
beta=6#*density_W#number of solute sites per atom (6 for W)
v_0=1e13#frequency factor s-1
E1=0.87 #in eV trap 1 activation energy
n1=1e-3#*density_W#trap 1 density
#n1=Expression('0+1*t',t=t,degree=1)#*density_W#trap 1 density

E2=1.0 #in eV activation energy
n2=4e-4#*density_W #trap 2 density
E3=1.5 #in eV activation energy
xp=1e-6
teta=Expression('x[0]<=xp ? 1/xp :0',xp=xp,degree=2)
n3amax=1e-1#*density_W
n3a=6e-4#*density_W
n3bmax=1e-2#*density_W
n3b=2e-4#*density_W
### Define variational problem
print('Defining the variational problem')
phi=Expression('t<=implantation_time ? flux: 0',implantation_time=implantation_time,flux=flux,t=0,degree=2)#this is the incident flux in D/m2/s      
r=0.56
fx = Expression('1/(width*pow(2*3.14,0.5))*exp(-0.5*(pow((x[0]-center)/width,2)))',center=4.5e-9,width=2.5e-9,degree=2)#  This is the tritium volumetric source term   -1/(1/3*e*pow(2*3.14,0.5))*exp(-0.5*(x[0]/pow(1/3*e,2)))   (x[0]<e/100 ? -2.5e19*(1-100*x[0]/e)
vc = TestFunction(V)
vc1 = TestFunction(V)
vc2 = TestFunction(V)

c_sol = Function(V)
c_trap=Function(V)
c_trap2=Function(V)
c_trap3=Function(V)
n3=Function(V)
A=Function(V)
B=Function(V)




c_sol = TrialFunction(V)#c is the tritium concentration
F1=((c_sol-c_sol_n)/dt)*vc*dx  +4.1e-7*exp(-0.39/(k_B*temp))/(2**0.5)*dot(grad(c_sol), grad(vc))*dx -(1-r)*phi*fx/density_W*vc*dx+ (c_trap-c_trap_n)/dt*vc*dx +(c_trap2-c_trap2_n)/dt*vc*dx#+(c_trap2-c_trap2_n)/dt*vc*dx+(c_trap3-c_trap3_n)/dt*vc*dx         4.1e-7*exp(-0.39/(k_B*temp))/(2**0.5)*dot(grad(c_sol), grad(vc))*dx    
ac1,Lc1= lhs(F1),rhs(F1)
c_sol=Function(V)

c_trap = TrialFunction(V)#c is the tritium concentration
F2=((c_trap-c_trap_n)/dt)*vc*dx -4.1e-7*exp(-0.39/(k_B*temp))/(2**0.5)/((alpha**2)*beta)*n1*c_sol_n*(1-c_trap/n1)*vc*dx  + c_trap*v_0*exp(-E1/(k_B*temp))*vc*dx#-1*vc1*dx#    4.1e-7*exp(-0.39/(k_B*300))/(2**0.5)/((alpha**2)*beta)*n1*
ac2,Lc2= lhs(F2),rhs(F2)
c_trap=Function(V)

print(4.1e-7*exp(-0.39/(k_B*300))/(2**0.5)/((alpha**2)*beta)*n1)


c_trap2 = TrialFunction(V)#c is the tritium concentration
F3=((c_trap2-c_trap2_n)/dt)*vc*dx - 4.1e-7*exp(-0.39/(k_B*temp))/(2**0.5)/(alpha**2*beta)*n2*c_sol_n*(1-c_trap2/n2)*vc*dx + c_trap2*v_0*exp(-E2/(k_B*temp))*vc*dx
ac3,Lc3= lhs(F3),rhs(F3)
c_trap2=Function(V)







#
#n3 = TrialFunction(V)#n3 is the density of trap 3
#Fn3=((n3-n3_n)/dt)*vc*dx-(1-r)*phi/density_W*((1-n3/n3amax)*n3a*fx+(1-n3/n3bmax)*n3b*teta)*vc*dx
#acn3,Lcn3=lhs(Fn3),rhs(Fn3)
#n3=Function(V)
#
#c_trap3 = TrialFunction(V)#c is the tritium concentration
#F4=((c_trap3-c_trap3_n)/dt)*vc*dx - 4.1e-7*exp(-0.39/(k_B*temp))/(2**0.5)/(alpha**2*beta)*c_sol_n*(n3_n-c_trap3)*vc*dx + c_trap3*v_0*exp(-E3/(k_B*temp))*vc*dx
#ac4,Lc4= lhs(F4),rhs(F4)
#c_trap3=Function(V)

### Time-stepping
Solute = File("Solutions/Trapping/Solute.pvd")
Trap1 = File("Solutions/Trapping/Trap1.pvd")
Trap2 = File("Solutions/Trapping/Trap2.pvd")
Trap3 = File("Solutions/Trapping/Trap3.pvd")
n3F = File("Solutions/Trapping/n3.pvd")
#fileT = File("solutions/Test/solutionT.pvd")
filedesorption="Solutions/Trapping/desorption.csv"
desorption=list()
total_sol_n=0
total_trap_n=0
total_n=0
tab=[]
for n in range(num_steps):
  # Update current time
  #print("t= "+str(t)+" s")
  #print("t= "+str(t/3600/24/365.25)+" years")
  #print(str(100*t/Time)+" %")  
  t += dt
  phi.t += dt
  temp.t +=dt

  # Compute solution concentration
  solve(ac1==Lc1,c_sol,bcs_c)  
  solve(ac2==Lc2,c_trap,bcs_c)  
  solve(ac3==Lc3,c_trap2,bcs_c)

  
  Solute << (c_sol,t)
  Trap1 << (c_trap,t)
  Trap2 << (c_trap2,t)
  

  print('c_trap = '+str(c_trap(size/2)))
  print('c_trap_n = '+str(c_trap_n(size/2)))
  print('var c_trap = '+str(c_trap(size/2)-c_trap_n(size/2)))

  print('c_sol = '+str(c_sol(size/2)))
  print('c_sol_n = '+str(c_sol_n(size/2)))
  print('var c_sol = '+str(c_sol(size/2)-c_sol_n(size/2)))


  A=(c_trap-c_trap_n)
  #print('Truc chelou = '+str((1-c_trap(size/2)/n1)))
  c_sol_n.assign(c_sol)
  c_trap_n.assign(c_trap)
  
  #

  #
  #solve(acn3==Lcn3,n3,bcs_c)
  #n3F << (n3,t)
  #
  #solve(ac4==Lc4,c_trap3,bcs_c)
  #Trap3 << (c_trap3,t)


  #Updating boundary conditions
  #bci_c=DirichletBC(V,inside_bc_c,boundary)
  #bcs_c=list()
  #bcs_c.append(bci_c)
  #bcs_c.append(bco_c)


  total_trap1=assemble(c_trap*dx)
  total_trap2=assemble(c_trap2*dx)
  total_trap3=assemble(c_trap3*dx)
  total_trap=total_trap1+total_trap2+total_trap3
  total_sol=assemble(c_sol*dx)
  total=total_trap+total_sol
  desorption_rate=[-(total-total_n)*density_W/dt,T_var(t),t]
  total_n=total
  print('Augmentation sol = '+str((total_sol-total_sol_n)/dt))
  print('Augmentation trap = '+str((total_trap-total_trap_n)/dt))
  

  total_sol_n=total_sol
  total_trap_n=total_trap
  if t==implantation_time+resting_time:
    valeur=(total)*density_W
  if t>implantation_time+resting_time:
    print(valeur)
    desorption.append(desorption_rate)
    #print('Total of D soluble = ' + str(total_sol))
    #print('Total of D traped 1= ' + str(total_trap1))
    #print('Total of D traped 2= ' + str(total_trap2))
    #print('Total of D traped 3= ' + str(total_trap3))
    #print("Total of D = " + str(total))
    #print("Desorption rate = " + str(desorption_rate))


  value=calculate_D(T_var(t),0)
  # Update previous solution
  for cell in cells(mesh):
    cell_no=cell.index()
    D.vector()[cell_no] = value 


  c_trap2_n.assign(c_trap2)
  c_trap3_n.assign(c_trap3)
  n3_n.assign(n3)
  print('T = ' + str(T_var(t))+' K')
  print('D = ' +str(calculate_D(T_var(t),0)))


with open(filedesorption, "w") as output:
    writer = csv.writer(output, lineterminator='\n')
    writer.writerows(['dTt'])
    for val in desorption:
        writer.writerows([val])#

print('Total retention before TDS = '+ str(valeur) + 'D/m2')