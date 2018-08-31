from fenics import *
from dolfin import *
import numpy as np
import csv
import sys
import os
import argparse
import json
import ast
from pprint import pprint


mesh=UnitSquareMesh(6,6)
V = FunctionSpace(mesh, 'P', 1)


class Bottom(SubDomain):
    def inside(self,x,on_boundary):
        if x[1]==0:
            return True
        else:
            return False
class Top(SubDomain):
    def inside(self,x,on_boundary):
        if x[1]==1:
            return True
        else:
            return False

class Left(SubDomain):
    def inside(self,x,on_boundary):
        if x[0]==0:
            return True
        else:
            return False
class Right(SubDomain):
    def inside(self,x,on_boundary):
        if x[0]==1:
            return True
        else:
            return False


domains = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
domains.set_all(0)

Bottom=Bottom()
Bottom.mark(domains,0)

Top=Top()
Top.mark(domains,0)

Left=Left()
Left.mark(domains,0)

Right=Right()
Right.mark(domains,1)

ds = Measure('ds', domain=mesh, subdomain_data = domains)

bcs=[]
bcs.append(DirichletBC(V,0,Bottom))
bcs.append(DirichletBC(V,0,Top))
bcs.append(DirichletBC(V,0,Left))



fileT=File("validationT.pvd")

T = TrialFunction(V) #T is the temperature
vT = TestFunction(V)

  
FT = 100*dot(grad(T), grad(vT))*dx - vT * 100*ds(1)
aT, LT = lhs(FT), rhs(FT)
T=Function(V)
solve(aT==LT,T,bcs)
fileT<<(T,1.0)
