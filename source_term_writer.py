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
from materials_properties import *
import inspect
#from tqdm import *
#from random import random, randint
#from time import sleep
import math
import matplotlib.pyplot as plt

mot = __import__('3-MOT')



data = mot.get_databases('MOT_parameters_breeder_blankets.json')
mesh, xdmf_in, n0 = mot.define_mesh(data)
volume_marker, dx = mot.get_volume_markers(mesh, xdmf_in)




V0 = FunctionSpace(mesh, 'DG', 0)
Q = Function(V0)

for cell in cells(mesh):
    material_id = mot.which_material_is_it(volume_marker.array()[cell.index()],data)
    coord = (cell.get_vertex_coordinates()[0], cell.get_vertex_coordinates()[1], cell.get_vertex_coordinates()[2])
    #print(coord[0])
    #print(material_id)
    r = ((coord[0]**2+coord[2]**2+coord[1]**2)**0.5-6.3)
    if material_id == 'eurofer':
        #Q.vector()[cell.index()] = 1      
        Q.vector()[cell.index()] = 16.08e6*np.exp(-11.3*r)
    elif material_id == 'lithium_lead':
        #Q.vector()[cell.index()] = 2      
        Q.vector()[cell.index()] = 9.96e6*np.exp(-5.23*r) + 19.27e6*np.exp(-20.56*r)
    elif material_id == 'tungsten':
        Q.vector()[cell.index()] = 0
    else:
        Q.vector()[cell.index()] = 10


File('volumetric_source_term.xml') << Q
File('volumetric_source_term_mesh.xml') << mesh
File('volumetric_source_term.pvd') << Q

mesh = Mesh('volumetric_source_term_mesh.xml')
plot(mesh)
plt.show()
