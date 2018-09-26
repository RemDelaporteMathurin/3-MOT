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
import inspect
#from tqdm import *
#from random import random, randint
#from time import sleep
import math
import matplotlib.pyplot as plt

mot = __import__('3-MOT')


def return_distance_to_surface(a, b, c, d, x, y, z):
    distance = abs(a * (x ) + b * (y) + c * (z) + d) / ((a*a+b*b+c*c)**0.5)
    #print(distance)
    return distance

r = return_distance_to_surface(a = 0.9908270,
                           b = 0.115811,
                           c = -0.06940,
                           d = -12.1610100263915,
                           x = 12.110082,
                           y = 1.411025,
                           z = 0.020152)

D = - 0.9908270 * 12.110082 - 0.1158111 * 1.411025 - (-0.06940) * 0.020152
print('D = ', D)

data = mot.get_databases('MOT_parameters_breeder_blankets (copy).json')
mesh, n0, volume_marker, dx, surface_marker, ds, mesh_fluid, volume_marker_fluid, dx_fluid, surface_marker_fluid, ds_fluid, n_fluid = mot.define_mesh(data, False)

V0 = FunctionSpace(mesh, 'DG', 0)
Q = Function(V0)
n = 0
for cell in cells(mesh):
    print(n , round(100*n/len(volume_marker),1), end='% \r')
    material_id = mot.which_material_is_it(volume_marker.array()[cell.index()],data)
    coord = (cell.get_vertex_coordinates()[0], cell.get_vertex_coordinates()[1], cell.get_vertex_coordinates()[2])
    #print(coord[0])
    #print(material_id)
    #print(coord)
    r_global = ((coord[0]**2+coord[2]**2+coord[1]**2)**0.5)
    r = return_distance_to_surface(a = 0.9908270,
                           b = 0.115811,
                           c = -0.06940,
                           d = -12.1610100263915,
                           x = coord[0],
                           y = coord[1],
                           z = coord[2])
    if material_id == 'eurofer':
        #Q.vector()[cell.index()] = 1      
        Q.vector()[cell.index()] = 7.534822e6*np.exp(-8.988724*r)
    elif material_id == 'lithium_lead':
        #Q.vector()[cell.index()] = 2      
        Q.vector()[cell.index()] = 9.460241e6*np.exp(-6.198768*r)
    elif material_id == 'tungsten':
        Q.vector()[cell.index()] = 23.20664e6*np.exp(-71.74517*r)
    else:
        Q.vector()[cell.index()] = 0
    n += 1


File('volumetric_source_term_quality1.xml') << Q
File('volumetric_source_term_mesh_quality1.xml') << mesh
File('volumetric_source_term_quality1.pvd') << Q