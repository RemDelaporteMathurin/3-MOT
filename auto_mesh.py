

import sys
import os

import argparse

parser = argparse.ArgumentParser(description='make a mesh')
parser.add_argument("thickness")
args = parser.parse_args()

print(args)

steel_thickness = float(args.thickness)

print(steel_thickness)


arg_string = '"thickness=' + "'" +str(steel_thickness) + "'" +'"'

#print(arg_string)

os.system('Trelis-16.5 -batch -nographics make_mesh.py '+arg_string)



sys.exit()
