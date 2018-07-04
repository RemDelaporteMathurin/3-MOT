#!python
#!python
#!python
#!python
#!python
#!python
#!python
#!python
#!python
#!python
#This Cubit?Trellis script generates a 3-layer cube mesh (tets)
#The dimension of the layers are specified below
#To run this script use the command Cubit [name_of_the_file].py

#The mesh generated is in the Abaqus (.inp) format
#You may have to change the output file 


#Trelis-16.5 make_mesh.py "thickness='2e-3'"


aprepro_vars = cubit.get_aprepro_vars()
steel_thickness = float(cubit.get_aprepro_value_as_string("thickness"))


concrete_thickness=240e-3
polymer_thickness=20e-3
internal_cavity_dimension= 1738e-3
quality=7
number_of_layer=1

cubit.cmd('reset')



#makes vol1
cubit.cmd('brick x ' + str(internal_cavity_dimension+steel_thickness+polymer_thickness+concrete_thickness) + \ 
' y ' + str(internal_cavity_dimension+steel_thickness+polymer_thickness+concrete_thickness) + \
' z ' + str(internal_cavity_dimension+steel_thickness+polymer_thickness+concrete_thickness))

#makes vol2
cubit.cmd('brick x ' + str(internal_cavity_dimension+steel_thickness+polymer_thickness) + \ 
' y ' + str(internal_cavity_dimension+steel_thickness+polymer_thickness) + \
' z ' + str(internal_cavity_dimension+steel_thickness+polymer_thickness))


for i in range(0,number_of_layer+1):
  cubit.cmd('brick x ' + str(internal_cavity_dimension+((number_of_layer-i)/float(number_of_layer))*steel_thickness)+ \ 
  ' y ' + str(internal_cavity_dimension+((number_of_layer-i)/float(number_of_layer))*steel_thickness) + \
  ' z ' + str(internal_cavity_dimension+((number_of_layer-i)/float(number_of_layer))*steel_thickness))



for i in range(1,number_of_layer+3):
  cubit.cmd('subtract volume '+str(i+1)+' from volume '+str(i)+' keep')

for i in range(1,number_of_layer+4):
  cubit.cmd('delete Volume '+ str(i))


cubit.cmd('compress id all')
cubit.cmd('imprint body all')
#cubit.cmd('merge tolerance 1.e-6')
cubit.cmd('merge all')


cubit.cmd('volume all scheme Tetmesh')
cubit.cmd('volume all scheme tetmesh proximity layers off geometry approximation angle 15 geometric sizing on ')
cubit.cmd('volume all size auto factor '+str(quality)) # replace with sizing function 10 times smaller than smallest size in model
cubit.cmd('volume all tetmesh growth_factor 1 ')
cubit.cmd('Trimesher surface gradation 1.3')
cubit.cmd('Trimesher volume gradation 1.3')
cubit.cmd('mesh volume all')

for i in range(0,6):
  cubit.cmd('refine min_through_thickness 3 source surface '+ str(19+i)+' target surface '+str(31+i)+' anisotropic')
  cubit.cmd('refine min_through_thickness 3 source surface '+ str(1+i)+' target surface '+str(7+i)+' anisotropic')
  cubit.cmd('refine min_through_thickness 3 source surface '+ str(7+i)+' target surface '+str(19+i)+' anisotropic')

cubit.cmd('export abaqus "/Users/rdelapor/Documents/tritium_diffusion/geo/mesh.inp" dimension 3 overwrite')
#cubit.cmd('export step "/Users/rdelapor/Desktop/fenics_fusion_solvers/fenics_fusion_solvers/geometry_files/rcb/air.stp" volume 1 overwrite')
#cubit.cmd('export step "/Users/rdelapor/Desktop/fenics_fusion_solvers/fenics_fusion_solvers/geometry_files/rcb/concrete.stp" volume 2 overwrite')
#cubit.cmd('export step "/Users/rdelapor/Desktop/fenics_fusion_solvers/fenics_fusion_solvers/geometry_files/rcb/polymer.stp" volume 3 overwrite')














