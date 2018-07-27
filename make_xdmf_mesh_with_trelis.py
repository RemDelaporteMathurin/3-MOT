import json
from pprint import pprint
import ast

# run this script with the following commands
# trelis -nographics -batch make_xdmf_mesh_with_trelis.py "quality='10'" "outputfile='mesh_and_markers.xdmf'" "structure='slice_armour_mod1.step,slice_back_lithium_lead_mod1.step,slice_back_plate_1_mod1.step,slice_back_plate_2_mod1.step,slice_back_plate_3_mod1.step,slice_cooling_plate_material_mod1.step,slice_first_wall_material_mod1.step,slice_lithium_lead_mod1.step,'" "coolant='slice_first_wall_coolant_mod1.step,slice_cooling_plate_coolant_mod1.step,slice_back_helium_mod1.step'"
# trelis make_xdmf_mesh_with_trelis.py "quality='10'" "outputfile='mesh_and_markers.xdmf'" "structure='slice_armour_mod1.step,slice_back_lithium_lead_mod1.step,slice_back_plate_1_mod1.step,slice_back_plate_2_mod1.step,slice_back_plate_3_mod1.step,slice_cooling_plate_material_mod1.step,slice_first_wall_material_mod1.step,slice_lithium_lead_mod1.step,'" "coolant='slice_first_wall_coolant_mod1.step,slice_cooling_plate_coolant_mod1.step,slice_back_helium_mod1.step'"

# trelis make_xdmf_mesh_with_trelis.py "json_input='name_of_the_json_file'"

def byteify(input):
    if isinstance(input, dict):
        return {byteify(key): byteify(value)
                for key, value in input.iteritems()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input

def read_arguments_from_json_file(aprepro_vars):
  print('reading arguments from json file')
  json_input = str(cubit.get_aprepro_value_as_string("json_input"))
  print('json_inputfile ='+str(json_input))
  with open(json_input) as f:
    data = json.load(f)
  data=byteify(data)
  pprint(data)
  quality= data["quality"]
  outputfile = data["mesh_file"]
  structure_and_materials=data["structure_and_materials"]
  structure_step_files=[]
  materials=[]
  material_ids=[]
  for entry in structure_and_materials:
    structure_step_files.append(entry['step_file'])
    print(entry['step_file'])
    materials.append(entry['material'])
    material_ids.append(entry['material_id'])
  return quality, outputfile, structure_step_files, materials, material_ids

def find_external_surfaces():
  print('looking for merged surfaces')
  surfaces_in_all_volumes = cubit.parse_cubit_list("surface"," in volume all ")  
  list_of_merged_surfaces = []
  list_of_unmerged_surfaces= []
  for surface in surfaces_in_all_volumes:
      is_merged = cubit.is_merged("surface", surface)
      if is_merged == True:
        list_of_merged_surfaces.append(surface)
      elif is_merged == False:
        cubit.cmd('color surface '+str(surface)+ ' red')
        list_of_unmerged_surfaces.append(surface)
  print('list_of_external_surfaces',list_of_unmerged_surfaces)
  return list_of_unmerged_surfaces


def mesh_and_remesh_till_done(quality):
  cubit.cmd('volume all scheme Tetmesh')
  cubit.cmd('volume all scheme tetmesh proximity layers off geometry approximation angle 15 geometric sizing on ')
  cubit.cmd('volume all size auto factor ' + str(quality)) # replace with sizing function 10 times smaller than smallest size in model
  cubit.cmd('volume all tetmesh growth_factor 1 ')
  cubit.cmd('Trimesher surface gradation 1.3')
  cubit.cmd('Trimesher volume gradation 1.3')
  cubit.cmd('mesh volume all')



aprepro_vars = cubit.get_aprepro_vars()

print("Found the following aprepro variables:")
print(aprepro_vars)
for var_name in aprepro_vars:
  val = cubit.get_aprepro_value_as_string(var_name)
  print("{0} = {1}".format(var_name, val))

quality, outputfile, structure_step_files, materials, material_ids=read_arguments_from_json_file(aprepro_vars)


cubit.cmd('reset')

for volume in structure_step_files:
  cubit.cmd('import step "'+volume+'" heal')

cubit.cmd('imprint body all')
cubit.cmd('merge tolerance 1.e-6')
cubit.cmd('merge all')


list_of_external_surfaces = find_external_surfaces()

mesh_and_remesh_till_done(quality)

nodes_in_volumes = sorted(cubit.parse_cubit_list("node"," in volume all "))
nodal_coordinates_list = []
# print('nodes_in_volumes',nodes_in_volumes)
for node_id in nodes_in_volumes:
    nodal_coordinates_list.append(cubit.get_nodal_coordinates(node_id))
print(nodal_coordinates_list)

tets_in_volumes = cubit.parse_cubit_list("tet"," in volume all ")
nodes_in_tets_list=[]
for tet_id in tets_in_volumes:
    nodes_in_tets = cubit.parse_cubit_list("node"," in tet "+str(tet_id))
    nodes_in_tets_list.append(nodes_in_tets)

f = open(outputfile, 'w') 

f.write('<?xml version="1.0"?>')
f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
f.write('<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
f.write('  <Domain>\n')
f.write('    <Grid Name="mesh" GridType="Uniform">\n')
f.write('      <Topology NumberOfElements="'+str(len(tets_in_volumes))+'" TopologyType="Tetrahedron" NodesPerElement="4">\n')
f.write('        <DataItem Dimensions="'+str(len(tets_in_volumes))+' 4" NumberType="UInt" Format="XML">\n')

string_to_write_to_file=''
for tet_id  in tets_in_volumes:
     nodes_in_tets = cubit.parse_cubit_list("node"," in tet "+str(tet_id))
     string_to_write_to_file = string_to_write_to_file+'            ' + ' '.join(str(i-1) for i in nodes_in_tets)+'\n'
f.write(string_to_write_to_file)

f.write('        </DataItem>\n')
f.write('      </Topology>\n')
f.write('    <Geometry GeometryType="XYZ">\n')
f.write('        <DataItem Dimensions="'+str(len(nodal_coordinates_list))+' 3" Format="XML">\n')

string_to_write_to_file=''
for coords in nodal_coordinates_list:
   string_to_write_to_file= string_to_write_to_file+ '            '+' '.join(str(i) for i in coords)+'\n'
f.write(string_to_write_to_file)

f.write('        </DataItem>\n')
f.write('    </Geometry>\n')



f.write('      <Attribute Name="volume_marker_material" AttributeType="Scalar" Center="Cell">\n')
f.write('      <DataItem Dimensions="'+str(len(tets_in_volumes))+' 1" Format="XML">\n')
all_volumes = cubit.parse_cubit_list("volume", ' all')
string_to_write_to_file=''
tets_in_each_volume = []
for tet_id in tets_in_volumes:
  volume_id = cubit.parse_cubit_list("volume", "in tet "+str(tet_id))[0]
  string_to_write_to_file += str(volume_id)+'\n'


f.write(string_to_write_to_file)
f.write('        </DataItem>\n')
f.write('    </Attribute>\n')


triangles_in_tets = sorted(cubit.parse_cubit_list("tri"," in tet all "))
# print('triangles_in_tets',triangles_in_tets)

f.write('    </Grid>\n')
f.write('    <Grid Name="mesh" GridType="Uniform">\n')

f.write('      <Topology NumberOfElements="'+str(len(triangles_in_tets))+'" TopologyType="Triangle" NodesPerElement="3">\n')
f.write('      <DataItem Dimensions="'+str(len(triangles_in_tets))+' 3" NumberType="UInt" Format="XML">\n')

string_to_write_to_file=''
for tri_ids in triangles_in_tets:
    nodes_in_triangles = cubit.parse_cubit_list("node"," in tri "+str(tri_ids))
    # print(nodes_in_triangles)
    string_to_write_to_file=string_to_write_to_file+'      '+' '.join(str(i-1) for i in nodes_in_triangles)+'\n'

f.write(string_to_write_to_file)
f.write('        </DataItem>\n')
f.write('    </Topology>\n')
f.write('      <Geometry Reference="XML">/Xdmf/Domain/Grid/Geometry</Geometry>\n')
f.write('      <Attribute Name="surface_marker" AttributeType="Scalar" Center="Cell">\n')
f.write('      <DataItem Dimensions="'+str(len(triangles_in_tets))+' 1" Format="XML">\n')


# print('tri_id_in_special_surfaces',tri_id_in_special_surfaces)


string_to_write_to_file = ''
      

for tri_id in triangles_in_tets:
  # print(tri_id)
  surface_id = cubit.parse_cubit_list("surface", " in tri "+str(tri_id))[0]
  string_to_write_to_file = string_to_write_to_file+'            ' + str(surface_id) + '\n'
  

f.write(string_to_write_to_file)
f.write('        </DataItem>\n')
f.write('    </Attribute>\n')
f.write('    </Grid>\n')
f.write('    </Domain>\n')
f.write('    </Xdmf>\n')
f.close()
