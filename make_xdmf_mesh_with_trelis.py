import json
from pprint import pprint
import ast

# run this script with the following commands
# trelis -nographics -batch make_xdmf_mesh_with_trelis.py "quality='10'" "outputfile='mesh_and_markers.xdmf'" "structure='slice_armour_mod1.step,slice_back_lithium_lead_mod1.step,slice_back_plate_1_mod1.step,slice_back_plate_2_mod1.step,slice_back_plate_3_mod1.step,slice_cooling_plate_material_mod1.step,slice_first_wall_material_mod1.step,slice_lithium_lead_mod1.step,'" "coolant='slice_first_wall_coolant_mod1.step,slice_cooling_plate_coolant_mod1.step,slice_back_helium_mod1.step'"
# trelis make_xdmf_mesh_with_trelis.py "quality='10'" "outputfile='mesh_and_markers.xdmf'" "structure='slice_armour_mod1.step,slice_back_lithium_lead_mod1.step,slice_back_plate_1_mod1.step,slice_back_plate_2_mod1.step,slice_back_plate_3_mod1.step,slice_cooling_plate_material_mod1.step,slice_first_wall_material_mod1.step,slice_lithium_lead_mod1.step,'" "coolant='slice_first_wall_coolant_mod1.step,slice_cooling_plate_coolant_mod1.step,slice_back_helium_mod1.step'"

# trelis make_xdmf_mesh_with_trelis.py "json_input='name_of_the_json_file'"


def find_number_of_volumes_in_each_step_file(input_locations):
    body_ids=''
    volumes_in_each_step_file=[]
    for i in range(0,len(input_locations)):
      current_vols =cubit.parse_cubit_list("volume", "all")
      print('input       ',    input_locations[i])
      if input_locations[i].endswith('.sat'):
        cubit.cmd('import acis "'+input_locations[i]+'" nofreesurfaces separate_bodies')
      if input_locations[i].endswith('.stp') or input_locations[i].endswith('.step'):
        cubit.cmd('import step "'+input_locations[i]+'" heal')
      #body_ids=body_ids+' '+str(i+1)
      all_vols =cubit.parse_cubit_list("volume", "all")
      new_vols = set(current_vols).symmetric_difference(set(all_vols))
      print('all_vols    ',str(all_vols))
      print('current_vols',str(current_vols))
      print('new_vols    ',str(new_vols))
      #volumes_in_each_step_file.append(new_vols)
      new_vols=map(str, new_vols)
      new_vols=' '.join(new_vols)
      volumes_in_each_step_file.append(new_vols.split())
    print('volumes_in_each_step_file')
    print(volumes_in_each_step_file)
    print('body_ids')
    print(body_ids)
    return volumes_in_each_step_file

def find_shared_surfaces_between_coolant_and_structure(volumes_in_coolant_step_file):
  print('looking for merged surfaces')
  surfaces_in_all_volumes = cubit.parse_cubit_list("surface"," in volume all ")
  surfaces_in_all_coolant_volumes = cubit.parse_cubit_list("surface"," in volume "+str(' '.join([item for sublist in volumes_in_coolant_step_file for item in sublist])))
  print('surfaces_in_all_coolant_volumes',surfaces_in_all_coolant_volumes)
  list_of_merged_surfaces =[]
  for surface in surfaces_in_all_volumes:
      is_merged = cubit.is_merged("surface", surface)
      #print(surface, is_merged)
      if is_merged == False:
          cubit.cmd('color surface '+str(surface)+ ' lightslateblue')
      if is_merged == True and surface in surfaces_in_all_coolant_volumes:
          cubit.cmd('color surface '+str(surface)+ ' orchid')
          # print('cooling surface id ',surface)
          list_of_merged_surfaces.append(surface)
          cubit.cmd('highlight surface '+str(surface))  
  return list_of_merged_surfaces

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
  pprint(data)
  data=byteify(data)
  quality= data["quality"]
  outputfile = data["mesh_file"]
  coolant_step_files= data["void"] #todo rename coolant by void
  neumann_surfaces=data["neumann_surfaces"]
  structure_and_materials=data["structure_and_materials"]
  structure_step_files=[]
  materials=[]
  material_ids=[]
  for entry in structure_and_materials:
    structure_step_files.append(entry['step_file'])
    materials.append(entry['material'])
    material_ids.append(entry['material_id'])
  return quality, outputfile, coolant_step_files, structure_step_files, materials, material_ids, neumann_surfaces

def read_arguments_from_terminal_input(aprepro_vars):
  print('reading arguments from terminal')
  if "quality" in aprepro_vars:
    quality = str(cubit.get_aprepro_value_as_string("quality"))
    print('quality ='+str(quality))  
  else:
    quality = '10'
  if "outputfile" in aprepro_vars:
    outputfile = str(cubit.get_aprepro_value_as_string("outputfile"))
    print('outputfile ='+str(outputfile))  
  else:
    outputfile = 'mesh_and_markers.xdmf'
  if "coolant" in aprepro_vars:
    coolant_step_files = cubit.get_aprepro_value_as_string("coolant").split(',')
    print('coolant geometry file ='+str(coolant_step_files))
  if "neumann_surfaces" in aprepro_vars:
    neumann_surfaces = cubit.get_aprepro_value_as_string("neumann_surfaces").split(',')
    print('coolant geometry file ='+str(coolant_step_files))
  if "structure_and_materials" in aprepro_vars:
    structure_and_materials = cubit.get_aprepro_value_as_string("structure").split(',')
    structure_and_materials_list_of_dicts = ast.literal_eval(structure_and_materials)
    structure_step_files=[]
    materials=[]
    material_ids=[]
    for entry in structure_and_materials:
      structure_step_files.append(entry['step_file'])
      materials.append(entry['material'])
      material_ids.append(entry['material_id'])
    print('structure geometry files ='+str(structure_step_files))
    print('number of structure geometry files ='+str(len(structure_step_files)))
    print('materials =', materials)
    print('material id numbers =', material_ids)
  return quality, outputfile, coolant_step_files, structure_step_files, materials, material_ids, neumann_surfaces

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
  print('list_of_external_surfaces',list_of_external_surfaces)
  return list_of_unmerged_surfaces





aprepro_vars = cubit.get_aprepro_vars()

print("Found the following aprepro variables:")
print(aprepro_vars)
for var_name in aprepro_vars:
  val = cubit.get_aprepro_value_as_string(var_name)
  print("{0} = {1}".format(var_name, val))

if "json_input" in aprepro_vars:
  quality, outputfile, coolant_step_files, structure_step_files, materials, material_ids , neumann_surfaces= read_arguments_from_json_file(aprepro_vars)
else:
  quality, outputfile, coolant_step_files, structure_step_files, materials, material_ids, neumann_surfaces= read_arguments_from_terminal_input(aprepro_vars)



print("quality",quality)
print("outputfile",outputfile)
print("coolant_step_files",coolant_step_files)
print("structure_step_files",structure_step_files)
print("materials",materials)
print("material_ids",material_ids)




cubit.cmd('reset')

volumes_in_structure_step_files = find_number_of_volumes_in_each_step_file(structure_step_files)
print('volumes_in_each_step_file',volumes_in_structure_step_files)



f1 = open('structure_step_file_and_volumes.txt', 'w') 
for structure_step_file,volumes_in_structure_step_file in zip(structure_step_files , volumes_in_structure_step_files):
  f1.write(structure_step_file +' ' + ' '.join(str(i) for i in volumes_in_structure_step_file) +'\n')
f1.close()

volumes_in_coolant_step_file = find_number_of_volumes_in_each_step_file(coolant_step_files)
print('volumes_in_each_step_file',volumes_in_coolant_step_file)

cubit.cmd('imprint body all')
cubit.cmd('merge tolerance 1.e-6')
cubit.cmd('merge all')


list_of_external_surfaces = find_external_surfaces()


list_of_merged_surfaces = find_shared_surfaces_between_coolant_and_structure(volumes_in_coolant_step_file)


cubit.cmd('delete volume '+str(' '.join([item for sublist in volumes_in_coolant_step_file for item in sublist])))


def mesh_and_remesh_till_done(quality):
  cubit.cmd('volume all scheme Tetmesh')
  cubit.cmd('volume all scheme tetmesh proximity layers off geometry approximation angle 15 geometric sizing on ')
  cubit.cmd('volume all size auto factor ' + str(quality)) # replace with sizing function 10 times smaller than smallest size in model
  cubit.cmd('volume all tetmesh growth_factor 1 ')
  cubit.cmd('Trimesher surface gradation 1.3')
  cubit.cmd('Trimesher volume gradation 1.3')
  cubit.cmd('mesh volume all')
  for step_file in volumes_in_structure_step_files:
    for volume in step_file:
      #print(volume)
      #print(cubit.is_meshed("volume", int(volume))) 
      if cubit.is_meshed("volume", int(volume)) == False:
        print(volume)
        print(cubit.is_meshed("volume", int(volume)))         
        print('full geometry ot meshed, try changing the quality')
        sys.exit()

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

tets_in_each_volume = []
for volume in all_volumes:
  tets_in_each_volume.append(cubit.get_volume_tets(volume))

volume_id_to_material_id_dict={}
for volume_ids, material_id in zip(volumes_in_structure_step_files,material_ids):
  for volume_id in volume_ids:
    volume_id_to_material_id_dict[int(volume_id)]=material_id

print('volume_id_to_material_id_dict',volume_id_to_material_id_dict)




string_to_write_to_file=''
for tet_id in tets_in_volumes:
    for tets_in_the_volume, volume_id in zip(tets_in_each_volume, all_volumes):#, all_volumes):
        if tet_id in tets_in_the_volume:
            string_to_write_to_file = string_to_write_to_file + '            '+str(volume_id_to_material_id_dict[volume_id])+'\n'
            break

#  f.write('            '+str(tet_id)+' '+str(volume)+'\n')    
f.write(string_to_write_to_file)
f.write('        </DataItem>\n')
f.write('    </Attribute>\n')




f.write('      <Attribute Name="volume_marker_source" AttributeType="Scalar" Center="Cell">\n')
f.write('      <DataItem Dimensions="'+str(len(tets_in_volumes))+' 1" Format="XML">\n')
all_volumes = cubit.parse_cubit_list("volume", ' all')

string_to_write_to_file=''
for tet_id  in tets_in_volumes:
    # print(tet_id)
    centre_of_mass = (cubit.get_center_point('tet',tet_id))
    string_to_write_to_file=string_to_write_to_file+'            '+str(abs(centre_of_mass[2])/1000)+'\n'
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


tri_id_in_cooling_surfaces = cubit.parse_cubit_list("tri"," in surface "+' '.join(str(i) for i in list_of_merged_surfaces))
identifyer_for_cooling_surfaces = 13

for surface_id in neumann_surfaces:
  cubit.cmd('color surface ' + str(surface_id) + ' green')

tri_id_in_heat_flux_surfaces = cubit.parse_cubit_list("tri"," in surface "+' '.join(str(i) for i in neumann_surfaces))
identifyer_for_heat_flux_surfaces = 24
# print('tri_id_in_special_surfaces',tri_id_in_special_surfaces)


string_to_write_to_file = ''
for tri_id in triangles_in_tets:
     # print(tri_id)
     if tri_id in tri_id_in_cooling_surfaces:
        string_to_write_to_file = string_to_write_to_file+'            ' + str(identifyer_for_cooling_surfaces) + '\n'
     elif tri_id in tri_id_in_heat_flux_surfaces:
        string_to_write_to_file = string_to_write_to_file+'            ' + str(identifyer_for_heat_flux_surfaces) + '\n'
     else:
        string_to_write_to_file = string_to_write_to_file+'            ' + str(0) + '\n'        

f.write(string_to_write_to_file)
f.write('        </DataItem>\n')
f.write('    </Attribute>\n')
f.write('    </Grid>\n')
f.write('    </Domain>\n')
f.write('    </Xdmf>\n')
f.close()
