import json
from pprint import pprint
import ast


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
  outputfile = data["outputfile"]
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


