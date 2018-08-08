from datetime import datetime
from pprint import pprint
import ast
import os
import json
import xml.etree.ElementTree as ET

# run this script with the following commands
# trelis -nographics -batch custom_shared_surface_identifier.py "json_input='MOT_parameters_breeder_blankets.sjon'" "void_input='voids.json'"
# trelis custom_shared_surface_identifier.py "json_input='MOT_parameters_breeder_blankets.json'" "void_input='voids.json'"

startTime = datetime.now()

def find_number_of_volumes_in_each_step_file(input_locations):
    body_ids=''
    volumes_in_each_step_file=[]
    volumes_in_each_step_file_list_of_dicts=[]
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
      volumes_in_each_step_file_list_of_dicts.append({input_locations[i]:new_vols.split()})
    print('volumes_in_each_step_file',volumes_in_each_step_file)
    print('body_ids',body_ids)
    print('volumes_in_each_step_file_list_of_dicts',volumes_in_each_step_file_list_of_dicts)
    return volumes_in_each_step_file_list_of_dicts

def find_shared_surfaces_between_coolant_and_structure(volumes_in_coolant_step_files):
    list_of_merged_surfaces=[]
    print(volumes_in_coolant_step_file)
    for entry in volumes_in_coolant_step_files:
        for key, value in entry.items():
            filename = key
            print('filename',filename)
            volumes_in_each_step_file = value
            print('volumes_in_each_step_file',volumes_in_each_step_file)
        surfaces_in_all_coolant_volumes = cubit.parse_cubit_list("surface"," in volume "+str(' '.join(volumes_in_each_step_file)))
        print('surfaces_in_all_coolant_volumes',surfaces_in_all_coolant_volumes)
        merged_surfaces_in_step_file_dict={'void':filename}
        list_of_shared_surfaces=[]
        for surface in surfaces_in_all_coolant_volumes:
            if cubit.is_merged("surface", surface) == True:
                list_of_shared_surfaces.append(surface)
                cubit.cmd('highlight surface '+str(surface))
        print('list_of_shared_surfaces',list_of_shared_surfaces)
        merged_surfaces_in_step_file_dict['shared_surfaces']=list_of_shared_surfaces
        list_of_merged_surfaces.append(merged_surfaces_in_step_file_dict)
    print('list_of_merged_surfaces',list_of_merged_surfaces)
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


def find_external_surfaces():
    print('looking for merged surfaces')
    surfaces_in_all_volumes = cubit.parse_cubit_list("surface", " in volume all ")
    list_of_merged_surfaces = []
    list_of_unmerged_surfaces = []
    for surface in surfaces_in_all_volumes:
        if cubit.is_merged("surface", surface) == True:
            list_of_merged_surfaces.append(surface)
        else == False:
            cubit.cmd('color surface '+str(surface)+ ' red')
            list_of_unmerged_surfaces.append(surface)
    print('list_of_external_surfaces',list_of_unmerged_surfaces)
    return {'external_surfaces':list_of_unmerged_surfaces}



aprepro_vars = cubit.get_aprepro_vars()

if "json_input" in aprepro_vars:
    json_input = cubit.get_aprepro_value_as_string("json_input")
    print('json_input found',json_input)
    with open(json_input) as f:
        data =byteify(json.load(f))
if "void_input" in aprepro_vars:
    void_input = cubit.get_aprepro_value_as_string("void_input")
    print('void_input found',void_input)
    with open(void_input) as f:
        void_data = byteify(json.load(f))
else:
    print('provide both json_input and void_input files')
    sys.exit()


volumes_in_structural_step_files = find_number_of_volumes_in_each_step_file(data["structure_and_materials"]['step_files'])

with open('volumes_in_structural_step_file.json', 'w') as outfile:
    json.dump(volumes_in_structural_step_files,outfile)


volumes_in_coolant_step_file = find_number_of_volumes_in_each_step_file(void_data['void'])

with open('volumes_in_coolant_step_file.json', 'w') as outfile:
    json.dump(volumes_in_coolant_step_file,outfile)

cubit.cmd('vol all scale '+ str(data['scaling']))
cubit.cmd('imprint body all')
cubit.cmd('merge tolerance 1.e-9')
cubit.cmd('merge all')


external_surfaces = find_external_surfaces()

with open('external_surfaces.json', 'w') as outfile:
    json.dump(external_surfaces, outfile)

print(volumes_in_coolant_step_file)

list_of_merged_surfaces = find_shared_surfaces_between_coolant_and_structure(volumes_in_coolant_step_file)

with open('void_structure_shared_surfaces.json', 'w') as outfile:
    json.dump(list_of_merged_surfaces,outfile,indent=4)


print('script exectution time =',datetime.now() - startTime)

