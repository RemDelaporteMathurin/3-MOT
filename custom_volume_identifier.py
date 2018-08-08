from datetime import datetime
from pprint import pprint
import ast
import os
import json
import xml.etree.ElementTree as ET

# run this script with the following commands
# trelis -nographics -batch custom_volume_identifier.py "json_input='MOT_parameters_breeder_blankets.sjon'" 
# trelis custom_volume_identifier.py "json_input='MOT_parameters_breeder_blankets.json'"

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
        new_vols = map(str, new_vols)
        new_vols = ' '.join(new_vols)
        volumes_in_each_step_file.append(new_vols.split())
        volumes_in_each_step_file_list_of_dicts.append({input_locations[i]:new_vols.split()})
    print('volumes_in_each_step_file', volumes_in_each_step_file)
    print('body_ids', body_ids)
    print('volumes_in_each_step_file_list_of_dicts', volumes_in_each_step_file_list_of_dicts)
    return volumes_in_each_step_file_list_of_dicts



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



aprepro_vars = cubit.get_aprepro_vars()

if "json_input" in aprepro_vars:
    json_input = cubit.get_aprepro_value_as_string("json_input")
    print('json_input found',json_input)
    with open(json_input) as f:
        data =byteify(json.load(f))
else:
    print('provide json_input filename')
    sys.exit()


volumes_in_step_files = find_number_of_volumes_in_each_step_file(data["structure_and_materials"]['step_files'])

with open('volumes_in_step_file.json', 'w') as outfile:
    json.dump(volumes_in_step_files,outfile,indent=4)


print('script exectution time =',datetime.now() - startTime)
print('script exectution time =',datetime.now() - startTime)

print(data["structure_and_materials"]['step_files'])


