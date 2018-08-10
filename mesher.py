import json
from pprint import pprint
import ast
import argparse

# run this script with the following commands in order to mesh from step files
# trelis -nographics -batch mesher.py "json_input='MOT_parameters_breeder_blankets.json'"
# trelis mesher.py "json_input='MOT_parameters_breeder_blankets.json'"


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

def get_json_input(aprepro_vars):
    if "json_input" in aprepro_vars:
        json_input = str(cubit.get_aprepro_value_as_string("json_input"))
        with open(json_input) as f:
            data = json.load(f)
        data = byteify(data)
        return data
    else:
        raise ValueError('provide a json_input file')

def auto_mesh(data):
    print('Meshing from STEP files')
    cubit.cmd('reset')
    for step_file in data['structure_and_materials']['project_files']:
        print('loading step file ', step_file)
        cubit.cmd('import step "' + step_file + '" heal')
    cubit.cmd('vol all scale ' + str(data['scaling']))
    cubit.cmd('imprint body all')
    cubit.cmd('merge tolerance 1.e-6')
    cubit.cmd('merge all')
    cubit.cmd('volume all scheme Tetmesh')
    cubit.cmd('volume all scheme tetmesh proximity layers off geometry approximation angle 15 geometric sizing on ')
    cubit.cmd('volume all size auto factor ' + str(data['quality'])) # replace with sizing function 10 times smaller than smallest size in model
    cubit.cmd('volume all tetmesh growth_factor 1 ')
    cubit.cmd('Trimesher surface gradation 1.3')
    cubit.cmd('Trimesher volume gradation 1.3')
    cubit.cmd('mesh volume all')
    return

def get_nodes_in_tets_string(tets_in_volumes):
    string = ''
    for tet_id in tets_in_volumes:
        nodes_in_tet = cubit.parse_cubit_list("node", " in tet " + str(tet_id))
        string = string + '            ' + ' '.join(str(i - 1) for i in nodes_in_tet)+ '\n'
    return string

def get_nodes_coordinates_string(nodes_in_volumes):
    string = ''
    for node in nodes_in_volumes:
        string = string + '            ' + ' '.join(str(i) for i in cubit.get_nodal_coordinates(node))+'\n'
    return string

def get_volumie_id_string(tets_in_volumes):
    string = ''
    for tet_id in tets_in_volumes:
        volume_id = cubit.parse_cubit_list("volume", "in tet " + str(tet_id))[0]
        string += str(volume_id) + '\n'
    return string

def get_nodes_in_tris_string(triangles_in_tets):
    string = ''
    for tri_ids in triangles_in_tets:
        nodes_in_triangles = cubit.parse_cubit_list("node", " in tri " + str(tri_ids))
        string += '      ' + ' '.join(str(i - 1) for i in nodes_in_triangles) + '\n'
    return string

def get_surface_id_string(triangles_in_tets):
    string = ''
    for tri_id in triangles_in_tets:
        surface_id = cubit.parse_cubit_list("surface", " in tri " + str(tri_id))[0]
        string += '            ' + str(surface_id) + '\n'
    return string

def write_file(data):
    tets_in_volumes = cubit.parse_cubit_list("tet", " in volume all ")
    triangles_in_tets = sorted(cubit.parse_cubit_list("tri", " in tet all "))
    f = open(data["mesh_file"], 'w')
    f.write('<?xml version="1.0"?>')
    f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
    f.write('<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
    f.write('  <Domain>\n')
    f.write('    <Grid Name="mesh" GridType="Uniform">\n')
    f.write('      <Topology NumberOfElements="' + str(len(tets_in_volumes)) + '" TopologyType="Tetrahedron" NodesPerElement="4">\n')
    f.write('        <DataItem Dimensions="' + str(len(tets_in_volumes)) + ' 4" NumberType="UInt" Format="XML">\n')
    f.write(get_nodes_in_tets_string(tets_in_volumes))
    f.write('        </DataItem>\n')
    f.write('      </Topology>\n')
    f.write('      <Geometry GeometryType="XYZ">\n')
    f.write('        <DataItem Dimensions="' + str(len(cubit.parse_cubit_list("node", " in volume all ")))+' 3" Format="XML">\n')
    f.write(get_nodes_coordinates_string(sorted(cubit.parse_cubit_list("node", " in volume all "))))
    f.write('        </DataItem>\n')
    f.write('      </Geometry>\n')
    f.write('      <Attribute Name="volume_marker_volume_id" AttributeType="Scalar" Center="Cell">\n')
    f.write('        <DataItem Dimensions="' + str(len(tets_in_volumes)) + ' 1" Format="XML">\n')
    f.write(get_volumie_id_string(tets_in_volumes))
    f.write('        </DataItem>\n')
    f.write('    </Attribute>\n')
    f.write('    </Grid>\n')
    f.write('    <Grid Name="mesh" GridType="Uniform">\n')
    f.write('      <Topology NumberOfElements="' + str(len(triangles_in_tets)) + '" TopologyType="Triangle" NodesPerElement="3">\n')
    f.write('        <DataItem Dimensions="' + str(len(triangles_in_tets)) + ' 3" NumberType="UInt" Format="XML">\n')
    f.write(get_nodes_in_tris_string(triangles_in_tets))
    f.write('        </DataItem>\n')
    f.write('      </Topology>\n')
    f.write('      <Geometry Reference="XML">/Xdmf/Domain/Grid/Geometry</Geometry>\n')
    f.write('      <Attribute Name="surface_marker" AttributeType="Scalar" Center="Cell">\n')
    f.write('        <DataItem Dimensions="' + str(len(triangles_in_tets)) + ' 1" Format="XML">\n')
    f.write(get_surface_id_string(triangles_in_tets))
    f.write('        </DataItem>\n')
    f.write('      </Attribute>\n')
    f.write('    </Grid>\n')
    f.write('  </Domain>\n')
    f.write('</Xdmf>\n')
    f.close()
    print('Cest fini')


data = get_json_input(cubit.get_aprepro_vars())
if type(data['structure_and_materials']['project_files'])==list and (data['structure_and_materials']['project_files'][0].endswith('.step') or data['structure_and_materials']['project_files'][0].endswith('.stp')):
    auto_mesh(data)
    write_file(data)
else:
    if data['structure_and_materials']['project_files'].endswith('.cub'):
        print('Opening from project', data['structure_and_materials']['project_files'])
        cubit.cmd('open "' + str(data['structure_and_materials']['project_files']) + '"')
        write_file(data)


#trailing empty line required by Trelis python 