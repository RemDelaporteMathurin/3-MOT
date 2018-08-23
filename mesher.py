import json

# run this script with the following commands in order to mesh from step files
# without the Trelis GUI
# trelis -nographics -batch mesher.py "json_input='input_parameters/breeder_blankets_parameters.json'"
# with the Trelis GUI
# trelis mesher.py "json_input='input_parameters/breeder_blankets_parameters.json'"


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


def get_nodes_in_tets_volume_id_string(volumes):
    string_ids = '          '
    string_nodes = '          '
    total_tets = 0
    for vol in volumes:
        tets_in_volume = 0
        tets = cubit.parse_cubit_list("tet", "in vol " + str(vol))
        #print('Writing tets in volume ' + str(vol))
        for tet in tets:
            tets_in_volume += 1
            nodes = cubit.parse_cubit_list("node", "in tet " + str(tet))
            for node in nodes:
                string_nodes += str(node - 1) + ' '
            string_nodes += '\n          '
            string_ids += str(vol) + '\n          '
        print('Found ' + str(tets_in_volume) + ' tets in volume ' + str(vol))
        total_tets += tets_in_volume
    print('Found ' + str(total_tets) + ' tets in total')
    return string_nodes, string_ids


def get_nodes_coordinates_string(nodes):
    string = ''
    for node in nodes:
        string = string + '              '
        coord = cubit.get_nodal_coordinates(node)
        string += str(coord[0]) + ' '
        string += str(coord[1]) + ' '
        string += str(coord[2]) + '\n'
    return string


def get_nodes_in_tris_surface_id_strings(surfaces):
    string_nodes = '          '
    string_ids = '          '
    for surf in surfaces:
        entity = "tri"
        cmd = " in surf " + str(surf)
        triangles = cubit.parse_cubit_list(entity, cmd)
        for tri in triangles:
            entity = "node"
            cmd = " in tri " + str(tri)
            nodes = cubit.parse_cubit_list(entity, cmd)
            string_nodes += ''
            for node in nodes:
                string_nodes += str(node - 1) + ' '
            string_nodes += '\n          '
            string_ids+=str(surf)+'\n          '
    return string_nodes, string_ids


def write_file(data):
    tets_in_volumes = cubit.parse_cubit_list("tet", " in volume all ")
    triangles_in_tets = sorted(cubit.parse_cubit_list("tri", " in tet all "))
    volumes = cubit.parse_cubit_list("volume", "in vol all")
    surfaces = cubit.parse_cubit_list("surface", "in surface all")
    string_tets_nodes, string_tets_ids = get_nodes_in_tets_volume_id_string(volumes)
    print('Tets connectivity : Done')
    print('Tets ids : Done')
    string_tris_nodes, string_tris_ids = get_nodes_in_tris_surface_id_strings(surfaces)
    print('Tris connectivity : Done')
    print('Tris ids : Done')
    #print(string_tris_ids)
    f = open(data["mesh_file"], 'w')
    f.write('<?xml version="1.0"?>')
    f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
    f.write('<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
    f.write('  <Domain>\n')
    f.write('    <Grid Name="mesh" GridType="Uniform">\n')
    f.write('      <Topology NumberOfElements="' + str(len(tets_in_volumes)) + '" TopologyType="Tetrahedron" NodesPerElement="4">\n')
    f.write('        <DataItem Dimensions="' + str(len(tets_in_volumes)) + ' 4" NumberType="UInt" Format="XML">\n')
    f.write(string_tets_nodes)
    f.write('</DataItem>\n')
    f.write('      </Topology>\n')
    f.write('      <Geometry GeometryType="XYZ">\n')
    f.write('        <DataItem Dimensions="' + str(len(cubit.parse_cubit_list("node", " in node all ")))+' 3" Format="XML">\n')
    f.write(get_nodes_coordinates_string(sorted(cubit.parse_cubit_list("node", " in node all "))))
    f.write('        </DataItem>\n')
    f.write('      </Geometry>\n')
    f.write('      <Attribute Name="volume_marker_volume_id" AttributeType="Scalar" Center="Cell">\n')
    f.write('<DataItem Dimensions="' + str(len(tets_in_volumes)) + ' 1" Format="XML">\n')
    f.write(string_tets_ids)
    f.write('        </DataItem>\n')
    f.write('    </Attribute>\n')
    f.write('    </Grid>\n')
    f.write('    <Grid Name="mesh" GridType="Uniform">\n')
    f.write('      <Topology NumberOfElements="' + str(len(triangles_in_tets)) + '" TopologyType="Triangle" NodesPerElement="3">\n')
    f.write('        <DataItem Dimensions="' + str(len(triangles_in_tets)) + ' 3" NumberType="UInt" Format="XML">\n')
    f.write(string_tris_nodes)
    f.write('</DataItem>\n')
    f.write('      </Topology>\n')
    f.write('      <Geometry Reference="XML">/Xdmf/Domain/Grid/Geometry</Geometry>\n')
    f.write('      <Attribute Name="surface_marker" AttributeType="Scalar" Center="Cell">\n')
    f.write('        <DataItem Dimensions="' + str(len(triangles_in_tets)) + ' 1" Format="XML">\n')
    f.write(string_tris_ids)
    f.write('</DataItem>\n')
    f.write('      </Attribute>\n')
    f.write('    </Grid>\n')
    f.write('  </Domain>\n')
    f.write('</Xdmf>\n')
    f.close()
    print('Cest fini')




data = get_json_input(cubit.get_aprepro_vars())
project_files = data['structure_and_materials']['project_files']
if type(project_files)==list and all(file.endswith('.step') or file.endswith('.stp') for file in project_files):
    auto_mesh(data)
    write_file(data)
elif project_files.endswith('.cub'):
    print('Opening from project', project_files)
    cubit.cmd('open "' + str(project_files) + '"')
    write_file(data)


#trailing empty line required by Trelis python 