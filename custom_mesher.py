#!python
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

f = open("custom_mesh.xdmf", 'w') 

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