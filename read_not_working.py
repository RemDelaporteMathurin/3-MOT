import argparse

from dolfin import *

import json

from pprint import pprint

parser = argparse.ArgumentParser()

# parser.add_argument("-m","--mesh", help="XDMF Mesh file input ", required=False, default='mesh_and_markers.xdmf')
# parser.add_argument("-om","--output_mesh", help="XDMF Mesh file output ", required=False, default='mesh_and_markers_from_fenics.xdmf')
# parser.add_argument("-ot","--output_temperature", help="XDMF Mesh file output with temperature", required=False, default='temperature_output.xdmf')
parser.add_argument("-j","--json_input", help="XDMF Mesh file output with temperature", required=True)

args = parser.parse_args()

def read_arguments_from_json_file(json_input):
  print('reading arguments from json file')
  print('json_inputfile ='+str(json_input))
  with open(json_input) as f:
    data = json.load(f)
  pprint(data)
  data=byteify(data)
  mesh_file = data["mesh_file"]
  field_file = data["field_file"]
  neumann_surfaces=data["neumann_surfaces"]
  structure_and_materials=data["structure_and_materials"]
  materials=[]
  material_ids=[]
  diffusivity=[]
  for entry in structure_and_materials:
    if entry['material_id'] not in material_ids:
        materials.append(entry['material'])
        material_ids.append(entry['material_id'])
        diffusivity.append(entry['diffusivity'])
  return mesh_file, field_file, material_ids, materials, diffusivity, neumann_surfaces

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

mesh_file, field_file, material_ids, materials, diffusivity, neumann_surfaces = read_arguments_from_json_file(args.json_input)

print(mesh_file, field_file, material_ids, materials, diffusivity, neumann_surfaces)



xdmf_encoding = XDMFFile.Encoding.ASCII
#xdmf_encoding = XDMFFile.Encoding.HDF5


# Read in Mesh and markers from file
mesh = Mesh()
xdmf_in = XDMFFile(MPI.comm_world, mesh_file)
xdmf_in.read(mesh)

# prepare output file for writing by writing the mesh to the file
xdmf_out = XDMFFile(MPI.comm_world, mesh_file.split('.')[1]+'_from_fenics.xdmf')
#xdmf_out.write(mesh, xdmf_encoding)

# read in volume markers from the mesh file
volume_marker_material = MeshFunction("int", mesh, mesh.topology().dim(), 0) # should change to int
xdmf_in.read(volume_marker_material, "volume_marker_material")
volume_marker_material.rename("volume_marker_material", "volume_marker_material")
xdmf_out.write(volume_marker_material, xdmf_encoding)


# read in volume markers from the mesh file
volume_marker_source = MeshFunction("double", mesh, mesh.topology().dim(), 0)
xdmf_in.read(volume_marker_source, "volume_marker_source")
volume_marker_source.rename("volume_marker_source", "volume_marker_source")
xdmf_out.write(volume_marker_source, xdmf_encoding)





# read in Surface markers for boundary conditions (first wal flux and cooling)
surface_marker_mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
xdmf_in.read(surface_marker_mvc, "surface_marker")

surface_marker_mvc.rename("surface_marker", "surface_marker")
xdmf_out.write(surface_marker_mvc, xdmf_encoding)

surface_marker = MeshFunction("size_t", mesh, surface_marker_mvc)


Q = FunctionSpace(mesh, "CG", 1) #change to P
u = TrialFunction(Q)
v = TestFunction(Q)



u0 = Constant(500) # temperature of coolant in kelvin 
htc = Constant(3500) # W per m2 Kelvin ,heat transfer coeff 1.5 is for natural convection, this is to low and would need fixing
q0 = Constant(0.5e6) # heat flux value fusion first wall = 0.5 MW/m2
u1 = Constant(400) #fixed temperature for dirich BC in kelvin, fixed temperature boundary for back wall


#assumes an order of for the volume markers, needs to be updated to find the order of materials in the volume markers
# int values are allowed for volume markers in xdmf format

#for material_id, material, diffusivity in material_ids, materials, diffusivities
#    diffusivity_dict[material_id]=(material,diffusivity)


diffusivity_dict = {1: ('tungsten',60e-6),
                    2: ('lithium_lead',4e-5),
                    3: ('eurofer', 1.8e-5),
                    }

print(diffusivity_dict)

# diffusivity = k / (pc)  
# units of diffusivity are m2 per second
# k is thermal conductivity (W/(mK))
# p is density (kg/m3)
# c is specific heat capacity (J/(kgK))




class Assign_volume_markers_with_dictionary_values(UserExpression): # used to assign diffusivity
    def __init__(self, meshfn, diffusivity_dict):
        self.meshfn = meshfn
        self.diff_dict = diffusivity_dict
        super().__init__()
    def eval_cell(self, values, x, ufc_cell):
        #try:
            cell = Cell(self.meshfn.mesh(), ufc_cell.index)
            values[0] = self.diff_dict[self.meshfn[cell]][1]
        #except: 
           # print(values, x, ufc_cell,self.meshfn[cell])
            #input()



diffusivity = Assign_volume_markers_with_dictionary_values(volume_marker_material, diffusivity_dict)


class MySimpleExpression(UserExpression): # used to assign
    def __init__(self, meshfn):
        self.meshfn = meshfn
        super().__init__()
    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.meshfn.mesh(), ufc_cell.index)
        values[0] = self.meshfn[cell]*1e-6 #scalled from W/cm3 to W/m3
        #print( self.meshfn[cell])

volumetric_heat_source = MySimpleExpression(volume_marker_source)


ds = Measure('ds', domain=mesh, subdomain_data = surface_marker)




#convective boundary condition
#robins bc
surface_cooling_fw_int_id = 13
print('cooling_surface area =',assemble(1*ds(surface_cooling_fw_int_id)))
coolant_fw_bc = htc*(u - u0) * v *  ds(surface_cooling_fw_int_id)


surface_flux_int_id = 24
print('heating_surface area =',assemble(1*ds(surface_flux_int_id)))
flux_bc = q0 * v * ds(surface_flux_int_id)

rear_wall_int_id = 36
#rear_wall_bc = 
back_wall = ds()


def k(T):
   return conditional(lt(T, 0.5), 1.0, 4*T - 1)



#eq = dot(diffusivity*grad(u), grad(v))*dx  - diffusivity*coolant_fw_bc + diffusivity*flux_bc + volumetric_heat_source*v*dx 
#eq = dot(diffusivity*grad(u), grad(v))*dx  - 1e-5*coolant_fw_bc + 1e-5*flux_bc #+ volumetric_heat_source*v*dx 


#eq = diffusivity * dot(grad(u), grad(v))*dx  - 1e-5*(-coolant_fw_bc + flux_bc) - volumetric_heat_source*v*dx 
eq = diffusivity * dot(grad(u), grad(v))*dx  - 1e-5*(-coolant_fw_bc) - 1*v*dx 



#eq = dot(diffusivity_better(T)*grad(u), grad(v))*dx  + coolant_fw_bc+ flux_bc +volumetric_heat_source*v*dx

a = lhs(eq)
L = rhs(eq)



# set the fixed temperature boundary condition using DirichletBC
#bcs = [DirichletBC(Q, u1, surface_marker, 0)]
bcs=[]


temperature = Function(Q)
solve(a == L, temperature, bcs)


temperature.rename("temperature", "temperature")
xdmf_out2 = XDMFFile(MPI.comm_world, field_file)
xdmf_out2.write(temperature, xdmf_encoding)


# for cell_no in diffusivity:
#     print(cell_no)
# print(max temperature for each different material)