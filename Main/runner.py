mot = __import__('3-MOT')
data = mot.get_databases('Problems/Square_Pipe/Parameters/MOT_parameters_square_pipe.json')
solve_transient, solve_heat_transfer, solve_laminar_flow, solve_diffusion, solve_diffusion_coefficient_temperature_dependent, solve_with_decay = mot.get_solvers(data)
Time, num_steps, dt = mot.get_solving_parameters(data, solve_transient)
mesh, n0, volume_marker, dx, surface_marker, ds, mesh_fluid, dx_fluid, surface_marker_fluid, ds_fluid, n_fluid = mot.define_mesh(data, solve_laminar_flow)
V, V0, U, Q = mot.define_functionspaces(data, mesh, mesh_fluid)
bcu, bcp = mot.define_BC_laminar_flow(data, solve_laminar_flow, U, Q, surface_marker_fluid)

c_n, T_n = mot.define_initial_values(solve_heat_transfer, solve_diffusion, data, V)

Source_c_diffusion, Source_T_diffusion = mot.define_source_terms(solve_heat_transfer, solve_diffusion, dx, data, V)

bcs_c, Neumann_BC_c_diffusion, Robin_BC_c_diffusion = mot.define_BC_diffusion(data, solve_diffusion, V, surface_marker, ds)

bcs_T, Neumann_BC_T_diffusion, Robin_BC_T_diffusion = mot.define_BC_heat_transfer(data, solve_heat_transfer, V, surface_marker, ds)


D, thermal_conductivity, specific_heat, density = mot.define_materials_properties(V0, data, volume_marker, solve_heat_transfer, solve_diffusion)

F = mot.define_variational_problem_diffusion(solve_diffusion, solve_transient, dt, solve_with_decay, V, Neumann_BC_c_diffusion, Robin_BC_c_diffusion, Source_c_diffusion)

FT = mot.define_variational_problem_heat_transfer(solve_heat_transfer, solve_transient, dt, V, specific_heat, density, thermal_conductivity, Neumann_BC_T_diffusion, Robin_BC_T_diffusion, Source_T_diffusion)

A1, L1, A2, L2, A3, L3 = mot.define_variational_problem_laminar_flow(solve_laminar_flow, solve_transient, U, Q, dx_fluid, ds_fluid, n_fluid, bcu, bcp, dt)

mot.solving(data, solve_transient, solve_heat_transfer, solve_diffusion, solve_laminar_flow, solve_diffusion_coefficient_temperature_dependent, Time, num_steps, dt, V, U, Q, False, False, F, Source_c_diffusion, bcs_c, FT, Source_T_diffusion, bcs_T, A1, L1, A2, L2, A3, L3, bcu, bcp, ds, dx, volume_marker, n0)