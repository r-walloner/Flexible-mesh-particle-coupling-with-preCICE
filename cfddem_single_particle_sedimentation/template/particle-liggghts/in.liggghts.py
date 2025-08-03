def generate(p):
    return f"""
# System settings
units           si
atom_style      granular
boundary        p f p 
newton          off
processors      1 1 1
communicate     single vel yes
region          domain block {-p["particle_diameter"] * 25 / 2} {p["particle_diameter"] * 25 / 2} 0 {p["particle_diameter"] * 75} {-p["particle_diameter"] * 25 / 2} {p["particle_diameter"] * 25 / 2} units box
create_box      1 domain
neigh_modify    delay 0
timestep        {p["particle_dt"]}
fix		        integrate all nve/sphere


# Collision model
# Does not matter for single particle case
fix m1 all property/global kn                   peratomtypepair 1 5000
fix m2 all property/global kt                   peratomtypepair 1 5000
fix m3 all property/global gamman               peratomtypepair 1 0.97
fix m4 all property/global gammat               peratomtypepair 1 0.97
fix m5 all property/global coefficientFriction  peratomtypepair 1 0.1
pair_style gran model hooke/stiffness tangential history
pair_coeff * *


# Create particles
create_atoms    1 single 0 {p["particle_diameter"] * 60} 0
set type        1 diameter {p["particle_diameter"]}
set type        1 density {p["particle_density"]}


# Set up coupling
precice_initialize Particle ../precice-config.xml {2 * p["write_mapping_radius"] if p["write_mapping"] == "coarse-graining" else p["particle_diameter"]} Fluid-Mesh
compute voro all voronoi/atom
fix cpl all fluid_coupling {p["particle_drag_model"]} {p["coupling"]} {p["fluid_density"]} {p["fluid_viscosity"]} 1 1 {p["write_mapping_radius"] if p["write_mapping"] == "coarse-graining" else "0"}


# Output to console
thermo_style    custom step dt time atoms ke 
thermo_modify   lost warn norm no
thermo          200


# Output to file
variable v_fluid_x atom f_cpl[1]
variable v_fluid_y atom f_cpl[2]
variable v_fluid_z atom f_cpl[3]
variable volume atom f_cpl[4]
variable alpha_p atom f_cpl[5]
variable f_drag_x atom f_cpl[6]
variable f_drag_y atom f_cpl[7]
variable f_drag_z atom f_cpl[8]
variable expl_momentum_x atom f_cpl[9]
variable expl_momentum_y atom f_cpl[10]
variable expl_momentum_z atom f_cpl[11]
variable impl_momentum atom f_cpl[12]
variable drag_coeff atom f_cpl[13]
variable v_rel_x atom f_cpl[14]
variable v_rel_y atom f_cpl[15]
variable v_rel_z atom f_cpl[16]
variable alpha_f atom f_cpl[17]
variable reynolds atom f_cpl[18]
variable f_buoyancy_x atom f_cpl[19]
variable f_buoyancy_y atom f_cpl[20]
variable f_buoyancy_z atom f_cpl[21]
variable f_fluid_total_x atom f_cpl[22]
variable f_fluid_total_y atom f_cpl[23]
variable f_fluid_total_z atom f_cpl[24]
variable f_gravity_x atom f_cpl[25]
variable f_gravity_y atom f_cpl[26]
variable f_gravity_z atom f_cpl[27]

dump dmp all custom/vtk {int(p["output_interval"] / p["particle_dt"])} out/particles_*.vtu &
    id &
    x y z &
    ix iy iz &
    omegax omegay omegaz &
    vx vy vz &
    v_v_fluid_x v_v_fluid_y v_v_fluid_z &
    v_v_rel_x v_v_rel_y v_v_rel_z &
    fx fy fz &
    v_f_drag_x v_f_drag_y v_f_drag_z &
    v_f_gravity_x v_f_gravity_y v_f_gravity_z &
    v_f_buoyancy_x v_f_buoyancy_y v_f_buoyancy_z &
    v_f_fluid_total_x v_f_fluid_total_y v_f_fluid_total_z &
    v_expl_momentum_x v_expl_momentum_y v_expl_momentum_z &
    v_impl_momentum &
    v_drag_coeff v_reynolds &
    radius v_volume v_alpha_p v_alpha_f
{"""
dump_modify dmp binary yes
dump_modify dmp compressor lz4"""
 if p["output_compression"] else ""}


# Run with coupling
run {int(p["end_time"] / p["particle_dt"])} pre no post no every 1 precice_advance


precice_finalize
"""