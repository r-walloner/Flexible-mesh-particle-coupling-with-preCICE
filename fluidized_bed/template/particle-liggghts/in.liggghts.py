def generate(p):
    return f"""
# System settings
units           si
atom_style      granular
boundary        p f p 
newton          off
processors      {p["particle_subdomains"]}
communicate     single vel yes
region          domain block -.075 .075 0 .75 -.0075 .0075 units box
create_box      1 domain
neigh_modify    delay 0
timestep        {p["particle_dt"]}
fix		        integrate all nve/sphere


# Domain walls
fix wall_x1 all wall/gran {p["particle_contact_model"]} primitive type 1 xplane -.075
fix wall_x2 all wall/gran {p["particle_contact_model"]} primitive type 1 xplane +.075
fix wall_y1 all wall/gran {p["particle_contact_model"]} primitive type 1 yplane -0
fix wall_y2 all wall/gran {p["particle_contact_model"]} primitive type 1 yplane +.75
fix wall_z1 all wall/gran {p["particle_contact_model"]} primitive type 1 zplane -.0075
fix wall_z2 all wall/gran {p["particle_contact_model"]} primitive type 1 zplane +.0075


# Collision model
fix m1 all property/global youngsModulus            peratomtype         {p["particle_youngs_modulus"]}
fix m2 all property/global poissonsRatio            peratomtype         {p["particle_poissons_ratio"]}
fix m3 all property/global coefficientRestitution   peratomtypepair 1   {p["particle_restitution"]}
fix m4 all property/global coefficientFriction      peratomtypepair 1   {p["particle_friction"]}
fix m5 all property/global characteristicVelocity   scalar              {p["particle_characteristic_velocity"]}
pair_style gran {p["particle_contact_model"]}
pair_coeff * *


# Output to console
thermo_style    custom step dt time atoms ke 
thermo_modify   lost warn norm no
thermo          200


# Create particles
region insert_region block -.075 .075 0 .75 -.0075 .0075 units box
fix p_template all particletemplate/sphere 15485863 &
    atom_type 1 &
    density constant {p["particle_density"]} &
    radius constant {p["particle_diameter"] / 2}
fix p_distribution all particledistribution/discrete 15485867 1 p_template 1.0
fix insert all insert/pack seed 32452843 distributiontemplate p_distribution &
    vel constant {p["particle_insert_velocity"][0]} {p["particle_insert_velocity"][1]} {p["particle_insert_velocity"][2]} &
	insert_every once &
    overlapcheck yes &
    all_in yes &
    particles_in_region {p["particle_count"]} &
    region insert_region
run 1
unfix insert


# Let particles settle
dump dmp_settle all custom/vtk {int(p["output_interval"] / p["particle_dt"])} out/settle/particles_*.vtu &
    id &
    x y z &
    ix iy iz &
    omegax omegay omegaz &
    vx vy vz &
    fx fy fz &
    radius
dump_modify dmp_settle binary yes
dump_modify dmp_settle compressor lz4

fix settling_gravity all gravity {p["particle_settling_gravity"]} vector 0 -1 0

run {int(p["particle_settling_time"] / p["particle_dt"])}

unfix settling_gravity
undump dmp_settle


# Set up coupling
precice_initialize Particle ../precice-config.xml Fluid-Mesh
compute voro all voronoi/atom
fix cpl all fluid_coupling zhao_shan {"force" if p["solver"] == "AndersonJacksonFoam" else "momentum_semi_implicit"} {p["fluid_density"]} {p["fluid_viscosity"]} 1 1


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
dump_modify dmp binary yes
dump_modify dmp compressor lz4


# Run with coupling
run {int(p["end_time"] / p["particle_dt"])} pre no post no every 1 precice_advance


precice_finalize
"""