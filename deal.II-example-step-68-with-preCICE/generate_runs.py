import pathlib

script_path = pathlib.Path(__file__).resolve()

def generate_run(path: pathlib.Path,
                 refinement=4, 
                 mapping="rbf-pum-direct", 
                 basis_function="compact_polynomial-c6", 
                 support_radius=0.5,
                 constraint="consisent"):
    """Generate a new run directory with the given parameters.
    
    This function creates the directory specified by `path` and writes the
    configuration files `precice-config.xml` and `parameters.prm` to it. The
    configuration files are generated with the given parameters."""

    if path.exists():
        print(f"Not generating run {path.relative_to(script_path.parent)}: Directory already exists")
        return
    
    path.mkdir(parents=True)

    # Write precice-config.xml
    with open(path / "precice-config.xml", "w") as file:
        file.write(f"""
<precice-configuration experimental="true">

  <data:vector name="Velocity" />

  <mesh name="Fluid-Mesh" dimensions="2">
    <use-data name="Velocity" />
  </mesh>

  <participant name="Fluid">
    <provide-mesh name="Fluid-Mesh" />
    <write-data name="Velocity" mesh="Fluid-Mesh" />
  </participant>

  <participant name="Particle">
    <receive-mesh name="Fluid-Mesh" from="Fluid" api-access="true" />
    <read-data name="Velocity" mesh="Fluid-Mesh" />
    <mapping:{mapping}
      direction="read"
      from="Fluid-Mesh"
      constraint="{constraint}">
      {f'<basis-function:{basis_function} support-radius="{support_radius}" />' if mapping == "rbf-pum-direct" else ""}
    </mapping:{mapping}>
  </participant>

  <m2n:sockets
    acceptor="Fluid"
    connector="Particle"
    exchange-directory="{path.resolve()}" />

  <coupling-scheme:serial-explicit>
    <participants first="Fluid" second="Particle" />
    <max-time value="4.0" />
    <time-window-size value="0.002" />
    <exchange data="Velocity" mesh="Fluid-Mesh" from="Fluid" to="Particle" />
  </coupling-scheme:serial-explicit>

</precice-configuration>
""")
        
    # Write parameters.prm
    with open(path / "parameters.prm", "w") as file:
        file.write(f"""
# Listing of Parameters
# ---------------------
subsection Particle Tracking Problem
  # End time of the simulation
  set Final time                    = 4

  # Refinement level of the fluid domain
  set Fluid refinement              = {refinement}
  set Output basename               = {path.name}_
  set Output directory              = {path.parent.resolve()}/

  # Iteration interval between which output results are written
  set Output interval               = 10

  # Refinement level of the particle domain
  set Particle grid refinement      = {refinement}

  # Refinement of the volumetric mesh used to insert the particles
  set Particle insertion refinement = 3

  # Iteration interval at which the mesh is load balanced
  set Repartition interval          = 5
  set Time step                     = 0.002
  set Velocity degree               = 1
end
""")

