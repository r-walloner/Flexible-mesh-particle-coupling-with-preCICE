import pathlib
import json
from typing import TypedDict
import shutil
import importlib.util

script_dir = pathlib.Path(__file__).parent.resolve()
runs_dir = script_dir / "runs"
template_dir = script_dir / "template"


class Parameters(TypedDict):
    solver: str
    end_time: float
    fluid_dt: float
    fluid_cells: tuple[int, int, int]
    fluid_subdomains: int
    fluid_viscosity: float
    fluid_density: float
    fluid_background_velocity: float
    fluid_spout_velocity: float
    particle_dt: float
    particle_subdomains: str
    particle_drag_model: str
    particle_diameter: float
    particle_density: float
    particle_contact_model: str
    particle_youngs_modulus: float
    particle_poissons_ratio: float
    particle_restitution: float
    particle_friction: float
    particl_characteristic_velocity: float
    particle_count: int
    particle_insert_velocity: tuple[float, float, float]
    particle_settling_gravity: float
    particle_settling_time: float
    read_mapping: str
    read_mapping_radius: float
    write_mapping: str
    write_mapping_radius: float
    output_interval: float


# Instantiate the case with the given parameters
def generate_run(p: Parameters):
    run_name = "_".join(
        [
            f"mesh-{p['fluid_cells'][0]}x{p['fluid_cells'][1]}x{p['fluid_cells'][2]}",
            f"particle-{p['particle_diameter']}",
            p["solver"],
            f"read-{p['read_mapping']}{'-' + str(p['read_mapping_radius']) if p['read_mapping_radius'] else ''}",
            f"write-{p['write_mapping']}{'-' + str(p['write_mapping_radius']) if p['write_mapping_radius'] else ''}",
        ]
    )

    # Create the run directory
    run_dir = runs_dir / run_name
    try:
        run_dir.mkdir(exist_ok=False, parents=True)
    except FileExistsError:
        print(f"{run_name} already exists, skipping")

    print(f"Generating {run_name}")

    # Write parameters to file
    with open(run_dir / "parameters.json", "w") as f:
        json.dump(p, f, indent=4)

    # Copy the template case
    shutil.copytree(template_dir, run_dir, dirs_exist_ok=True)

    # Instantiate the template files with parameters
    for template_file in run_dir.rglob("*.py"):
        # Load the template file as a module
        spec = importlib.util.spec_from_file_location(template_file.name, template_file)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        # Execute the generate() funcion from the template and write the output
        output = module.generate(p)
        if output is not None:
            output_file = template_file.with_suffix("")
            with open(output_file, "w") as f:
                f.write(output)

            # Make output file executable, if its a script
            if output_file.suffix in [".sh", ".bash"]:
                output_file.chmod(output_file.stat().st_mode | 0o111)

        # Remove the template file
        template_file.unlink()

    


# Set default parameters
p = Parameters(
    solver="AndersonJacksonFoam",
    end_time=2,
    fluid_dt=1e-4,
    fluid_cells=(30, 250, 1),
    fluid_subdomains=1,
    fluid_viscosity=1.8e-5,
    fluid_density=1,
    fluid_background_velocity=1.5,
    fluid_spout_velocity=30,
    particle_dt=5e-5,
    particle_subdomains="8 1 1",
    particle_drag_model="gidaspow",
    particle_diameter=3e-3,
    particle_density=2505,
    particle_contact_model="model hooke tangential history",
    particle_youngs_modulus=1e7, # Is this right? Could not find it in Kloss or Link
    particle_poissons_ratio=0.23, # Is this right? Could not find it in Kloss or Link
    particle_restitution=0.97,
    particle_friction=0.1,
    particle_characteristic_velocity=10, # Is this right? Could not find it in Kloss or Link
    particle_count=24500,
    particle_insert_velocity=(0, -0.5, 0),
    particle_settling_gravity=100, # Amplified gravity to speed up settling
    particle_settling_time=0.3,
    read_mapping="nearest-neighbor",
    read_mapping_radius=None,
    write_mapping="coarse-graining",
    write_mapping_radius=12e-3,
    output_interval=1e-3,
)

generate_run(p)

# Generate runs with varying parameters

# Vary solver
# for solver in ["AndersonJacksonFoam", "pimpleFoam"]:
#     p["solver"] = solver

#     if solver == "AndersonJacksonFoam":
#         p["write_mapping"] = "coarse-graining"
#         p["write_mapping_radius"] = 4 * p["particle_diameter"]
#     elif solver == "pimpleFoam":
#         p["write_mapping"] = "nearest-neighbor"
#         p["write_mapping_radius"] = None

#     # Vary read mapping
#     for read_mapping in ["nearest-neighbor", "rbf"]:
#         p["read_mapping"] = read_mapping

#         if read_mapping == "nearest-neighbor":
#             read_radii = [None]
#         else:
#             read_radii = [n * p["particle_diameter"] for n in [1, 3, 6, 12]]
#         for read_radius in read_radii:
#             p["read_mapping_radius"] = read_radius

#             generate_run(p)
