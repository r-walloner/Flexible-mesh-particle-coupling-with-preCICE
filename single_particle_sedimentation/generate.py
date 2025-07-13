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
    coupling: str
    end_time: float
    fluid_dt: float
    fluid_cells: tuple[int, int, int]
    fluid_subdomains: int
    fluid_viscosity: float
    fluid_density: float
    particle_dt: float
    particle_diameter: float
    particle_density: float
    particle_drag_model: str
    read_mapping: str
    read_mapping_radius: float
    write_mapping: str
    write_mapping_radius: float
    output_interval: float
    output_compression:bool
    coupling_scheme: str
    precice_debug_log: bool
    precice_profiling: str



# Instantiate the case with the given parameters
def generate_run(p: Parameters, run_name: str = None):
    if run_name is None:
        mapping_abbreviations = {
            "nearest-neighbor": "NN",
            "rbf": "RBF",
            "coarse-graining": "CG",
        }
        run_name = "_".join(
            [
                f"mesh-{p['fluid_cells'][0]}x{p['fluid_cells'][1]}x{p['fluid_cells'][2]}",
                p["solver"],
                p["coupling"].replace("_", "-"),
                p["particle_drag_model"].replace("_", "-"),
                f"read-{mapping_abbreviations[p["read_mapping"]]}{"-" + str(p["read_mapping_radius"]) if p["read_mapping_radius"] else ""}",
                f"write-{mapping_abbreviations[p["write_mapping"]]}{"-" + str(p["write_mapping_radius"]) if p["write_mapping_radius"] else ""}",
                f"d-{p['particle_diameter']}",
            ]
        )

    # Create the run directory
    run_dir = runs_dir / run_name
    try:
        run_dir.mkdir(exist_ok=False, parents=True)
    except FileExistsError:
        print(f"{run_name} already exists, skipping")
        return

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
    coupling="explicit",
    end_time=0.25,
    fluid_dt=1e-3,
    fluid_cells=(6, 18, 6),
    fluid_subdomains=1,
    fluid_viscosity=1.002e-3,
    fluid_density=998.25,
    particle_dt=5e-5,
    particle_diameter=2e-3,
    particle_density=2463,
    particle_drag_model="zhao_shan",
    read_mapping="rbf",
    read_mapping_radius=0.5,
    write_mapping="coarse-graining",
    write_mapping_radius=8e-3,
    output_interval=1e-3,
    output_compression=False,
    coupling_scheme="serial-explicit",
    precice_debug_log=False,
    precice_profiling="off",
)

# generate_run(p, "debug")

# Generate runs with varying parameters

for fluid_cells in [(6, 18, 6), (25, 75, 25)]:
    p["fluid_cells"] = fluid_cells

    for particle_drag_model in ["zhao_shan", "gidaspow"]:
        p["particle_drag_model"] = particle_drag_model

        for read_mapping in ["rbf"]:
            p["read_mapping"] = read_mapping
            p["read_mapping_radius"] = 0.5 if read_mapping == "rbf" else None

            generate_run(p)