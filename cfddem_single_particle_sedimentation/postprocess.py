import pathlib
from sys import argv
import pyvista as pv
import numpy as np
import json
import subprocess

script_dir = pathlib.Path(__file__).parent.resolve()
runs_dir = script_dir / "runs"
out_dir = script_dir / "data"


# Extract the particle velocity over time from the simulation data
def extract_particle_velocity(run_dir: pathlib.Path):
    # Read run parameters
    with open(run_dir / "parameters.json", "r") as f:
        parameters = json.load(f)

    time_list = []
    velocity_list = []

    # Check if output file already exists
    output_file = out_dir / "particle_velocity" / f"{run_dir.name}.json"
    if output_file.exists():
        print(f"Particle velocity {output_file.name} allready exists, skipping")
        return
    
    print(f"Extracting particle velocity for {run_dir.name}")

    # Check if particle data exists
    data_dir = run_dir / "particle-liggghts" / "out"
    if not data_dir.exists():
        print(f"{run_dir.name} does not contain data, skipping")
        return


    # Iterate over timesteps
    for timestep_file in data_dir.glob("particles_*.vtu"):
        timestep: pv.UnstructuredGrid = pv.read(timestep_file)

        if timestep.number_of_points < 1:
            continue  # We have lost our particle :(

        velocity = np.linalg.norm(timestep.point_data["v"][0])
        velocity_list.append(velocity)

        timestep = int(timestep_file.stem.split("_")[-1])
        time = timestep * parameters["particle_dt"]
        time_list.append(time)

    # Write results to file
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(
            {
                "time": time_list,
                "particle_velocity": velocity_list,
            },
            f,
            indent=4,
        )

    print()


# Merge and extract profiling events
def extract_profiling_events(run_dir: pathlib.Path):
    output_file = out_dir / "profiling" / f"{run_dir.name}.csv"
    if output_file.exists():
        print(f"Profiling data {output_file.name} already exists, skipping")
        return
    output_file.parent.mkdir(parents=True, exist_ok=True)

    print(f"Extracting profiling events for {run_dir.name}")

    subprocess.call(
        ["precice-profiling", "merge", "fluid-openfoam", "particle-liggghts"],
        cwd=run_dir
    )
    subprocess.call(
        ["precice-profiling", "export", "-o", output_file.as_posix()],
        cwd=run_dir
    )

    print()


if __name__ == "__main__":
    if len(argv) == 1 or argv[1] == "all":
        # postprocess all runs
        runs = list(run for run in runs_dir.iterdir() if run.is_dir())
        runs.sort(key=lambda x: x.name)
    else:
        runs = [pathlib.Path(arg) for arg in argv[1:]]

    for run_dir in runs:
        extract_particle_velocity(run_dir)
        # extract_profiling_events(run_dir)
