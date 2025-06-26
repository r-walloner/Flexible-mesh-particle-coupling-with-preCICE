import pathlib
import pyvista as pv
import numpy as np
import json
import subprocess
import tqdm
from math import pi

script_dir = pathlib.Path(__file__).parent.resolve()
runs_dir = script_dir / "runs"
out_dir = script_dir / "data"


# domain boundaries and measurement parameters
# TODO get these from parameters.json (need to update generate.py first)
x_min = -.075 # m
x_max = +.075 # m
y_min = -.0075 # m
y_max = +.0075 # m
plane_of_measurement = 0.13 # z-coordinate [m] of the plane where flux is measured
number_of_bins = 32 # number of measurement bins along the x-axis


# Find runs
runs = list(run for run in runs_dir.iterdir() if run.is_dir())
runs.sort(key=lambda x: x.name)


# Extract the averaged particle flux across the measurement plane
def extract_particle_flux(run_dir: pathlib.Path):
    # Read run parameters
    with open(run_dir / "parameters.json", "r") as f:
        parameters = json.load(f)

    # Check if output file already exists
    output_file = out_dir / "particle_flux" / f"{run_dir.name}.json"
    if output_file.exists():
        print(f"Particle flux {output_file.name} allready exists, skipping")
        return
    
    print(f"Extracting particle flux for {run_dir.name}")

    # Check if particle data exists
    data_dir = run_dir / "particle-liggghts" / "out"
    if not data_dir.exists():
        print(f"{run_dir.name} does not contain data, skipping")
        return

    positive_flux = np.zeros(number_of_bins)
    negative_flux = np.zeros(number_of_bins)

    previous_positions: dict[int, np.array] = None

    particle_volume = 4/3 * pi * (parameters["particle_diameter"] / 2)**3
    particle_mass = particle_volume * parameters["particle_density"]

    # Iterate over timesteps
    for timestep_file in tqdm(data_dir.glob("particles_*.vtu")):
        timestep: pv.UnstructuredGrid = pv.read(timestep_file)

        # Store particle positions by their IDs
        positions_by_id: dict[int, np.array] = {}
        for position, particle_id in zip(timestep.points, timestep.point_data["id"]):
            positions_by_id[particle_id] = position

        # Keep track of particle positions in the previous timestep
        if previous_positions is None:
            previous_positions = positions_by_id

        # Iterate over all particles
        for id, position in positions_by_id.items():
            previous_position = previous_positions[id]

            # Determine which bin the particle is in based on its x-coordinate
            bin_index = int((position[0] - x_min) / (x_max - x_min) * number_of_bins)
            if bin_index < 0 or bin_index >= number_of_bins:
                print(f"Warning: Particle at {position[0]} is out of bounds for bins [{x_min}, {x_max}]")
                continue
            
            # Check if the particle crossed the plane of measurement
            if (position[1] > plane_of_measurement and previous_position[1] <= plane_of_measurement):
                # Crossing in positive y-direction
                positive_flux[bin_index] += particle_mass
            elif (position[1] < plane_of_measurement and previous_position[1] >= plane_of_measurement):
                # Crossing in negative y-direction
                negative_flux[bin_index] += particle_mass

        previous_positions = positions_by_id

    # Normalize the flux by the total time and the area of the plane
    area = (x_max - x_min) * (y_max - y_min) / number_of_bins # area of each bin
    positive_flux /= parameters["end_time"] * area
    negative_flux /= parameters["end_time"] * area
    net_flux = positive_flux - negative_flux

    # Write results to file
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump({
            'run_parameters': parameters,
            'x_min': x_min,
            'x_max': x_max,
            'y_min': y_min,
            'y_max': y_max,
            'plane_of_measurement': plane_of_measurement,
            'number_of_bins': number_of_bins,
            'positive_flux': positive_flux.tolist(),
            'negative_flux': negative_flux.tolist(),
            'net_flux': net_flux.tolist(),
        }, f, indent=4)

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


# Postprocess each run
for run_dir in runs:
    extract_particle_flux(run_dir)
    extract_profiling_events(run_dir)
