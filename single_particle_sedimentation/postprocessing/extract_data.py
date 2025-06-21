from math import pi
import pathlib
import pyvista as pv
import numpy as np
from tqdm import tqdm
import json

# Input and output paths
script_dir = pathlib.Path(__file__).parent
particle_path = script_dir.parent / "particle-liggghts" / "out"
timestep_files = list(particle_path.glob("particles_*.vtu"))
output_path = script_dir / "data" / "particles.json"

timestep_size = 3e-3 # time [s] between two timestep files (caution: this is not the solver timestep)
max_time = timestep_size * len(timestep_files) # total time of the simulation

time = np.arange(0, max_time, timestep_size)
particle_velocity = []

# Iterate over timesteps
for timestep_file in tqdm(timestep_files, desc="Processing timesteps"):
    timestep: pv.UnstructuredGrid = pv.read(timestep_file)

    if timestep.number_of_points < 1:
        continue # We have lost our particle :(

    velocity = np.linalg.norm(timestep.point_data["v"][0])
    particle_velocity.append(velocity)

# Write results to file
with open(output_path, 'w') as f:
    json.dump({
        "time": time.tolist(),
        "particle_velocity": particle_velocity,
    }, f, indent=4)
