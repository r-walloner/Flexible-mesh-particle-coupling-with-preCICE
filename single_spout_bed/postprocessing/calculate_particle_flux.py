from math import pi
import pathlib
import pyvista as pv
import numpy as np
from tqdm import tqdm
import json

particle_path = pathlib.Path(__file__).parent.parent / "particle-liggghts" / "post"
x_min = -.075 # minimum and maximum x-coordinates of the domain
x_max = +.075
y_min = -.0075 # minimum and maximum y-coordinates of the domain
y_max = +.0075
plane_of_measurement = 0.13  # z-coordinate of the plane where flux is measured
number_of_bins = 32 # number of measurement bins along the x-axis
max_time = 2.617 # duration of the simulation
particle_radius = 1.5e-3 # particle properties used for flux calculation
particle_density = 2505
particle_mass = 4/3 * pi * particle_radius**3 * particle_density

positive_flux = np.zeros(number_of_bins)
negative_flux = np.zeros(number_of_bins)

previous_positions: dict[int, np.array] = None

# Iterate over timesteps
timestep_files = list(particle_path.glob("mdb_*.vtu"))
for timestep_file in tqdm(timestep_files, desc="Processing timesteps"):
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
positive_flux /= max_time * (y_max - y_min) * (x_max - x_min)
negative_flux /= max_time * (y_max - y_min) * (x_max - x_min)
net_flux = positive_flux - negative_flux

# Write results to file
with open(particle_path / 'particle_flux.json', 'w') as f:
    json.dump({
        'x_min': x_min,
        'x_max': x_max,
        'y_min': y_min,
        'y_max': y_max,
        'plane_of_measurement': plane_of_measurement,
        'number_of_bins': number_of_bins,
        'max_time': max_time,
        'number_of_timesteps': len(timestep_files),
        'particle_radius': particle_radius,
        'particle_density': particle_density,
        'particle_mass': particle_mass,
        'positive_flux': positive_flux.tolist(),
        'negative_flux': negative_flux.tolist(),
        'net_flux': net_flux.tolist(),
    }, f, indent=4)
