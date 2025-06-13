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
number_of_bins = 100 # number of measurement bins along the x-axis
max_time = 2.617 # duration of the simulation
particle_radius = 1.5e-3 # particle properties used for flux calculation
particle_density = 2505
particle_mass = 4/3 * pi * particle_radius**3 * particle_density

flux = np.zeros(number_of_bins)
previous_points: np.ndarray = None

# Iterate over timesteps
timestep_files = list(particle_path.glob("mdb_*.vtu"))
for timestep_file in tqdm(timestep_files, desc="Processing timesteps"):
    timestep: pv.UnstructuredGrid = pv.read(timestep_file)

    # Keep track of each particle's position in the previous timestep
    if previous_points is None:
        previous_points = timestep.points

    for point, previous_point in zip(timestep.points, previous_points):
        # Determine which bin the particle is in based on its x-coordinate
        bin_index = int((point[0] - x_min) / (x_max - x_min) * number_of_bins)
        if bin_index < 0 or bin_index >= number_of_bins:
            print(f"Warning: Particle at {point[0]} is out of bounds for bins [{x_min}, {x_max}]")
            continue
        

        # Check if the particle crossed the plane of measurement
        if (point[1] > plane_of_measurement and previous_point[1] <= plane_of_measurement):
            # Crossing in positive y-direction
            flux[bin_index] += particle_mass
        elif (point[1] < plane_of_measurement and previous_point[1] >= plane_of_measurement):
            # Crossing in negative y-direction
            flux[bin_index] -= particle_mass

    previous_points = timestep.points

# Normalize the flux by the total time and the area of the plane
flux /= max_time * (y_max - y_min) * (x_max - x_min)

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
        'average_flux': flux.tolist()
    }, f, indent=4)
