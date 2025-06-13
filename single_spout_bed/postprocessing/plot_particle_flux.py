import numpy as np
import pathlib
import matplotlib.pyplot as plt
import json

particle_path = pathlib.Path(__file__).parent.parent / "particle-liggghts" / "post"

with open(particle_path / "particle_flux.json", "r") as f:
    data = json.load(f)

x_min = data['x_min']
x_max = data['x_max']
number_of_bins = data['number_of_bins']
positive_flux = np.array(data['positive_flux'])
negative_flux = np.array(data['negative_flux'])
net_flux = np.array(data['net_flux'])

plt.figure(figsize=(10, 6), dpi=600)
plt.title("Average Particle Z-Flux at z=0.13")
plt.grid(True)
plt.plot(np.linspace(x_min, x_max, number_of_bins), positive_flux, label='Flux in positive z-direction', color='red')
plt.plot(np.linspace(x_min, x_max, number_of_bins), negative_flux, label='Flux in negative z-direction', color='blue')
plt.plot(np.linspace(x_min, x_max, number_of_bins), net_flux, label='Net flux', linestyle='--', color='black')
plt.legend()
plt.show()