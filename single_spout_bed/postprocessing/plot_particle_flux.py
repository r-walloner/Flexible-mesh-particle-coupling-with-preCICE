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
flux = np.array(data['average_flux'])

plt.figure(figsize=(10, 6), dpi=600)
plt.plot(np.linspace(x_min, x_max, number_of_bins), flux, label='Average Particle Flux')
plt.show()