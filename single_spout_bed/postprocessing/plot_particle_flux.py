import numpy as np
import pathlib
import matplotlib.pyplot as plt
import json

script_dir = pathlib.Path(__file__).parent
particle_path = script_dir.parent / "figures" / "flux_data"

# Set up plot
plt.figure(figsize=(10, 6), dpi=250)
plt.title("Average Particle Z-Flux at z=0.13")
plt.grid(True)

xlim_max = None

# Load and plot flux data
for file in particle_path.glob("*.json"):
    with open(file, "r") as f:
        data = json.load(f)

    x_min = data['x_min']
    x_max = data['x_max']
    number_of_bins = data['number_of_bins']
    bins = np.linspace(0, x_max - x_min, number_of_bins)

    positive_flux = np.array(data['positive_flux'])
    negative_flux = np.array(data['negative_flux'])
    net_flux = np.array(data['net_flux'])

    if xlim_max is None:
        xlim_max = x_max - x_min
    else:
        xlim_max = max(x_max - x_min, xlim_max)

    # plt.plot(bins, positive_flux, label='Flux in positive z-direction', color='red')
    # plt.plot(bins, negative_flux, label='Flux in negative z-direction', color='blue')
    plt.plot(bins, net_flux, label=file.stem)

# Load reference data
reference_data_path = script_dir / "reference_data" / "Link_2005_experiment.csv"
reference_data = np.loadtxt(reference_data_path, delimiter=",")
plt.scatter(reference_data[:,1], reference_data[:,0], label='Link 2005', color='black')

# Output plot
plt.xlim(0, xlim_max)
plt.xticks(np.linspace(0, xlim_max, 7))
plt.yticks(np.linspace(-800, 800, 9))
plt.legend()
plt.show()