import numpy as np
import pathlib
import matplotlib.pyplot as plt
import json

script_dir = pathlib.Path(__file__).parent.resolve()
data_dir = script_dir / "data" / "particle_flux"
reference_data_path = script_dir / "reference_data" / "Link 2005.csv"

# Set up plot
plt.figure(figsize=(10, 6), dpi=250)
plt.title("Average Particle Z-Flux at z=0.13")
plt.grid(True)

xlim_max = None

# Load and plot flux data
for file in data_dir.glob("*.json"):
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

# Load and plot reference data
reference_data = np.genfromtxt(reference_data_path, delimiter=",", skip_header=2)
link_exp = reference_data[:, 0:2][~np.isnan(reference_data[:, 0])]
link_sim_min = reference_data[:, 2:4][~np.isnan(reference_data[:, 2])]
link_sim_old = reference_data[:, 4:6][~np.isnan(reference_data[:, 4])]
link_sim_koch = reference_data[:, 6:8][~np.isnan(reference_data[:, 6])]

link_exp = link_exp[link_exp[:, 0].argsort()]
link_sim_min = link_sim_min[link_sim_min[:, 0].argsort()]
link_sim_old = link_sim_old[link_sim_old[:, 0].argsort()]
link_sim_koch = link_sim_koch[link_sim_koch[:, 0].argsort()]

plt.scatter(link_exp[:,0], link_exp[:,1], label='Link 2005 - exp', color="black")
# plt.plot(link_sim_min[:,0], link_sim_min[:,1], label='Link 2005 - sim - min', color="black", linestyle="-.")
plt.plot(link_sim_old[:,0], link_sim_old[:,1], label='Link 2005 - sim - gidaspow', color="black", linestyle=":")
plt.plot(link_sim_koch[:,0], link_sim_koch[:,1], label='Link 2005 - sim - koch', color="black", linestyle="--")

# Output plot
plt.xlim(0, xlim_max)
plt.xticks(np.linspace(0, xlim_max, 7))
plt.yticks(np.linspace(-800, 800, 9))
plt.legend()
plt.show()