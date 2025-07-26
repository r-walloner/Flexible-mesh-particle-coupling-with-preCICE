import pathlib
from sys import argv
import matplotlib.pyplot as plt
import numpy as np

script_dir = pathlib.Path(__file__).parent.resolve()
reference_data_path = script_dir / "reference_data" / "Link 2005.csv"

bins = 30
rho_p = 2505

#Set up plot
plt.figure(figsize=(10, 7))
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Liberation Serif"],
        "mathtext.fontset": "stix",
        "font.size": 16,
        "axes.labelsize": 16,
        "axes.titlesize": 16,
        "xtick.labelsize": 16,
        "ytick.labelsize": 16,
        "legend.fontsize": 16,
        "xtick.color": "black",
        "ytick.color": "black",
        "axes.labelcolor": "black",
    }
)
# plt.title("Settling velocity of a single particle")
plt.grid(True)
plt.xlabel("x [m]")
plt.ylabel(r"Average particle z-flux [$\text{kg}/(\text{m}^2 \text{s})$]")
# plt.xlim(0, 0.25)
# plt.ylim(0, 0.28)


# Determine what runs to read
if len(argv) > 1:
    # If a runs are specified on the command line, use the data from those
    runs = [pathlib.Path(r) for r in argv[1:]]
else:
    # Otherwise, use all runs in the runs directory
    runs = list(script_dir.glob("runs/*"))


for run in runs:
    velocity_file = run / "particle-liggghts" / "out" / "average_velocity"
    if not velocity_file.is_file():
        print(f"Skipping {run.name} as it does not contain a velocity file")
        continue

    # Read in average velocity data
    velocity: dict[int, list[float]] = {}
    count: dict[int, list[float]] = {}
    x_coords = []
    with open(velocity_file, "r") as f:
        lines = f.readlines()

        # Skip header
        lines = lines[3:]

        while lines:
            line = lines.pop(0)
            values = line.split(" ")
            timestep = int(values[0])
            velocity[timestep] = []
            count[timestep] = []

            for bin in range(bins):
                line = lines.pop(0)
                values = line.strip().split(" ")
                velocity[timestep].append(float(values[3]))
                count[timestep].append(float(values[2]))

                # Only read x-coordinates from the first timestep
                if len(x_coords) < bins:
                    x_coords.append(float(values[1]))


    # Calculate flux on each timestep
    flux: dict[int, list[float]] = {}
    for timestep in velocity.keys():
        flux[timestep] = []
        total_count = 0
        for bin in range(bins):
            flux[timestep].append(velocity[timestep][bin] * count[timestep][bin] * rho_p)
            total_count += count[timestep][bin]
        
        for bin in range(bins):
            flux[timestep][bin] /= (total_count / bins) if total_count > 0 else 1


    # Average over timestep range and weight with count
    flux_avg = [0] * bins

    averaging_range = velocity.keys()  # Use all timesteps
    averaging_range = range(16100, 116101, 100) # Select specific range

    for bin in range(bins):
        for timestep in averaging_range:
            flux_avg[bin] += flux[timestep][bin]
        flux_avg[bin] /= len(averaging_range)

    # Plot flux
    plt.plot(x_coords, flux_avg, label=run.name)


# Load and plot reference data
reference_data = np.genfromtxt(reference_data_path, delimiter=",", skip_header=2)
link_exp = reference_data[:, 0:2][~np.isnan(reference_data[:, 0])]
link_sim_min = reference_data[:, 2:4][~np.isnan(reference_data[:, 2])]
link_sim_old = reference_data[:, 4:6][~np.isnan(reference_data[:, 4])]
link_sim_koch = reference_data[:, 6:8][~np.isnan(reference_data[:, 6])]

# Offset x-coordinates of reference data to match simulation
link_exp[:, 0] -= 0.075
link_sim_min[:, 0] -= 0.075
link_sim_old[:, 0] -= 0.075
link_sim_koch[:, 0] -= 0.075

link_exp = link_exp[link_exp[:, 0].argsort()]
link_sim_min = link_sim_min[link_sim_min[:, 0].argsort()]
link_sim_old = link_sim_old[link_sim_old[:, 0].argsort()]
link_sim_koch = link_sim_koch[link_sim_koch[:, 0].argsort()]

plt.scatter(link_exp[:,0], link_exp[:,1], label='Link 2005 - exp', color="black")
# plt.plot(link_sim_min[:,0], link_sim_min[:,1], label='Link 2005 - sim - min', color="black", linestyle="-.")
plt.plot(link_sim_old[:,0], link_sim_old[:,1], label='Link 2005 - sim - gidaspow', color="black", linestyle=":")
plt.plot(link_sim_koch[:,0], link_sim_koch[:,1], label='Link 2005 - sim - koch', color="black", linestyle="--")

plt.legend()
plt.show()
    

                

    

