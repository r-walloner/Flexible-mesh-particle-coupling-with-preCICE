import numpy as np
import pathlib
import matplotlib.pyplot as plt
import json
import pyvista as pv

script_dir = pathlib.Path(__file__).parent.resolve()
runs_dir = script_dir / "runs"
reference_file = script_dir / "reference_data" / "Song Park 2020 - fine.csv"


# Set up plot
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
plt.xlabel("Time [s]")
plt.ylabel("Particle velocity [m/s]")
plt.xlim(0, 0.25)
plt.ylim(0, 0.28)


# Load and plot reference data
reference_data = np.genfromtxt(reference_file, delimiter=",", skip_header=2)
song_sim = reference_data[:, 0:2][~np.isnan(reference_data[:, 0])]
song_theoretical = reference_data[:, -2:][~np.isnan(reference_data[:, -2])]

song_sim = song_sim[song_sim[:, 0].argsort()]
song_theoretical = song_theoretical[song_theoretical[:, 0].argsort()]

plt.plot(song_theoretical[:,0], song_theoretical[:,1], label='Empirical Correlation', color="black", linestyle="--")
plt.plot(song_sim[:,0], song_sim[:,1], label='Song and Park', color="black")


# Load and plot particle velocity data
runs = sorted(run for run in runs_dir.iterdir() if run.is_dir())
for run in runs:

    # Read run parameters
    with open(run / "parameters.json", "r") as f:
        parameters = json.load(f)

    time_list = []
    velocity_list = []

    # Check if particle data exists
    data_dir = run / "particle-liggghts" / "out"
    if not data_dir.exists():
        print(f"{run.name} does not contain data, skipping")
        continue

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

    # Plot
    plt.plot(time_list, velocity_list, label=run.name)

plt.legend()

plt.savefig(script_dir / "figures" / "sps_particle_velocity.pdf", bbox_inches="tight")
plt.show()
