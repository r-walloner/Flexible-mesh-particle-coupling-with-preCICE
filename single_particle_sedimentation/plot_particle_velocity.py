import numpy as np
import pathlib
import matplotlib.pyplot as plt
import json

script_dir = pathlib.Path(__file__).parent.resolve()
data_dir = script_dir / "data" / "particle_velocity"
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
plt.plot(song_sim[:,0], song_sim[:,1], label='Song Park', color="black")


# Load and plot flux data
files = list(data_dir.glob("*.json"))
files.sort()
for file in files:
    with open(file, "r") as f:
        data = json.load(f)

    time = data["time"]
    particle_velocity = data["particle_velocity"]

    plt.plot(time, particle_velocity, label=file.stem)



plt.legend()

plt.savefig(script_dir / "figures" / "sps_particle_velocity.pdf", bbox_inches="tight")
plt.show()
