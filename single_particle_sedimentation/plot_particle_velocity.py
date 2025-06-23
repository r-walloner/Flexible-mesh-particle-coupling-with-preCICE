import numpy as np
import pathlib
import matplotlib.pyplot as plt
import json

script_dir = pathlib.Path(__file__).parent.resolve()
data_dir = script_dir / "data" / "particle_velocity"
reference_file = script_dir / "data" / "particle_velocity_reference" / "Song Park 2020 - coarse.csv"


# Set up plot
plt.figure(figsize=(10, 6), dpi=250)
plt.title("Settling velocity of a single particle")
plt.grid(True)
plt.xlabel("Time [s]")
plt.ylabel("Velocity [m/s]")
plt.xlim(0, 0.25)
plt.ylim(0, 0.3)


# Load and plot reference data
reference_data = np.genfromtxt(reference_file, delimiter=",", skip_header=2)
song_sim = reference_data[:, 0:2][~np.isnan(reference_data[:, 0])]
song_theoretical = reference_data[:, -2:][~np.isnan(reference_data[:, -2])]

song_sim = song_sim[song_sim[:, 0].argsort()]
song_theoretical = song_theoretical[song_theoretical[:, 0].argsort()]

plt.plot(song_sim[:,0], song_sim[:,1], label='Song Park 2020', color="black")
plt.plot(song_theoretical[:,0], song_theoretical[:,1], label='Empirical Correlation', color="black", linestyle="--")


# Load and plot flux data
files = list(data_dir.glob("*.json"))
files.sort()
for file in files:
    with open(file, "r") as f:
        data = json.load(f)

    time = data["time"]
    particle_velocity = data["particle_velocity"]

    plt.plot(time, particle_velocity, label=file.stem)



# Output plot
plt.legend()
plt.show()
