import pathlib
import matplotlib.pyplot as plt
from tqdm import tqdm

script_dir = pathlib.Path(__file__).parent.resolve()
log_path = script_dir / "runs" / "Ours_Koch_and_Hill" / "particle-liggghts" / "log.liggghts"


# Extract particle counts from log file
with open(log_path, "r") as f:
    lines = f.readlines()

    # Skip lines until the coupling is initialized
    while "fluid_coupling" not in lines[0]:
        lines.pop(0)

    times: list[float] = []
    nlocal_max: list[float] = []
    nlocal_min: list[float] = []

    for i, line in tqdm(enumerate(lines), total=len(lines)):
        if "Nlocal:" in line:
            # Parse relevant data
            timestep_line = lines[i - 10]
            fields = timestep_line.split(" ")
            fields = [f for f in fields if f]
            time = float(fields[2])
            times.append(time)
            
            fields = line.split(" ")
            fields = [f for f in fields if f]
            nlocal_max.append(float(fields[3]))
            nlocal_min.append(float(fields[5]))


# Subtract 0.8001 from all times to account for initial settling phase
times = [t - 0.8001 for t in times]


# Plot data
plt.figure(figsize=(10, 3.5))
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

plt.grid(axis="y")

# Ideal load balancing for 54 ranks and 24500 particles
plt.plot(
    [0, 5],
    [24500 / 54, 24500 / 54],
    linestyle="--",
    color="black",
    label="Ideal load balancing",
)

# Only keep every 100th element to reduce the number of points
times = times[::100]
nlocal_max = nlocal_max[::100]
nlocal_min = nlocal_min[::100]

plt.plot(times, nlocal_max, label="Maximum", color="tab:blue", alpha=1, linewidth=1)
plt.plot(times, nlocal_min, label="Minimum", color="tab:blue", alpha=1, linewidth=1)
plt.fill_between(times, nlocal_min, nlocal_max, color="tab:blue", alpha=0.2)

plt.xlabel("Time [s]")
plt.ylabel("Particle count per rank")

plt.xlim(0, 5)
plt.ylim(0, 600)



plt.savefig(script_dir.parents[1] / "figures" / "fb_load_balancing.pdf", bbox_inches="tight")
plt.show()
