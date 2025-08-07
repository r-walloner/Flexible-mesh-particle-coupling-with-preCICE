import pathlib
import matplotlib.pyplot as plt

script_dir = pathlib.Path(__file__).parent.resolve()
log_path = script_dir / "particle-deal.II" / "particle.log"


# Extract particle counts from log file
with open(log_path, "r") as f:
    lines = f.readlines()

    particles_in_rank: dict[int, dict[float, int]] = {}

    for line in lines:
        if "load balancing info" in line:
            # Parse relevant data
            time = float(line.split("\ttime: ")[1].split("\t")[0])
            rank = int(line.split("\trank: ")[1].split("\t")[0])
            local_particles = int(line.split("local particles: ")[1])

            if rank not in particles_in_rank:
                particles_in_rank[rank] = {}
            particles_in_rank[rank][time] = local_particles


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

# plt.grid()

# Add vertical lines for repartitioning events
dt = 1e-3
repartition_interval = 5
repartition_times = [(i * repartition_interval + 1) * dt for i in range(1, 50)]
for t in repartition_times:
    plt.axvline(x=t, color="gray", linestyle="dotted", linewidth=1)

# Ideal load balancing for 4 ranks and 12288 particles
plt.plot(
    [0, 0.22],
    [12288 / 4, 12288 / 4],
    linestyle="--",
    color="black",
    label="Ideal load balancing",
)

# Set line color pallete to custom colors

# Define a custom color palette
custom_colors = ["#4747db", "#00d4a9", "#ff9500", "#e04d4d"]

for idx, (rank, particles) in enumerate(sorted(particles_in_rank.items())):
    times = particles.keys()
    particles_counts = particles.values()
    color = custom_colors[idx % len(custom_colors)]
    plt.plot(times, particles_counts, label=f"Rank {rank}", color=color)

plt.xlabel("Time [s]")
plt.ylabel("Particle count per rank")

plt.xlim(1e-3, 0.22)
plt.ylim(0, 10800)



plt.savefig(script_dir / "figures" / "channel_tracing_load_balancing.pdf", bbox_inches="tight")
plt.show()
