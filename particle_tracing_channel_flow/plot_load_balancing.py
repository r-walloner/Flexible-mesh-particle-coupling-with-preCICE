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
plt.figure(figsize=(10, 4))
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

# Ideal load balancing for 4 ranks and 768 particles
plt.plot(
    [0, 0.22],
    [768 / 4, 768 / 4],
    linestyle="--",
    color="black",
    label="Ideal load balancing",
)

for rank, particles in particles_in_rank.items():
    times = particles.keys()
    particles_counts = particles.values()
    plt.plot(times, particles_counts, label=f"Rank {rank}")

plt.xlabel("Time [s]")
plt.ylabel("Particle count per rank")

plt.xlim(0, 0.22)
plt.ylim(0, 680)


plt.savefig(script_dir / "figures" / "channel_tracing_load_balancing.pdf", bbox_inches="tight")
plt.show()
