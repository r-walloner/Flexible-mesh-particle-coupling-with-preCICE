import pathlib
import matplotlib.pyplot as plt

script_dir = pathlib.Path(__file__).parent.resolve()
log_path = script_dir / "particle-deal.II" / "particle.log"


# Extract particle counts from log file
with open(log_path, "r") as f:
    lines = f.readlines()

    times: list[int] = [0]
    particles_in_rank: dict[int, list[int]] = {}

    for line in lines:
        if "repartitioning results" in line:
            # Parse relevant data
            time = float(line.split("\ttime: ")[1].split("\t")[0])
            rank = int(line.split("\trank: ")[1].split("\t")[0])
            local_particles = int(line.split("local particles: ")[1])

            if rank not in particles_in_rank:
                particles_in_rank[rank] = []
            particles_in_rank[rank].append(local_particles)
            
            if time not in times:
                times.append(time)

# Plot data
plt.figure(figsize=(10, 6))
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

for rank, particles in particles_in_rank.items():
    plt.plot(times, particles, label=f"Rank {rank}")

plt.xlabel("Time [s]")
plt.ylabel("Particle count")

plt.savefig(script_dir / "channel_tracing_load_balancing.pdf", bbox_inches="tight")
plt.show()
