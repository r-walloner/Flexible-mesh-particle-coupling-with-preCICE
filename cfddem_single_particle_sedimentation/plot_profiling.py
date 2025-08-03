import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pathlib


script_dir = pathlib.Path(__file__).parent.resolve()
data_dir = script_dir / "data" / "profiling"

# # Plot
# ax = pivot.plot(kind="bar", stacked=True, figsize=(12, 6))
# plt.ylabel("Total Duration")
# plt.title("Total Runtime per Event (Stacked by Solver)")
# plt.legend(title="Event", bbox_to_anchor=(1.05, 1), loc="upper left")
# plt.tight_layout()

# Set up plot
plt.figure(figsize=(10, 6), dpi=250)
plt.title("Runtime breakdown")
plt.grid(True, axis="y")
plt.ylabel("Runtime [s]")


# Load and plot profiling data
files = list(data_dir.glob("*.csv"))
files.sort()
for file in files:
    # Read CSV
    df = pd.read_csv(file)

    # Group by participant and event, sum durations
    grouped = df.groupby(["participant", "event"])["duration"].sum().reset_index()
    pivot = grouped.pivot(index="participant", columns="event", values="duration").fillna(0)

    bottom = 0
    # Fluid
    for event in ["advance", "solver.advance"]:
        plt.bar(pivot.index, pivot[event], bottom=bottom, label=event)
        bottom += pivot[event]

# Output plot
plt.legend()
plt.show()