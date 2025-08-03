import pathlib
from sys import argv
import pyvista as pv
import numpy as np
import json
import subprocess
from tqdm import tqdm
from math import pi

script_dir = pathlib.Path(__file__).parent.resolve()
runs_dir = script_dir / "runs"
out_dir = script_dir / "data"


# Merge and extract profiling events
def extract_profiling_events(run_dir: pathlib.Path):
    output_file = out_dir / "profiling" / f"{run_dir.name}.csv"
    if output_file.exists():
        print(f"Profiling data {output_file.name} already exists, skipping")
        return
    output_file.parent.mkdir(parents=True, exist_ok=True)

    print(f"Extracting profiling events for {run_dir.name}")

    subprocess.call(
        ["precice-profiling", "merge", "fluid-openfoam", "particle-liggghts"],
        cwd=run_dir
    )
    subprocess.call(
        ["precice-profiling", "export", "-o", output_file.as_posix()],
        cwd=run_dir
    )

    print()


if __name__ == "__main__":
    if len(argv) == 1 or argv[1] == "all":
        # postprocess all runs
        runs = list(run for run in runs_dir.iterdir() if run.is_dir())
        runs.sort(key=lambda x: x.name)
    else:
        runs = [pathlib.Path(arg) for arg in argv[1:]]

    for run_dir in runs:
        # extract_profiling_events(run_dir)
