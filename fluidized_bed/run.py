import pathlib
import subprocess
from sys import argv
import json
from tqdm import tqdm
from math import ceil


script_dir = pathlib.Path(__file__).parent.resolve()
runs_dir = script_dir / "runs"


def run(run_dir: pathlib.Path):
    openfoam_log = run_dir / "fluid-openfoam" / "fluid-openfoam.log"
    liggghts_log = run_dir / "particle-liggghts" / "log.liggghts"

    if openfoam_log.exists() or liggghts_log.exists():
        # If there are logs, it means the run was already executed, so we skip it
        print(f"{run_dir.name} was allready run, skipping")
        return
    
    # Load parameters
    with open(run_dir / "parameters.json", "r") as f:
        parameters = json.load(f)

    progress = tqdm(desc=run_dir.name, total=int(ceil(parameters["end_time"] / parameters["fluid_dt"])))

    # Start the solvers
    fluid_command = ["openfoam2412", "./run.sh"]
    n_subdomains = parameters["fluid_subdomains"][0] * parameters["fluid_subdomains"][1] * parameters["fluid_subdomains"][2]
    if n_subdomains > 1:
        fluid_command.append("-parallel")
    fluid_process = subprocess.Popen(
        fluid_command,
        cwd=run_dir / "fluid-openfoam",
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
    )

    particle_process = subprocess.Popen(
        ["bash", "./run.sh"],
        cwd=run_dir / "particle-liggghts",
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        text=True,
    )

    # Continuously listed to the output of both processes
    while True:
        # Print the precice status messages
        fluid_output = fluid_process.stdout.readline()
        if "[impl::ParticipantImpl]:460" in fluid_output:
            progress.update()

        # Check if both processes have finished
        fluid_returncode = fluid_process.poll()
        particle_returncode = particle_process.poll()
        if fluid_returncode is not None and particle_returncode is not None:
            progress.close()
            if fluid_returncode != 0:
                print("WARNING: Fluid solver exited with non-zero return code")
            if particle_returncode != 0:
                print("WARNING: Particle solver exited with non-zero return code")
            break

def run_all():
    # Find all runs
    runs = list(run for run in runs_dir.iterdir() if run.is_dir())
    runs.sort(key=lambda x: x.name)

    # Run all runs
    for run_dir in runs:
        run(run_dir)


if __name__ == "__main__":
    if len(argv) == 1 or argv[1] == "all":
        run_all()
    else:
        # Get paths from command line arguments
        runs = [pathlib.Path(arg) for arg in argv[1:]]
        for run_dir in runs:
            if not run_dir.is_absolute():
                run_dir = script_dir / run_dir
            run(run_dir)
