import pathlib
import subprocess
import concurrent.futures

script_path = pathlib.Path(__file__).resolve()

executable_path = script_path.parent / "build" / "step-68-with-precice"
if not executable_path.is_file():
    raise FileNotFoundError(f"Executable {executable_path} not found")


def run(path: pathlib.Path):
    """Execute one run with the config files in the given directory."""
    try:
        print(f"Running {path.relative_to(script_path.parent)}")

        # Check if config files exist
        solver_config_path = path / "parameters.prm"
        if not solver_config_path.is_file():
            raise FileNotFoundError(
                f"Solver config {solver_config_path.relative_to(script_path.parent)} not found"
            )

        precice_config_path = path / "precice-config.xml"
        if not precice_config_path.is_file():
            raise FileNotFoundError(
                f"Solver config {precice_config_path.relative_to(script_path.parent)} not found"
            )

        # Create log files
        fluid_log_path = path / "fluid.log"
        if fluid_log_path.exists():
            print(
                f"Log file {fluid_log_path.relative_to(script_path.parent)} already exists. Skipping run."
            )
            return
        fluid_log_path.touch()

        particle_log_path = path / "particle.log"
        if particle_log_path.exists():
            print(
                f"Log file {particle_log_path.relative_to(script_path.parent)} already exists. Skipping run."
            )
            return
        particle_log_path.touch()

        # Run participants
        with (
            open(fluid_log_path, "w") as fluid_log,
            open(particle_log_path, "w") as particle_log,
        ):
            fluid_process = subprocess.Popen(
                [
                    str(executable_path),
                    "Fluid",
                    str(solver_config_path),
                    str(precice_config_path),
                ],
                cwd=path,
                stdout=fluid_log,
                stderr=fluid_log,
            )
            particle_process = subprocess.Popen(
                [
                    str(executable_path),
                    "Particle",
                    str(solver_config_path),
                    str(precice_config_path),
                ],
                cwd=path,
                stdout=particle_log,
                stderr=particle_log,
            )

            # Wait for both processes to finish and check return codes
            fluid_process.wait()
            if fluid_process.returncode != 0:
                raise RuntimeError(
                    f"Fluid process failed with code {fluid_process.returncode}"
                )

            particle_process.wait()
            if particle_process.returncode != 0:
                raise RuntimeError(
                    f"Particle process failed with code {particle_process.returncode}"
                )

            print(f"Finished {path.relative_to(script_path.parent)}")

    except Exception as e:
        print(f"Error in {path.relative_to(script_path.parent)}: {e}")


def find_runs(path: pathlib.Path) -> list[pathlib.Path]:
    """Recursively finds all runs in the given directory.
    Returns a list of paths to the individual runs."""
    runs = []

    for path in path.iterdir():
        if path.name in ["__pycache__"]:
            continue

        if (path / "precice-config.xml").is_file() and (
            path / "parameters.prm"
        ).is_file():
            runs.append(path)

        elif path.is_dir():
            runs.extend(find_runs(path))

    return runs


def run_all(path: pathlib.Path, threads: int = 1):
    """Recursively find all runs in the given directory and execute them.

    If threads is greater than 1, runs are executed in parallel using the given
    number of threads."""

    run_paths = find_runs(path)

    if threads > 1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            executor.map(run, run_paths)
    else:
        for path in run_paths:
            run(path)
