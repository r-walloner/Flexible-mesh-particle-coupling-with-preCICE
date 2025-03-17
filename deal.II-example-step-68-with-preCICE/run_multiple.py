import pathlib
import subprocess

script_path = pathlib.Path(__file__).resolve()

executable_path = script_path.parent / "build" / "step-68-with-precice"
if not executable_path.is_file():
    raise FileNotFoundError(f"Executable {executable_path} not found")


def run(path: pathlib.Path):
    """Execute one run with the config files in the given directory."""

    print(f"Running {path}")

    # Check if config files exist
    solver_config_path = path / "parameters.prm"
    if not solver_config_path.is_file():
        raise FileNotFoundError(f"Solver config {solver_config_path} not found")

    precice_config_path = path / "precice-config.xml"
    if not precice_config_path.is_file():
        raise FileNotFoundError(f"Solver config {precice_config_path} not found")

    # Create log files
    fluid_log_path = path / "fluid.log"
    if fluid_log_path.exists():
        raise FileExistsError(f"Log file {fluid_log_path} already exists")
    fluid_log_path.touch()

    particle_log_path = path / "particle.log"
    if particle_log_path.exists():
        raise FileExistsError(f"Log file {particle_log_path} already exists")
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

    print(f"Finished {path}\n")


def run_all(path: pathlib.Path):
    """Recursively find all runs in the given directory and execute them."""
    for path in path.iterdir():
        if path.name in ["__pycache__"]:
            continue

        if (path / "precice-config.xml").is_file() and (
            path / "parameters.prm"
        ).is_file():
            try:
                run(path)
            except Exception as e:
                print(f"Error in {path}. Skipping.")
                print(e)
                continue

        elif path.is_dir():
            run_all(path)
