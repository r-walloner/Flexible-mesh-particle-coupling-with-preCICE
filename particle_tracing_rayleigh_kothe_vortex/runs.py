"""This file contains utility functions to automate the process of running
multiple instances (called runs) of the simulation with different parameters.
The functions in this file are used by the scripts in the `runs` directory to
generate simulation runs, run them, and postprocess the results."""

import pathlib
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor

script_path = pathlib.Path(__file__).resolve()

executable_path = script_path.parent / "build" / "step-68-with-precice"
if not executable_path.is_file():
    raise FileNotFoundError(f"Executable {executable_path} not found")


def generate_run(
    path: pathlib.Path,
    refinement=4,
    mapping="rbf-pum-direct",
    basis_function="compact_polynomial-c6",
    support_radius=0.5,
    constraint="consistent",
    method="euler_explicit",
    time_step=0.002,
    final_time=4.0,
    output_interval=10,
):
    """Generate a new run directory with the given parameters.

    This function creates the directory specified by `path` and writes the
    configuration files `precice-config.xml` and `parameters.prm` to it. The
    configuration files are generated with the given parameters."""

    if path.exists():
        print(
            f"Not generating run {path.relative_to(script_path.parent)}: Directory already exists"
        )
        return

    print(f"Generating run {path.relative_to(script_path.parent)}")

    path.mkdir(parents=True)
    (path / "solution").mkdir()

    # Write parameters.json (used by postprocessing scripts)
    with open(path / "parameters.json", "w") as file:
        json.dump(
            {
                "refinement": refinement,
                "mapping": mapping,
                "basis_function": basis_function,
                "support_radius": support_radius,
                "constraint": constraint,
                "method": method,
                "time_step": time_step,
                "final_time": final_time,
                "output_interval": output_interval,
            },
            file,
            indent=2,
        )

    # Write precice-config.xml
    with open(path / "precice-config.xml", "w") as file:
        file.write(f"""
<precice-configuration experimental="true">

  <data:vector name="Velocity" />
  <data:scalar name="Time" />

  <mesh name="Fluid-Mesh" dimensions="2">
    <use-data name="Velocity" />
    <use-data name="Time" />
  </mesh>

  <participant name="Fluid">
    <provide-mesh name="Fluid-Mesh" />
    <write-data name="Velocity" mesh="Fluid-Mesh" />
    <write-data name="Time" mesh="Fluid-Mesh" />
  </participant>

  <participant name="Particle">
    <receive-mesh name="Fluid-Mesh" from="Fluid" api-access="true" />
    <read-data name="Velocity" mesh="Fluid-Mesh" />
    <read-data name="Time" mesh="Fluid-Mesh" />
    <mapping:{mapping}
      direction="read"
      from="Fluid-Mesh"
      constraint="{constraint}">
      {f'<basis-function:{basis_function} support-radius="{support_radius}" />' if mapping == "rbf-pum-direct" else ""}
    </mapping:{mapping}>
  </participant>

  <m2n:sockets
    acceptor="Fluid"
    connector="Particle"
    exchange-directory="{path.resolve()}" />

  <coupling-scheme:serial-explicit>
    <participants first="Fluid" second="Particle" />
    <max-time value="{final_time}" />
    <time-window-size value="{time_step}" />
    <exchange data="Velocity" mesh="Fluid-Mesh" from="Fluid" to="Particle" />
    <exchange data="Time" mesh="Fluid-Mesh" from="Fluid" to="Particle" />
  </coupling-scheme:serial-explicit>

</precice-configuration>
""")

    # Write parameters.prm
    with open(path / "parameters.prm", "w") as file:
        file.write(f"""
# Listing of Parameters
# ---------------------
subsection Particle Tracking Problem
  # End time of the simulation
  set Final time                    = {final_time}

  # Refinement level of the fluid domain
  set Fluid refinement              = {refinement}
  set Output basename               = {path.name}_
  set Output directory              = {(path / "solution").resolve()}/

  # Iteration interval between which output results are written
  set Output interval               = {output_interval}

  # Refinement level of the particle domain
  set Particle grid refinement      = {refinement}

  # Refinement of the volumetric mesh used to insert the particles
  set Particle insertion refinement = 3

  # Iteration interval at which the mesh is load balanced
  set Repartition interval          = 5
  set Time step                     = {time_step}
  set Velocity degree               = 1

  set Method = {method}
end
""")


def run(path: pathlib.Path):
    """Execute the run in the given directory."""
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
    run_paths = []

    for path in sorted(path.iterdir(), key=lambda p: p.name):
        if path.name in ["__pycache__"]:
            continue

        if (path / "precice-config.xml").is_file() and (
            path / "parameters.prm"
        ).is_file():
            run_paths.append(path)

        elif path.is_dir():
            run_paths.extend(find_runs(path))

    return run_paths


def run_all(path: pathlib.Path, threads: int = 1):
    """Recursively find all runs in the given directory and execute them.

    If threads is greater than 1, runs are executed in parallel using the given
    number of threads. Else, runs are executed sequentially."""

    run_paths = find_runs(path)

    if threads > 1:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            executor.map(run, run_paths)
    else:
        for path in run_paths:
            run(path)
