import numpy as np
import pyvista as pv
import pathlib


solution_dir = pathlib.Path("./solution")

# Iterate over all solution subdirectories
for path in solution_dir.iterdir():
    if not path.is_dir():
        continue

    print(f"Postprocessing {path}")

    # Iterate over all time steps
    for file in path.glob("*_particles_*.pvtu"):
        if not file.is_file():
            continue

        mesh = pv.read(file)

        # Compute error between analytical and coupled solution
        mesh["location_error"] = mesh["analytical_location"] - mesh.points
        mesh.field_data["location_error_mean_L2"] = np.mean(
            np.linalg.norm(mesh["location_error"], axis=1)
        )

        output_file = solution_dir / "postprocessed" / f"{file.stem}.vtu"
        output_file.parent.mkdir(exist_ok=True)
        mesh.save(output_file)
