from pathlib import Path
import pyvista as pv
import numpy as np

# Create output directory
out_dir = Path(__file__).resolve().parent / "meshes"
out_dir.mkdir(parents=True, exist_ok=True)

# Define grid dimensions
nx, ny, nz = 10, 10, 10  # Number of points in x, y, z directions

# Define the number of time steps
timesteps = 10000

# Define the physical dimensions of the grid
x = np.linspace(0, 0.004, nx)  # 0 to 0.004 in x-direction
y = np.linspace(0, 0.004, ny)  # 0 to 0.004 in y-direction
z = np.linspace(0, 0.006, nz)  # 0 to 0.006 in z-direction

# Create a structured grid
x, y, z = np.meshgrid(x, y, z)

points = np.vstack((x.ravel(), y.ravel(), z.ravel())).T
cells = []
for i in range(nx - 1):
    for j in range(ny - 1):
        for k in range(nz - 1):
            # Define the vertices of each hexahedron
            v0 = i * ny * nz + j * nz + k
            v1 = v0 + 1
            v2 = v0 + nz
            v3 = v0 + nz + 1
            v4 = (i + 1) * ny * nz + j * nz + k
            v5 = v4 + 1
            v6 = v4 + nz
            v7 = v4 + nz + 1
            cells.append([8, v0, v1, v3, v2, v4, v5, v7, v6])

celltypes = [pv.CellType.HEXAHEDRON] * (nx - 1) * (ny - 1) * (nz - 1)

# Create the grid object
grid = pv.UnstructuredGrid(cells, celltypes, points)

for i in range(timesteps):
    # Add a vector field to each point
    vectors = np.zeros((nx * ny * nz, 3))  # Initialize a zero vector field
    vectors[:, 0] = 0
    vectors[:, 1] = 0
    vectors[:, 2] = 1

    # Attach the vector field to the grid
    grid.point_data["Force"] = vectors

    # Save the grid to a VTK file
    out_path = out_dir / f"Fluid-Mesh.dt{str(i)}.vtu"
    print(f"Writing {out_path}")
    grid.save(out_path.as_posix())
