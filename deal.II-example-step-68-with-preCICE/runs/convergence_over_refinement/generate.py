import pathlib
import sys

# Add top-level directory to path
script_path = pathlib.Path(__file__).resolve()
sys.path.append(str(script_path.parents[2]))
from runs import generate_run  # noqa: E402

# Set constant parameters
constraint = "consistent"
support_radius = 0.5  # only used for rbf

# Iterate over variable parameters
for mapping in ["nearest-neighbor", "rbf-pum-direct"]:
    for refinement in range(0, 7, 1):
        for basis_function in [
            "compact-polynomial-c0",
            "compact-polynomial-c2",
            "compact-polynomial-c4",
            "compact-polynomial-c6",
            "compact-polynomial-c8",
        ]:
            path = script_path.parent / mapping
            if mapping == "rbf-pum-direct":
                path = path / basis_function / str(support_radius)
            path = path / f"refinement-{refinement}"

            generate_run(
                path,
                refinement=refinement,
                mapping=mapping,
                basis_function=basis_function,
                support_radius=support_radius,
                constraint=constraint,
            )

            if mapping == "nearest-neighbor":
                # NN does not support different basis functions. Therefore,
                break
