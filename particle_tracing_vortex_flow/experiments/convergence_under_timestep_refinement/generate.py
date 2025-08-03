import pathlib
import sys

script_dir = pathlib.Path(__file__).resolve().parent
sys.path.append(str(script_dir.parent))
import util  # noqa: E402

# Set constant parameters
# final_time = 0.14 # timestep 70, which we evaluate in the error plots
final_time = 4  # one full period of the velocity function
constraint = "consistent"
support_radius = 0.5  # only used for rbf

method = "euler_explicit"

for time_step in [10**-i for i in range(1, 7)]:
    output_interval = int(max(1, final_time / time_step / 100))

    for mapping in [
        "nearest-neighbor",
        # "rbf-pum-direct",
    ]:
        for refinement in range(4, 9, 2):
            # Only generate different refinements for nearest-neighbor
            if refinement != 7 and mapping != "nearest-neighbor":
                continue

            for basis_function in [
                "compact-polynomial-c0",
                # "compact-polynomial-c2",
                # "compact-polynomial-c4",
                "compact-polynomial-c6",
                # "compact-polynomial-c8",
            ]:
                path = script_dir / method / mapping
                if mapping == "rbf-pum-direct":
                    path = path / basis_function / str(support_radius)
                path = path / f"timestep-{time_step}" / f"refinement-{refinement}"

                util.generate_run(
                    path,
                    refinement=refinement,
                    mapping=mapping,
                    basis_function=basis_function,
                    support_radius=support_radius,
                    constraint=constraint,
                    method=method,
                    time_step=time_step,
                    final_time=final_time,
                    output_interval=output_interval,
                )

                if mapping == "nearest-neighbor":
                    # NN does not support different basis functions. Therefore,
                    break
