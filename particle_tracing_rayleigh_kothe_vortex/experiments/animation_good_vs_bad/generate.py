import pathlib
import sys

script_dir = pathlib.Path(__file__).resolve().parent
sys.path.append(str(script_dir.parent))
import util  # noqa: E402

util.generate_run(
    script_dir / "good",
    refinement=7,
    mapping="rbf-pum-direct",
    basis_function="compact-polynomial-c6",
    support_radius=0.5,
    time_step=1e-3,
    final_time=4.0,
    output_interval=2,
)

util.generate_run(
    script_dir / "bad",
    refinement=3,
    mapping="nearest-neighbor",
    time_step=1e-3,
    final_time=4.0,
    output_interval=2,
)
