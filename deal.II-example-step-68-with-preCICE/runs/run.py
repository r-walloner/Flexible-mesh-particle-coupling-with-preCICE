import sys
import pathlib

# Add top-level directory to path
script_path = pathlib.Path(__file__).resolve()
sys.path.append(str(script_path.parents[1]))
from runs import run_all  # noqa: E402

run_all(script_path.parent, threads=12)
