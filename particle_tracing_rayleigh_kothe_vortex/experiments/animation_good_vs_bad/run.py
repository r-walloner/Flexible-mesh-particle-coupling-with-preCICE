import sys
import pathlib

script_dir = pathlib.Path(__file__).resolve().parent
sys.path.append(str(script_dir.parent))
import util  # noqa: E402

util.run_all(script_dir, threads=12)
