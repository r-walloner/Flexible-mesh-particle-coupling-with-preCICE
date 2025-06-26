from argparse import ArgumentParser
import logging
import sys

from liggghts import liggghts
from mpi4py import MPI
from precice import Participant

# TODO get these parameters from a config file
MESH_NAME = "Fluid-Mesh"
DATA_NAME = "Force"


# Parse command line arguments
arg_parser = ArgumentParser()
arg_parser.add_argument(
    "participantName",
    type=str,
    help="Name of this participant",
)
arg_parser.add_argument(
    "liggghtsSkript",
    type=str,
    help="Path of the LIGGGHTS input script",
)
arg_parser.add_argument(
    "preciceConfig",
    type=str,
    help="Path of the preCICE configuration file",
)
try:
    args = arg_parser.parse_args()
except SystemExit:
    sys.exit(1)

# Set up MPI
mpi_comm = MPI.COMM_WORLD

# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format=f"({mpi_comm.Get_rank()}) - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)
logger.debug(f"I am process {mpi_comm.Get_rank()} of {mpi_comm.Get_size()}")
mpi_comm.Barrier()

# Initialize preCICE and LIGGGHTS
precice = Participant(
    args.participantName,
    args.preciceConfig,
    mpi_comm.Get_rank(),
    mpi_comm.Get_size(),
)

lig = liggghts()

# TODO this is hacky
lig.file(args.liggghtsSkript)

bounding_box = [
    float(lig.extract_global("subxlo", 1)),
    float(lig.extract_global("subxhi", 1)),
    float(lig.extract_global("subylo", 1)),
    float(lig.extract_global("subyhi", 1)),
    float(lig.extract_global("subzlo", 1)),
    float(lig.extract_global("subzhi", 1)),
]
precice.set_mesh_access_region(MESH_NAME, bounding_box)

precice.initialize()

# Main simulation loop
while precice.is_coupling_ongoing():
    if precice.requires_writing_checkpoint():
        # TODO write checkpoint
        raise NotImplementedError("Checkpoint writing not implemented yet")

    # Figure out the time step size
    max_dt = precice.get_max_time_step_size()
    solver_dt = lig.extract_global("dt", 0)
    dt = min(solver_dt, max_dt)
    logger.info(f"stepping dt = {dt}\n\tsolver_dt = {solver_dt}\n\tmax_dt = {max_dt}")

    # TODO do required timesteps
    n_atoms = lig.extract_global("nlocal", 0)
    logging.debug(f"local attoms = {n_atoms}")

    # Extract atom properties
    atom_positions = lig.extract_atom("x", 3)
    atom_forces = lig.extract_atom("f", 3)

    # TODO Can we do this with just one map_and_read_data call?
    for i in range(n_atoms):
        atom_position = (
            float(atom_positions[i][0]),
            float(atom_positions[i][1]),
            float(atom_positions[i][2]),
        )
        force = precice.map_and_read_data(MESH_NAME, DATA_NAME, [atom_position], dt)

        atom_forces[i * 3 + 0] += force[0]
        atom_forces[i * 3 + 1] += force[1]
        atom_forces[i * 3 + 2] += force[2]

    # TODO replace this with a smarter way of running the liggghts script
    lig.command("run 1")

    precice.advance(dt)

    if precice.requires_reading_checkpoint():
        # TODO read checkpoint
        raise NotImplementedError("Checkpoint reading not implemented yet")

logger.info("Finalizing")
lig.close()
precice.finalize()
