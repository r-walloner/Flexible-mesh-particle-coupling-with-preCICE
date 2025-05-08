from argparse import ArgumentParser
import logging
import sys

from liggghts import liggghts
from mpi4py import MPI
from precice import Participant


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

mpi_comm = MPI.COMM_WORLD

logging.basicConfig(
    level=logging.DEBUG,
    format=f"({mpi_comm.Get_rank()}) - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)
logger.debug(f"I am process {mpi_comm.Get_rank()} of {mpi_comm.Get_size()}")
mpi_comm.Barrier()

precice = Participant(
    args.participantName,
    args.preciceConfig,
    mpi_comm.Get_rank(),
    mpi_comm.Get_size(),
)

lig = liggghts()

# TODO set bounding box

precice.initialize()

while precice.is_coupling_ongoing():
    if precice.requires_writing_checkpoint():
        # TODO write checkpoint
        raise NotImplementedError("Checkpoint writing not implemented yet")

    max_dt = precice.get_max_time_step_size()
    solver_dt = lig.extract_global("dt", 0)
    dt = min(solver_dt, max_dt)
    logger.info(
        f"stepping dt = {dt}\n",
        f"\tsolver_dt = {solver_dt}",
        f"\tmax_dt = {max_dt}",
    )

    # TODO do required timesteps

    precice.advance(dt)

    if precice.requires_reading_checkpoint():
        # TODO read checkpoint
        raise NotImplementedError("Checkpoint reading not implemented yet")

logger.info("Finalizing")    
lig.close()
precice.finalize()
