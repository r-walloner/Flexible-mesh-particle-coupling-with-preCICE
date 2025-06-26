from sys import argv
from mpi4py import MPI
from liggghts import liggghts

mpi_comm = MPI.COMM_WORLD
print(f"Process {mpi_comm.Get_rank()} of {mpi_comm.Get_size()} running")

lig = liggghts()

in_file = argv[1]

lig.file(in_file)

n_atoms = lig.extract_global("nlocal", 0)

atom_ids = lig.extract_atom("id", 0)
atom_id_list = [atom_ids[i] for i in range(n_atoms)]

atom_radiuses = lig.extract_atom("radius", 2)
atom_radiuses[0] = 0.0001

lig.command("run 1")



lig.close()

MPI.Finalize()
