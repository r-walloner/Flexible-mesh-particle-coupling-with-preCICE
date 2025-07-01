def generate(p):
    subdomains = p["particle_subdomains"].split(" ")
    total_subdomains = int(subdomains[0]) * int(subdomains[1]) * int(subdomains[2])

    return f"""#!/bin/bash

{"srun" if p["slurm"] else "mpirun"} -n {total_subdomains} ../../../../../LIGGGHTS-PUBLIC/build/liggghts < in.liggghts
"""