def generate(p):
    if p["particle_total_subdomains"] is None:
        subdomains = p["particle_subdomains"].split(" ")
        p["particle_total_subdomains"] = int(subdomains[0]) * int(subdomains[1]) * int(subdomains[2])

    return f"""#!/bin/bash

{"srun" if p["slurm"] else "mpirun"} -n {p["particle_total_subdomains"]} ../../../../../LIGGGHTS-PUBLIC/build/liggghts < in.liggghts
"""