#!/bin/bash

# shellcheck disable=SC1091
. ../../liggghts_adapter/.venv/bin/activate

mpirun python ../../liggghts_adapter/main.py \
    Particle \
    in.liggghts \
    ../precice-config.xml
