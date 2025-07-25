#!/bin/sh

mkdir -p ./solution

mpirun -n 4 ../build/particle-deal.II
