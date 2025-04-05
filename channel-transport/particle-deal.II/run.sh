#!/bin/sh

mkdir -p ../solution

mpirun -np 6 ../build/particle-deal.II
