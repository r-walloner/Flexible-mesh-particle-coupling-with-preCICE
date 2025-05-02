#!/bin/sh

mkdir -p ../solution

mpirun -np 2 ../build/particle-deal.II
