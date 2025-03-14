#!/bin/bash

# iterate over refinements
for refine in $(seq 2 2 6); do
    # iterate over basis functions
    for c in $(seq 0 2 8); do
        #iterate over support radius
        for radius in $(seq 0.2 0.3 0.8); do

            echo "Running rbf with refinement $refine, c $c, radius $radius"
            
            # modify precice config
            sed -i ../precice-config.xml -E \
                -e "s/(basis-function:compact-polynomial-c)[0-9]+/\1${c}/g" \
                -e "s/(support-radius=\")[^\"]+/\1${radius}/g"

            # modify solver config
            sed -i ../parameters.prm -E \
                -e "s/(refine)[0-9]+(-rbf-c)[0-9]+(-r)[0-9.]+(_)/\1${refine}\2${c}\3${radius}\4/g" \
                -e "s/(set Fluid refinement[[:space:]]*=[[:space:]]*)[0-9]+/\1${refine}/g" \
                -e "s/(set Particle grid refinement[[:space:]]*=[[:space:]]*)[0-9]+/\1${refine}/g" \

            # make solution directory
            mkdir -p ../solution/refine"$refine"-rbf-c"$c"-r"$radius"_

            # run
            ./step-68-with-precice Fluid &
            ./step-68-with-precice Particle &
            wait

        done
    done
done
