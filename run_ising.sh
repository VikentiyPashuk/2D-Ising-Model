#!/bin/bash

for size in 5 10 20 50 100 125 250
do
    for mcs in 50 500 2000
    do
        echo "Running the 2D Ising Model with size of $size and mcs of $mcs"
        time ./ising-2d "$size" "$mcs" > "ising_${size}_${mcs}.dat"
        echo "The results are in ising_${size}_${mcs}.dat"
    done
done

echo "Writing the Initial and Final Matrices for the 2D Ising Model"
time ./matrix_ising
echo "The resutls are in spins_final_matrix.txt and spins_initial_matrix.txt"

echo "Computing the algorithm for Critical Temperature for 5000 MCS"
echo "Strap in because this will take a while"
time ./Tc_ising  > tc_ising.dat
echo "The resutls are in tc_ising.dat <Tc Magnetic Sus.> <Tc Heat Cap.> <Lattice Size>"



