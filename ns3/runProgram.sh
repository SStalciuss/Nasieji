#!/bin/bash
#SBATCH -p main
#SBATCH -n64
module load openmpi
mpic++ lab_3.cpp -o lab_3.o
mpirun lab_3.o
