#!/bin/bash
NUM_PROCS=$1
#SBATCH -p main
#SBATCH -n $NUM_PROCS
mpic++ lab_3.cpp -o lab_3.o
mpirun -n $NUM_PROCS ./lab_3.o