#!/bin/bash
#SBATCH -N 1
#SBATCH -J kaon
###SBATCH --mail-user=yz3210@columbia.edu
###SBATCH --mail-type=ALL
#SBATCH -t 00:30:00
#SBATCH -o output-%j.txt
#SBARCH -p debug

#OpenMP settings:
export OMP_NUM_THREADS=256
### export OMP_NUM_THREADS=128
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
# srun -n 1 -c 64 --cpu_bind=cores ./convert_wall --mpi 1.1.1.1
mpirun -n 1 ./typeIII --mpi 1.1.1.2 --decomposition 

