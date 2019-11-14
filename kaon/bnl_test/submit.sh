#!/bin/bash
#SBATCH -N 4
#SBATCH -J kaon
###SBATCH --mail-user=yz3210@columbia.edu
###SBATCH --mail-type=ALL
#SBATCH -t 00:30:00
###SBATCH -t 06:00:00
#SBATCH -o output-%j.txt

source /hpcgpfs01/software/intel/bin/compilervars.sh -arch intel64

#OpenMP settings:
export OMP_NUM_THREADS=256
### export OMP_NUM_THREADS=128
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
# srun -n 4 -c 64 --cpu_bind=cores ./convert_wall --mpi 1.1.1.4 --decomposition 
srun -n 4 -c 64 --cpu_bind=cores ./test --mpi 1.1.1.4 --decomposition 

