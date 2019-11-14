#!/bin/bash
#SBATCH -N 2
#SBATCH -J kaon
###SBATCH --mail-user=yz3210@columbia.edu
###SBATCH --mail-type=ALL
#SBATCH -t 00:30:00
###SBATCH -t 6:00:00
#SBATCH -o output-%j.txt

# source /hpcgpfs01/software/intel/bin/compilervars.sh -arch intel64
module list

#OpenMP settings:
export OMP_NUM_THREADS=256
### export OMP_NUM_THREADS=128
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


#run the application:
# srun -n 2 -c 64 --cpu_bind=cores ./convert_wall --mpi 1.1.1.2 --decomposition 
# srun -n 2 -c 64 --cpu_bind=cores ./convert_point --mpi 1.1.1.2 --decomposition 
# srun -n 2 -c 64 --cpu_bind=cores ./test2 --mpi 1.1.1.2 --decomposition 
# srun -n 2 -c 64 --cpu_bind=cores ./test --mpi 1.1.1.2 --decomposition 
mpirun -n 2 ./test --mpi 1.1.1.2 --decomposition 

