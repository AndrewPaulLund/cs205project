#!/bin/bash
#SBATCH -p mpi #partition
#SBATCH -t 0-24:01 #time days-hr:min
#SBATCH -n 5 #number of tasks
#SBATCH --mem=8G #memory per job (all cores), GB
#SBATCH -o %j.out #out file
#SBATCH -e %j.err #error file

source activate cs205
#module load openmpi/2.0.1

#python /home/kt184/cs205/mpi_examples/mpi4py-examples/01-hello-world
#mpiexec -n 50 python /home/kt184/cs205/mpi_examples/mpi4py-examples/03-scatter-gather
#mpiexec -n 50 python /home/kt184/cs205/mpi_examples/mpi4py-examples/09-task-pull.py

extension=/home/kt184/cs205/mpi_examples/mpi4py-examples/results/RNA1.core4
/home/kt184/.conda/envs/cs205/bin/mpirun -n 5 python ./runMPIpileup.py 1000000 /n/scratch2/kt184/data/RNA/HG00096.1.M_111124_6.bam /n/scratch2/kt184/data/genome/KT_package/hg19_kt/hg19.fa 5 $extension chr1 249250621

