#!/bin/bash
#SBATCH -p mpi #partition
#SBATCH -t 0-24:01 #time days-hr:min
#SBATCH -n 3 #number of tasks
#SBATCH --mem=8G #memory per job (all cores), GB
#SBATCH -o %j.out #out file
#SBATCH -e %j.err #error file

module load conda2/4.2.13
source activate cs205
#module load openmpi/2.0.1

#python /home/kt184/cs205/mpi_examples/mpi4py-examples/01-hello-world
#mpiexec -n 50 python /home/kt184/cs205/mpi_examples/mpi4py-examples/03-scatter-gather
#mpiexec -n 50 python /home/kt184/cs205/mpi_examples/mpi4py-examples/09-task-pull.py

extension=/home/kt184/cs205/mpi_examples/mpi4py-examples/results/DNA1.core2
/home/kt184/.conda/envs/cs205/bin/mpirun -n 3 python ./runMPIpileup.py 1000000 /n/scratch2/kt184/data/DNA/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam /home/kt184/scratch/data/genome/KT_package/icgc_genome_broadvariant/Homo_sapiens_assembly19.fasta 3 $extension 1 249250621

