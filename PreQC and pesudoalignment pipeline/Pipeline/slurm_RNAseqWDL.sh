#!/bin/bash
#SBATCH -p skylake
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=1:00:00
#SBATCH --mem=20GB

# Notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=f.tang@adelaide.edu.au

module load arch/arch/haswell
module load fastqc
module load kallisto/0.43.1-foss-2017a
module load Java
module load Python

java -jar cromwell-54.jar run RNAseq.wdl --inputs ${1}
