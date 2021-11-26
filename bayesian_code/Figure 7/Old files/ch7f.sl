#!/bin/bash -e
#SBATCH --job-name=rishikas_job
## Your NeSI project:
#SBATCH --account=uoa03374
## Maximum time for each individual job on each individual core. If
## your job runs for longer it will be TERMINATED. I am asking for ten
## hours below (HH:MM:SS).
#SBATCH --time=8:00:00
## Amount of memory (RAM) for each individual job on each core.
#SBATCH --mem-per-cpu=20G
## You can include the below two lines if you want emails when your
## jobs start/end.
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rcho900@aucklanduni.ac.nz
## Code for array job.
module load R/3.6.2-gimkl-2020a
## Running an R script on each job.
srun R --vanilla < ch7f.R
