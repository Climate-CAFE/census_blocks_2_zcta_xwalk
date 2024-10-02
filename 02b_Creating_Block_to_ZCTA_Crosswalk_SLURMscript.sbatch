#!/bin/bash -l

#SBATCH --job-name=Block2ZCTA      ## Name the job
#SBATCH --output=Block2ZCTA_%A.out ## Output file
#SBATCH --error=Block2ZCTA_%A.err  ## Error file
#SBATCH --cpus-per-task=8          ## Request 8 CPU cores
#SBATCH --mem=16G                  ## Request memory (adjust as needed)
#SBATCH --array=1-51               ## Submit as a job array with 51 tasks

module load R/4.2.1
Rscript 02a_Creating_Block_to_ZCTA_Crosswalk.R $SLURM_ARRAY_TASK_ID

## NOTE: if R script is not saved in the same directory as the bash script, enter
##       the full directory pathway to the R file in the line beginning with "Rscript"
##
## In Terminal, cd to the directory in which this bash script is located. 
##
## Then, submit the job with this command in Terminal:
##
## sbatch 02b_Creating_Block_to_ZCTA_Crosswalk_SLURMscript.sbatch