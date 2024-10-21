#!/bin/bash -l

#$ -N Block2ZCTA      ## Name the job
#$ -j y               ## Merge error & output files
#$ -pe omp 8          ## Request between 16 GB and 32 GB of RAM

module load R/4.2.1
Rscript 02a_Creating_Block_to_ZCTA_Crosswalk.R $SGE_TASK_ID 

## NOTE: if R script is not saved in the same directory as the bash script, enter
##       the full directory pathway to the R file in the line beginning with "Rscript"
##
## In Terminal, cd to the directory in which this bash script is located. 
##
## Then, submit the job as an array *for each state* with this command in Terminal:
##
## qsub -P projectName -t 1-51 02b_Creating_Block_to_ZCTA_Crosswalk_BSGEscript.sh