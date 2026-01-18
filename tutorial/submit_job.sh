#### submit_job.sh START ####

#!/bin/bash
#$ -cwd
#$ -o joblog.$JOB_ID              #error = Merged with joblog
#$ -j y

## Edit the line below as needed:
#$ -l h_rt=23:59:59

## Modify the parallel environment and the number of cores as needed:
#$ -pe shared 15

#$ -M $USER@g.ucla.edu            #Email address to notify
#$ -m bea                         #Notify when

# load the job environment:
. /u/local/Modules/default/init/modules.sh

module purge
module load gcc/11.3.0
module load mpich/3.4

mpirun -n 15 /path/to/spf_W-Cr/bin/spf_3d.x
