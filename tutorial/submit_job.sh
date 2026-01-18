#########################################################################
#    SPF - Stochastic Phase Field
#    Copyright (C) 2025 
#    Peng Geng <penggeng@g.ucla.edu>
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#    See the README file in the top-level SPF directory.
#########################################################################
# File: submit_job.sh
# Purpose: This script is an example of an HPC submit script.

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
