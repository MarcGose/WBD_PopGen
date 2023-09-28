#!/bin/sh
# Grid Engine options
#$ -N parallel_stru
#$ -cwd
#$ -l h_rt=10:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 16
#$ -R y
#$ -e e_files/
#$ -o o_files/
#$ -t 1-8

# Jobscript to run Parallel Structure

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Dependencies and setup env variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Initialise the Modules environment
. /etc/profile.d/modules.sh

# Load R

module load R

# Run

infile="/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/data/out/structure/WBD_GLs.stru"

echo ${infile}

Rscript scripts/run_structure.R ${SGE_TASK_ID} ${infile}
