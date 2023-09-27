#!/bin/sh
# Grid Engine options
#$ -N ngsadmix
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=4G
#$ -pe sharedmem 4
#$ -o o_files
#$ -e e_files
#$ -t 1-8

# Jobscript to NGSadmix on genotype likelihoods

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Dependencies and setup env variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Initialise the Modules environment
. /etc/profile.d/modules.sh

SOFTWARE_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/software/angsd/misc/
TARGET_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/data/out/angsd
OUTPUT_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/data/out/angsd/ngsadmix

for i in 1 2 3 4 5 6 7 8 9 10
do
        $SOFTWARE_DIR/NGSadmix \
        -likes $TARGET_DIR/WBD_GLs.beagle.gz \
        -K ${SGE_TASK_ID} \
	-minMaf 0.01 \
        -P 4 \
        -seed $i \
        -o $OUTPUT_DIR/WBD_GLs_NGSadmix_K${SGE_TASK_ID}_run${i}_out
done


