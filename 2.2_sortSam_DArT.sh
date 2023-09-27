#!/bin/sh
# Grid Engine options
#$ -N sort_sam
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -R y
#$ -t 1-193
#$ -o o_files
#$ -e e_files
#$ -m beas

# Jobscript to sort bam files

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Dependencies and setup env variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Initialise the Modules environment
. /etc/profile.d/modules.sh

# Load java
module add java 

TARGET_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/data/out/bwa_ref/mapped
OUTPUT_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/data/out/bwa_ref/sorted
SCRATCH=/exports/eddie/scratch/s2113685

# Get list of files in target directory

BAM=$(ls -1 ${TARGET_DIR}/*mapped.bam)

# Get file to be processed by *this* task
# Extract the Nth file in the list of files, $bam, where N == $SGE_TASK_ID

THIS_BAM=$(echo "${BAM}" | sed -n ${SGE_TASK_ID}p)
echo Processing file: ${THIS_BAM} on $HOSTNAME

BASE=$(echo "$THIS_BAM" | cut -f 14 -d '/' | cut -f 1 -d '.')

echo Saving file ${THIS_BAM} as $SCRATCH/${BASE%.bam}_sorted.bam

java -Xmx16g -jar /exports/cmvm/eddie/eb/groups/ogden_grp/software/picard/picard.jar SortSam \
	I=${THIS_BAM} \
	O=${OUTPUT_DIR}/${BASE%.bam}_sorted.bam \
	SORT_ORDER=coordinate \
	TMP_DIR=${SCRATCH}
