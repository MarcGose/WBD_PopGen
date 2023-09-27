#!/bin/sh
# Grid Engine options
#$ -N add_RG
#$ -cwd
#$ -t 1-193
#$ -l h_rt=16:00:00
#$ -l h_vmem=6G
#$ -pe sharedmem 4
#$ -o o_files
#$ -e e_files

# Jobscript to add read groups

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Dependencies and setup env variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Initialise the Modules environment
. /etc/profile.d/modules.sh

# Load java
module add java 

SCRATCH=/exports/eddie/scratch/s2113685
INPUT_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/data/out/bwa_ref/sorted/
OUTPUT_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/data/out/bwa_ref/rg_bams/

# Get list of files in target directory

bam=$(ls -1 ${INPUT_DIR}/*_sorted.bam)

# Get file to be processed by *this* task
# Extract the Nth file in the list of files, $bam, where N == $SGE_TASK_ID

this_bam=$(echo "${bam}" | sed -n ${SGE_TASK_ID}p)
base=$(echo "${bam}" | sed -n ${SGE_TASK_ID}p | cut -f 13 -d '/' | cut -f 1 -d '.')
echo Processing file: ${this_bam} on $HOSTNAME

java -Xmx4g -jar /exports/cmvm/eddie/eb/groups/ogden_grp/software/picard/picard.jar AddOrReplaceReadGroups \
	I=$this_bam \
       	O=${this_bam%.bam}_RG.bam \
       	RGID=$base \
        RGPL=illumina \
        RGLB=$base \
        RGPU=$base \
        RGSM=$base \
        VALIDATION_STRINGENCY=SILENT \
        SORT_ORDER=coordinate \
        TMP_DIR=$OUTPUT_DIR

