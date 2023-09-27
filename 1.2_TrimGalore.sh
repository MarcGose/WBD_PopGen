#!/bin/sh
# Grid Engine options
#$ -N trim_galore
#$ -cwd
#$ -l h_rt=16:00:00
#$ -l h_vmem=12G
#$ -R y
#$ -e e_files
#$ -o o_files
#$ -t 1-3
#$ -pe sharedmem 4

# Jobscript to run trim_galore

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Dependencies and setup env variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#initialise modules
. /etc/profile.d/modules.sh

# Load modules
module load igmm/apps/TrimGalore/0.6.6
module load igmm/apps/FastQC/0.11.9
module load igmm/apps/pigz/2.3.3
module load igmm/apps/cutadapt/1.16

TARGET_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/WBD_reseq_2022/data/raw/low_pass
SAMPLE_SHEET="/exports/cmvm/eddie/eb/groups/ogden_grp/marc/WBD_reseq_2022/file_lists/fastq_list.txt"
OUPUT_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/WBD_reseq_2022/data/out/TrimGalore

# Get list of files
base=`sed -n "$SGE_TASK_ID"p $SAMPLE_SHEET | awk '{print $1}'`
r1=`sed -n "$SGE_TASK_ID"p $SAMPLE_SHEET | awk '{print $2}'`
r2=`sed -n "$SGE_TASK_ID"p $SAMPLE_SHEET | awk '{print $3}'`


# Process
echo Processing sample: ${base} ${r1} ${r2}

trim_galore --fastqc \
	-q 30 \
	-j 4 \
	--length 35 \
	--paired \
	$TARGET_DIR/$r1 $TARGET_DIR/$r2



