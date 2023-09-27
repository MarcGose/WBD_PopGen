#!/bin/sh
# Grid Engine options
#$ -N bwa
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1-3
#$ -pe sharedmem 6
#$ -R y
#$ -o o_files
#$ -e e_files

# Jobscript to align reads to reference

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Dependencies and setup env variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Initialise the Modules environment
. /etc/profile.d/modules.sh

# Load blast
#module load roslin/bwa/2.1.0
module add roslin/bwa/0.7.17
module load roslin/samtools/1.9 

TARGET_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/WBD_reseq_2022/data/out/TrimGalore
OUTPUT_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/WBD_reseq_2022/data/out/bwa/mapped_bams
SAMPLE_SHEET="/exports/cmvm/eddie/eb/groups/ogden_grp/marc/WBD_reseq_2022/file_lists/trimmed_reads.txt"
REFERENCE="/exports/cmvm/eddie/eb/groups/ogden_grp/marc/assembly/WBD/ncbi_dataset/data/GCA_949774975.1/GCA_949774975.1_mLagAlb1.1_genomic.fna"

# Get list of files

base=`sed -n "$SGE_TASK_ID"p $SAMPLE_SHEET | awk '{print $1}'`
r1=`sed -n "$SGE_TASK_ID"p $SAMPLE_SHEET | awk '{print $2}'`
r2=`sed -n "$SGE_TASK_ID"p $SAMPLE_SHEET | awk '{print $3}'`


# Process and filter unmapped reads with samtools

echo Processing sample: ${base} ${r1} ${r2} on $HOSTNAME

bwa mem -t 6 $REFERENCE $TARGET_DIR/$r1 $TARGET_DIR/$r2 | samtools view -bF 4 - > $OUTPUT_DIR/${base}_mapped.bam
