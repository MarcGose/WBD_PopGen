#!/bin/sh
# Grid Engine options
#$ -N IndexRef
#$ -cwd
#$ -l h_rt=1:00:00
#$ -l h_vmem=6G
#S -pe sharedmem 4
#$ -o o_files
#$ -e e_files

# Jobscript to index reference

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Dependencies and setup env variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Initialise the Modules environment
. /etc/profile.d/modules.sh

# Load bwa
module load roslin/samtools/1.9
module load roslin/bwa/0.7.17

REFERENCE=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/assembly/WBD/ncbi_dataset/data/GCA_949774975.1/GCA_949774975.1_mLagAlb1.1_genomic.fna

# index ref with  bwa for mapping SNPs

samtools faidx $REFERENCE
bwa index $REFERENCE
