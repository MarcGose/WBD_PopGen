#!/bin/sh
# Grid Engine options
#$ -N angsd_gl
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 4
#$ -o o_files
#$ -e e_files
#$ -o xtrace

# Jobscript to:
# Filter SNPs by read depth, callrate and Maf
# Genotype calling
# calculate allele frequency

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Dependencies and setup env variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Initialise the Modules environment
. /etc/profile.d/modules.sh

# Load java
module add java

SOFTWARE_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/emily/software/angsd-v0.935-53
BAM_LIST=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/file_lists/bamlist.txt
OUTPUT_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/data/out/angsd
REFERENCE=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/assembly/WBD/ncbi_dataset/data/GCA_949774975.1/GCA_949774975.1_mLagAlb1.1_genomic.fna
SNP_LIST=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/file_lists/snp_list_tab.txt

# Run process

# Genotype likelihoods and SNP / allele frequencies

$SOFTWARE_DIR/angsd \
        -b $BAM_LIST \
        -ref $REFERENCE \
        -P 4 \
        -out $OUTPUT_DIR/WBD_GLs \
        -uniqueOnly 1 \
        -remove_bads 1 \
        -only_proper_pairs 1 \
        -trim 0 \
	-minInd 126 \
        -C 50 \
        -baq 1 \
        -minMapQ 30 \
        -minQ 30 \
        -doCounts 1 \
        -GL 2 \
        -doGlf 2 \
        -doMajorMinor 4 \
        -doMaf 1 \
        -SNP_pval 1e-6 \
	-setMinDepth 785 \
	-setMaxDepth 3140 \
	-sites $SNP_LIST \
        -doGeno 2 \
        -doPost 1 \
        -doPlink 2
