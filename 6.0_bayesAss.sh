#!/bin/sh
# Grid Engine options
#$ -N bayesass
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=32G
#$ -pe sharedmem 4
#$ -o o_files
#$ -e e_files
#$ -o xtrace

# Jobscript to:
# estimate contemporary migration rates

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Load Dependencies and setup env variables #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Initialise the Modules environment
. /etc/profile.d/modules.sh

# Define Variables

SOFTWARE=/exports/cmvm/eddie/eb/groups/ogden_grp/emily/software/BayesAss3-SNPs-master/BA3-SNPS.exe
INPUT_FILE=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/data/out/bayesass/final_run/WBD_GLs.inp
OUTPUT_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/data/out/bayesass/final_run

	$Software \
	--file ${INPUT_FILE} \
	--loci 1086 \
	--iterations 20000000 \
	--seed 6789 \
	--trace \
	--verbose \
	--genotypes \
	--out ${OUTPUT_DIR}/BA3_out.txt
