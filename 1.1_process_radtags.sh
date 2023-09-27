#!/bin/sh
# Grid Engine options
#$ -N process_radtags
#$ -cwd
#$ -l h_rt=4:00:00
#$ -l h_vmem=4G
#$ -R y
#$ -t 1-193
#$ -e e_files
#$ -o o_files
#$ -m beas

# Jobscript to remove barcodes and low-quality reads from DArTseq libraries

# Initialise the Modules environment
. /etc/profile.d/modules.sh

# Load stacks and set up vars

module load roslin/stacks/2.5.4

INPUT_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/data/raw
BARCODES=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/file_lists/barcode_file.txt
OUTPUT_DIR=/exports/cmvm/eddie/eb/groups/ogden_grp/marc/wbd_dart_2022/data/out/barcode_removal

# Run process_radtags

SAMPLE=`ls $INPUT_DIR/*.FASTQ.gz`

THIS_SAMPLE=$(echo "${SAMPLE}" | sed -n ${SGE_TASK_ID}p)

mkdir $OUTPUT_DIR/${SGE_TASK_ID}

process_radtags \
        -f $THIS_SAMPLE \
        -b $BARCODES \
        -o $OUTPUT_DIR/${SGE_TASK_ID} \
        --inline-null \
        --renz_1 pstI \
        --renz_2 sphI \
        -c \
        -q \
        -r \
        --filter_illumina \
        -s 20 \
        -D \

# -c remove reads with uncalled bases
# -q discard reads with low quality scores
# -r rescue barcodes and radtags