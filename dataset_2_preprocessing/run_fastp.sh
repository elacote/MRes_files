#!/bin/bash
# Runs fastp in PE mode for all .fastq.gz files in a directory
# Run as:
# ./run_fastp.sh directory_of_samples 

timestamp=`date "+%Y%m%d-%H%M%S"`
logfile="run_$timestamp.log"
exec > $logfile 2>&1  #all output will be logged to logfile

FASTP_EXEC="/d/in7/s/fastp/fastp"
DIR=$1
ADAPTER_FILE="./adapters_all.fa"

echo "Running fastp using executable: $FASTP_EXEC"

for file in `ls $DIR/*_1.fastq.gz`;
do
  sample=${file/$DIR\//}
  sample=${sample/_1.fastq.gz/}
  echo "Sample= $sample"
  $FASTP_EXEC -i  "$DIR/$sample"_1.fastq.gz  \
              -I  "$DIR/$sample"_2.fastq.gz \
         -o "$sample"_1_trimmed.fastq.gz \
         -O "$sample"_2_trimmed.fastq.gz \
         --adapter_fasta $ADAPTER_FILE \
         -l 25 --trim_poly_x \
         -h "$sample"_fastp.html -j "$sample"_fastp.json
done
