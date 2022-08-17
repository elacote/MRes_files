#!/bin/bash

# runs kallisto alignment-free mapping to homo sapiens cDNA+ncRNA_norRNA tanscripts 

# Run as:
# run_kallisto_unstranded.sh directory_of_fastq_files samples

timestamp=`date "+%Y%m%d-%H%M%S"`
logfile="run_$timestamp.log"
exec > $logfile 2>&1



dir=$1
shift

#set location of executables
KALLISTO_EXEC="/d/in7/s/kallisto/kallisto_linux-v0.46.1/kallisto"

#set parameters for kallisto 
numProc=8   #number of processors/threads to be used
transcriptomeIndex="/d/in18/u/we001/kallisto_norRNA/with_gencode/gencode.v34.transcripts_norRNA.index" #index must be created ahead of the run
numBootstraps=100 #how many times to bootstrap the sample

for sample in "$@";
do
   echo "Running kallisto on sample $sample paired reads..."
   #this needs to be adapted according to how the sample files are named

   r1PairedFile="$dir$sample"_1_trimmed.fastq.gz
   r2PairedFile="$dir$sample"_2_trimmed.fastq.gz


   $KALLISTO_EXEC quant \
          -i $transcriptomeIndex \
          -t $numProc \
          -b $numBootstraps \
          --fusion \
          -o $sample \
          $r1PairedFile \
          $r2PairedFile

done
