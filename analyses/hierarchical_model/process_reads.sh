#!/bin/bash
# Script to process sequencing reads
# Jordi Sevilla Fortuny 2025
# conda create -n process_reads -c conda-forge -c bioconda sra-tools==3.2.1 fastqc==0.12.1 snippy==4.6.0 shovill==1.1.0

# Check if all tools are installed
if ! command -v fasterq-dump &> /dev/null; then
    echo "fasterq-dump could not be found, please install it."
    exit 1
fi
if ! command -v fastqc &> /dev/null; then
    echo "FastQC could not be found, please install it."
    exit 1
fi
if ! command -v snippy &> /dev/null; then
    echo "snippy could not be found, please install it."
    exit 1
fi
if ! command -v shovill &> /dev/null; then
    echo "shovill could not be found, please install it."
    exit 1
fi 



logthis () { #Log things
    echo $(date) "|" $1
}

# Read parameters using optargs: accession ID, reads path, mapping path, assembly path, reports path, reference genome
while getopts "a:r:m:s:p:f:" opt; do
  case $opt in
    a) ACCESSION_ID="$OPTARG" ;;
    r) READS_PATH="$OPTARG" ;;
    m) MAPPING_PATH="$OPTARG" ;;
    s) ASSEMBLY_PATH="$OPTARG" ;;
    p) REPORTS_PATH="$OPTARG" ;;
    f) REFERENCE="$OPTARG" ;;
    *) echo "Invalid option: -$OPTARG" >&2 ;;
    esac
done
shift $((OPTIND -1))    

# Log the parameters:
echo "Accession ID: $ACCESSION_ID"
echo "Reads Path: $READS_PATH"
echo "Mapping Path: $MAPPING_PATH"
echo "Assembly Path: $ASSEMBLY_PATH"
echo "Reports Path: $REPORTS_PATH"
echo "Reference Genome: $REFERENCE"

# Create directories if they do not exist
mkdir -p $READS_PATH
mkdir -p $MAPPING_PATH
mkdir -p $ASSEMBLY_PATH
mkdir -p $REPORTS_PATH

# Step 1: Download reads using fasterq-dump
logthis "Downloading reads for accession ID: $ACCESSION_ID"
fasterq-dump $ACCESSION_ID --outdir $READS_PATH --split-files -p --threads 6
if [ $? -ne 0 ]; then
    logthis "Error downloading reads for $ACCESSION_ID"
    exit 1
fi
logthis "Successfully downloaded reads for $ACCESSION_ID"

## Compress the downloaded fastq files
logthis "Compressing downloaded reads"
gzip $READS_PATH/${ACCESSION_ID}_1.fastq
gzip $READS_PATH/${ACCESSION_ID}_2.fastq

# Step 2: Quality check using FastQC
logthis "Running FastQC on downloaded reads"
fastqc $READS_PATH/${ACCESSION_ID}_1.fastq.gz $READS_PATH/${ACCESSION_ID}_2.fastq.gz -outdir $REPORTS_PATH --threads 6
if [ $? -ne 0 ]; then
    logthis "Error running FastQC for $ACCESSION_ID"
    exit 1
fi
logthis "FastQC completed for $ACCESSION_ID"

# Step 3: Mapping reads to reference genome using snippy
logthis "Mapping reads to reference genome using snippy"
snippy --cpus 6 --outdir $MAPPING_PATH/$ACCESSION_ID \
    --ref $REFERENCE --R1 $READS_PATH/${ACCESSION_ID}_1.fastq.gz --R2 $READS_PATH/${ACCESSION_ID}_2.fastq.gz \
    --mincov 10 --minfrac 0.9 --cleanup

if [ $? -ne 0 ]; then
    logthis "Error mapping reads for $ACCESSION_ID"
    exit 1
fi
logthis "Successfully mapped reads for $ACCESSION_ID"

# Step 4: De novo assembly using shovill
logthis "Assembling reads using shovill"
shovill --R1 $READS_PATH/${ACCESSION_ID}_1.fastq.gz --R2 $READS_PATH/${ACCESSION_ID}_2.fastq.gz \
    --outdir $ASSEMBLY_PATH/$ACCESSION_ID --force 

if [ $? -ne 0 ]; then
    logthis "Error assembling reads for $ACCESSION_ID"
    exit 1
fi
mv $ASSEMBLY_PATH/$ACCESSION_ID/contigs.fa $ASSEMBLY_PATH/${ACCESSION_ID}_assembly.fasta
rm -rf $ASSEMBLY_PATH/$ACCESSION_ID
logthis "Successfully assembled reads for $ACCESSION_ID"

logthis "Pipeline completed successfully for $ACCESSION_ID"
exit 0