#!/bin/bash
# Script to process sequencing reads
# Jordi Sevilla Fortuny 2025
# Check if all tools are installed
# conda create -n phylo -c conda-forge -c bioconda iqtree==2.4.0 gubbins==3.4.3 snp-sites==2.5.1 seqtk==1.5

if ! command -v iqtree &> /dev/null; then
    echo "fasterq-dump could not be found, please install it."
    exit 1
fi
if ! command -v snp-sites &> /dev/null; then
    echo "FastQC could not be found, please install it."
    exit 1
fi
if ! command -v run_gubbins.py &> /dev/null; then
    echo "snippy could not be found, please install it."
    exit 1
fi
if ! command -v seqtk &> /dev/null; then
    echo "snippy could not be found, please install it."
    exit 1
fi


logthis () { #Log things
    echo $(date) "|" $1
}

# Read parameters using optargs: accession ID, reads path, mapping path, assembly path, reports path, reference genome
while getopts "c:f:o:" opt; do
  case $opt in
    c) CLUSTER_FILE="$OPTARG" ;;
    f) FULL_ALIGNMENT="$OPTARG" ;;
    o) OUTPUT_FOLDER="$OPTARG" ;;
    *) echo "Invalid option: -$OPTARG" >&2 ;;
    esac
done
shift $((OPTIND -1))    

# Log the parameters:
echo "Cluster file: $CLUSTER_FILE"
echo "Full_alignment: $FULL_ALIGNMENT"
echo "Output: $OUTPUT_FOLDER"

# Create directories if they do not exist
mkdir -p $OUTPUT_FOLDER

# Create a folder for storing outputs
cluster_name=$(basename $CLUSTER_FILE .clust)
outdir="$OUTPUT_FOLDER/$cluster_name"

mkdir -p  $outdir

# Step 1: Subset the cluster alignment from the full alignment
logthis "Subsetting cluster: $cluster_name"

#seqtk subseq $FULL_ALIGNMENT $CLUSTER_FILE > $outdir/${cluster_name}.alignment.fasta

if [ $? -ne 0 ]; then
    logthis "Error subsetting $alignment_name"
    exit 1
fi
logthis "Successfully subsetted $alignment_name"

# Step 2: Use gubbins
curdir=$(pwd)
cd $outdir
run_gubbins.py --prefix $cluster_name  ${cluster_name}.alignment.fasta
cd $curdir 

logthis "Run gubbins on the alignment for cluster: $cluster_name"
if [ $? -ne 0 ]; then
    logthis "Error running gubbins for cluster: $cluster_name"
    exit 1
fi
logthis "Gubbins completed for cluster: $cluster_name"

# Step 3: Infering phylogeny for  cluster
logthis "Infering phylogeny for cluster: $cluster_name"
fconst=$(snp-sites -C $outdir/${cluster_name}.alignment.fasta)
snp-sites -c $outdir/${cluster_name}.filtered_polymorphyc_sites.fasta > $outdir/${cluster_name}.gubbins.snps
iqtree -s $outdir/${cluster_name}.gubbins.snps -fconst $fconst

logthis "Phylogeny infered for cluster: $cluster_name"

logthis "Pipeline completed successfully for $ACCESSION_ID"
exit 0
