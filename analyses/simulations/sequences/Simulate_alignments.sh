#! /usr/bin/env bash
# Simulate alignments for different trees
# Usage: Simulate_alignments.sh <xml_folder> <output_dir> <reference_fasta> 


xml_folder=$1
output_dir=$2
reference_fasta=$3

# Check inputs
if [ ! -d $output_dir ]; then
    mkdir -p $output_dir
fi
if [ ! -d $xml_folder ]; then
    echo "XML folder $xml_folder not found!"
    exit 1
fi
if [ ! -f $reference_fasta ]; then
    echo "Reference fasta file $reference_fasta not found!"
    exit 1
fi
if [ ! -f "$(dirname $0)/rescale_tree.py" ]; then
    echo "rescale_tree.py script not found!"
    exit 1
fi

# Parameters
model="HKY{6.721}+I{0.9991}+G{0.5002}"
seq_name=$(grep ">" $reference_fasta | sed -e 's/>//;s/ .*//' -e 's/ .*$//')
group=$(basename $xml_folder)
substitution_rate=0.0000045

# Rescale trees
python "$(dirname $0)/rescale_tree.py" -f $xml_folder/*.xml -s $substitution_rate -o $output_dir/$group.rescaled.trees


# run alisim for each tree
tmp_dir=$(mktemp -d)
mkdir -p $output_dir/$group
i=0
for tree in $(cat $output_dir/$group.rescaled.trees); do
    echo $tree > $tmp_dir/tree.nwk
    iqtree3 -T 8 --alisim $output_dir/$group/tree_${i} -m $model -t $tmp_dir/tree.nwk --length $("$(dirname "$0")/../../../src/alnlen" -i $reference_fasta) --redo --out-format fasta --root-seq $reference_fasta,$seq_name
    ((i++))
done
rm -rf $tmp_dir

