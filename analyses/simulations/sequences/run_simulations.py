#!/usr/bin/env python3
# This script is used to update the ReMaster master script and BDMMprime master script
# Jordi Sevilla, 31/07/2025

from argparse import ArgumentParser
from Bio import SeqIO, Phylo
import os
from io import StringIO

# Constants
RATE = 0.0000045


def read_arguments():
    parser = ArgumentParser(description="Generate xml files")
    parser.add_argument(
        "-a",
        "--aln_folder",
        type=str,
        required=True,
        help="Path to the folder containing the alignments.",
    )
    parser.add_argument(
        "-bf",
        "--bdmmprime_file",
        type=str,
        required=True,
        help="Path to the BDMMprime master script file.",
    )
    parser.add_argument(
        "-d",
        "--dir",
        help="Directory to store the output files",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-t", "--trees", help="Path to the trees file", type=str, required=True
    )

    return parser.parse_args()


def parse_alignments(aln_folder):
    alignments = []
    for filename in os.listdir(aln_folder):
        filepath = os.path.join(aln_folder, filename)
        records = list(SeqIO.parse(filepath, "fasta"))
        name = "_".join([os.path.basename(aln_folder), filename.split(".")[0]])
        alignments.append((name, records))

    return alignments


def update_bdmmprime_script_no_tiptypes(alignment, template, new_filepath, tree):
    common_alignment_line = (
        '<sequence id="dir_{}" spec="Sequence" taxon="{}" totalcount="4" value="{}"/>'
    )

    alignment_str = ""
    tiptypes = ""
    dates = ""
    for seq in alignment:
        taxon = seq.id
        height = float(seq.id.split("_")[2])

        alignment_str += common_alignment_line.format(taxon, taxon, seq.seq) + "\n"
        tiptypes += f"{taxon}=?,\n"
        dates += f"{taxon}={height},"

    tiptypes = tiptypes.rstrip(",\n")  # Remove trailing comma and newline
    dates = dates.rstrip(",")  # Remove trailing comma

    # Read tree as string
    tree_str = StringIO()
    Phylo.write(tree, tree_str, "newick")
    tree_str = tree_str.getvalue().strip()

    with open(template, "r") as file:
        content = file.read()
        content = content.replace("[ALIGNMENT]", alignment_str)
        content = content.replace("[PHENOTYPES]", tiptypes)
        content = content.replace("[DATES]", dates)
        content = content.replace("[CLOCKRATE]", str(RATE))
        content = content.replace("[TREE]", str(tree_str))

    with open(new_filepath, "w") as file:
        file.write(content)


def update_bdmmprime_script_tiptypes(alignment, template, new_filepath, tree):
    common_alignment_line = (
        '<sequence id="dir_{}" spec="Sequence" taxon="{}" totalcount="4" value="{}"/>'
    )

    alignment_str = ""
    tiptypes = ""
    dates = ""
    for seq in alignment:
        taxon = seq.id
        height = float(seq.id.split("_")[2])
        tip_type = seq.id.split("_")[3]

        alignment_str += common_alignment_line.format(taxon, taxon, seq.seq) + "\n"
        tiptypes += f"{taxon}={tip_type},\n"
        dates += f"{taxon}={height},"

    tiptypes = tiptypes.rstrip(",\n")  # Remove trailing comma and newline
    dates = dates.rstrip(",")  # Remove trailing comma

    # Read tree as string
    tree_str = StringIO()
    Phylo.write(tree, tree_str, "newick")
    tree_str = tree_str.getvalue().strip()

    with open(template, "r") as file:
        content = file.read()
        content = content.replace("[ALIGNMENT]", alignment_str)
        content = content.replace("[PHENOTYPES]", tiptypes)
        content = content.replace("[DATES]", dates)
        content = content.replace("[CLOCKRATE]", str(RATE))
        content = content.replace("[TREE]", str(tree_str))

    with open(new_filepath, "w") as file:
        file.write(content)


def process_multiple_alignments(args, alignments, trees):

    for alignment_name, records in alignments:
        new_bdmmprime_file = os.path.join(args.dir, f"{alignment_name}_no_tiptypes.xml")
        new_bdmmprime_file_tiptypes = os.path.join(args.dir, f"{alignment_name}.xml")
        tree_index = int(alignment_name.split("_")[-1].split(".")[0])
        tree = trees[tree_index]
        update_bdmmprime_script_no_tiptypes(
            records, args.bdmmprime_file, new_bdmmprime_file, tree
        )
        update_bdmmprime_script_tiptypes(
            records, args.bdmmprime_file, new_bdmmprime_file_tiptypes, tree
        )


def main():
    args = read_arguments()
    if not os.path.exists(args.dir):
        os.makedirs(args.dir, exist_ok=True)

    alignments = parse_alignments(args.aln_folder)
    trees = list(Phylo.parse(args.trees, "newick"))
    if len(alignments) != len(trees):
        print(
            f"Error: The number of alignments ({len(alignments)}) does not match the number of trees ({len(trees)})."
        )
        exit(1)
    process_multiple_alignments(args, alignments, trees)


if __name__ == "__main__":
    main()
