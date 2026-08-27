#! /usr/bin/env python
# Rescale a tree by a given factor (with high-precision branch lengths)
# Jordi Sevilla, 03/09/2025 – modified for precision

from argparse import ArgumentParser
from io import StringIO
import re
from decimal import Decimal, getcontext
from Bio import Phylo

# Increase precision
getcontext().prec = 20


def read_args():
    parser = ArgumentParser(
        description="Gret a time tree from a beast .xml file and rescale it by a given factor (with high-precision branch lengths)"
    )
    parser.add_argument(
        "-f",
        "--files",
        type=str,
        required=True,
        nargs="+",
        help="Beast .xml file with tree in newick tag",
    )
    parser.add_argument(
        "-s", "--sub_rate", type=float, required=True, help="Rescaling factor"
    )
    parser.add_argument(
        "-o",
        "--outfile",
        type=str,
        required=True,
        help="Output file for the rescaled tree",
    )
    parser.add_argument(
        "--precision",
        type=int,
        default=12,
        help="Decimal precision for branch lengths in output (default=12)",
    )
    return parser.parse_args()


def rescale_tree(tree, factor):
    """Rescale the branch lengths of the trees by the given factor."""
    for clade in tree.find_clades():
        if clade.branch_length is None:
            clade.branch_length = 0.0
        else:
            clade.branch_length = float(Decimal(clade.branch_length) * Decimal(factor))
    return tree


def rename_leaves(tree, labels):
    """Rename each leaf as name_{height}_{label} where height = distance from root."""
    root = tree.root
    for clade in tree.get_terminals():  # only the tips
        height = tree.distance(root, clade)
        clade.name = f"{clade.name}_{height:.6f}_{labels[clade.name]}"
    return tree


def get_tree(filename):
    with open(filename) as f:
        text = f.read()
    tree_regex = re.compile(r"<newick[^>]*>(.*?)</newick>", re.DOTALL)
    match = tree_regex.search(text)
    tree_str = match.group(1).strip() if match else text.strip()
    return Phylo.read(StringIO(tree_str), "newick")


def get_labels(filename):
    with open(filename) as f:
        text = f.read()
    leaf_regex = re.compile(r"\b(leaf_\d+)\s*=\s*(SS|NS)\b", re.IGNORECASE)
    return {name: g.upper() for name, g in leaf_regex.findall(text)}


def write_high_precision(trees, outfile, precision=12):
    """Write trees with high-precision branch lengths (requires Biopython ≥1.81)."""
    Phylo.write(trees, outfile, "newick", format_branch_length=f"%1.{precision}f")


def main():
    args = read_args()

    renamed_trees = []

    for file in args.files:
        tree = get_tree(file)
        labels = get_labels(file)
        tree = rename_leaves(tree, labels)
        renamed_trees.append(tree)

    # Rescale trees
    rescaled_trees = [rescale_tree(tree, args.sub_rate) for tree in renamed_trees]
    # Write with high precision
    write_high_precision(rescaled_trees, args.outfile, precision=args.precision)
    print(f"Rescaled trees written to {args.outfile} with {args.precision} decimals")


if __name__ == "__main__":
    main()
