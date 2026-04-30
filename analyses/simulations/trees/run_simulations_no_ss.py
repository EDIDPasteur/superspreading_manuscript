#!/bin/python
# This script is used to update the ReMaster master script and BDMMprime master script
# Jordi Sevilla, 31/07/2025

from argparse import ArgumentParser
import os
import dendropy
import time
from concurrent.futures import ThreadPoolExecutor


def read_arguments():
    parser = ArgumentParser(description="Update ReMaster and BDMMprime master scripts.")
    parser.add_argument(
        "-rf",
        "--re_master_file",
        type=str,
        required=True,
        help="Path to the ReMaster master script file.",
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
    return parser.parse_args()


def update_re_master_script(file_path, new_filepath):
    with open(file_path, "r") as file:
        content = file.read()

    with open(new_filepath, "w") as file:
        file.write(content)


def update_bdmmprime_script(file_path, new_filepath, tree):
    common_alignment_line = (
        '<sequence id="dir_{}" spec="Sequence" taxon="{}" totalcount="4" value="nnn"/>'
    )
    newick = tree._as_newick_string()
    alignment = ""
    tiptypes = ""
    dates = ""
    for node in tree.leaf_nodes():
        taxon = node.taxon.label.replace(" ", "_")
        alignment += common_alignment_line.format(taxon, taxon) + "\n"
        tiptypes += f"{taxon}={node.annotations['type'].value},\n"
        dates += f"{taxon}={node.annotations['time'].value},"
    tiptypes = tiptypes.rstrip(",\n")  # Remove trailing comma and newline
    dates = dates.rstrip(",")  # Remove trailing comma
    with open(file_path, "r") as file:
        content = file.read()
        content = content.replace("[ALIGNMENT]", alignment)
        content = content.replace("[TREE]", newick)
        content = content.replace("[PHENOTYPES]", tiptypes)
        content = content.replace("[DATES]", dates)

    with open(new_filepath, "w") as file:
        file.write(content)


def update_bdmmprime_script_no_tiptypes(file_path, new_filepath, tree):
    common_alignment_line = (
        '<sequence id="dir_{}" spec="Sequence" taxon="{}" totalcount="4" value="nnn"/>'
    )
    newick = tree._as_newick_string()
    alignment = ""
    tiptypes = ""
    dates = ""
    for node in tree.leaf_nodes():
        taxon = node.taxon.label.replace(" ", "_")
        alignment += common_alignment_line.format(taxon, taxon) + "\n"
        tiptypes += f"{taxon}=?,\n"
        dates += f"{taxon}={node.annotations['time'].value},"
    tiptypes = tiptypes.rstrip(",\n")  # Remove trailing comma and newline
    dates = dates.rstrip(",")  # Remove trailing comma
    with open(file_path, "r") as file:
        content = file.read()
        content = content.replace("[ALIGNMENT]", alignment)
        content = content.replace("[TREE]", newick)
        content = content.replace("[PHENOTYPES]", tiptypes)
        content = content.replace("[DATES]", dates)

    with open(new_filepath, "w") as file:
        file.write(content)


def process_combination_of_parameters(args, commands):
    # Update ReMaster master script
    re_master_file = args.re_master_file
    outdir = args.dir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    name = f"ReMaster_0_1"
    newoutdir = os.path.join(outdir, name)
    if not os.path.exists(newoutdir):
        os.makedirs(newoutdir)

    # Update the ReMaster master script with the new parameters
    new_re_master_file = os.path.join(newoutdir, f"{name}.xml")
    update_re_master_script(re_master_file, new_re_master_file)

    # Correct trees
    correct_tree = False
    while not correct_tree:
        # Run BEAST on the ReMaster master script
        beast_command = f"{args.beast} -overwrite -working {new_re_master_file}"
        print(f"Running BEAST command: {beast_command}")
        os.system(beast_command)

        # Update BDMMprime master script
        bdmmprime_file = args.bdmmprime_file
        bdmmprime_parent_dir = os.path.join(newoutdir, "BDMMprime")
        if not os.path.exists(bdmmprime_parent_dir):
            os.makedirs(bdmmprime_parent_dir)

        # Load the tree from the ReMaster master script
        # First modify the tree file
        file = f"{os.path.join(newoutdir, f'{name}_untyped.trees')}"
        with open(file, "r") as fi:
            content = fi.read()
            content = content.replace("500 leaf_499,", "500 leaf_499")
        with open(file, "w") as fi:
            fi.write(content)

        try:
            # Load the tree using dendropy
            trees = list(
                dendropy.TreeList.get(
                    path=f"{os.path.join(newoutdir, f'{name}_untyped.trees')}",
                    schema="nexus",
                )
            )
            correct_tree = True
        except:
            print("Retrying to run BEAST on the ReMaster master script...")
            time.sleep(5)

    for i, tree in enumerate(trees):
        bdmmprime_dir = os.path.join(bdmmprime_parent_dir, f"tree_{i}")
        if not os.path.exists(bdmmprime_dir):
            os.makedirs(bdmmprime_dir)

        if not os.path.exists(
            os.path.join(bdmmprime_dir, f"BDMMprime_0_1_tree{i}_no_tiptypes.log")
        ):
            # Without tip types
            new_bdmmprime_file_no_tiptypes = os.path.join(
                bdmmprime_dir, f"BDMMprime_0_1_tree{i}_no_tiptypes.xml"
            )
            update_bdmmprime_script_no_tiptypes(
                bdmmprime_file, new_bdmmprime_file_no_tiptypes, tree
            )


def main():
    args = read_arguments()

    processes = []
    commands = []
    # Iterate over the combinations of parameters to get all the commands
    process_combination_of_parameters(args, commands)

    print(f"XML files have been created")


if __name__ == "__main__":
    main()
