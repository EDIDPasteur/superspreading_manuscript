#!/usr/bin/env python
# Add dates and type ? to sequence names

from Bio import SeqIO
from argparse import ArgumentParser
import pickle
from datetime import datetime

# COnstants
DATE_FORMAT = "%m/%d/%y"

def parse_args():
    parser = ArgumentParser(description="Rename sequences in Fasta file to append the date at the end")
    parser.add_argument("-f","--fasta_file", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-o","--output_file", required=True, help="Path to the output FASTA file.")
    parser.add_argument("-d", "--dates_object", required=True, help="Path to the dates object file.")
    return parser.parse_args()

def load_dates(dates_file):
    with open(dates_file, 'rb') as f:
        dates_dict = pickle.load(f)
    return dates_dict

def date_to_decimal(date_str, date_format=DATE_FORMAT):
    """Convert a date string in the format %m/%d/%y to a decimal date"""
    dt = datetime.strptime(date_str, date_format)
    year = dt.year

    # Year day
    year_day = dt.timetuple().tm_yday

    # Days in year
    days_in_year = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365

    # Convert to decimal
    decimal_date = year + (year_day - 1) / days_in_year
    return decimal_date

def rename_sequences(fasta_file, output_file, dates_dict):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    renamed_records = []
    
    for record in records:
        seq_id = record.id
        if seq_id in dates_dict:
            date = date_to_decimal(dates_dict[seq_id])
            new_id = f"{seq_id}_N_{date}"
            record.id = new_id
            record.description = ""
            renamed_records.append(record)
        else:
            print(f"Warning: No date found for sequence ID {seq_id}. Skipping.")
    
    SeqIO.write(renamed_records, output_file, "fasta")


def main():
    args = parse_args()
    dates_dict = load_dates(args.dates_object)
    rename_sequences(args.fasta_file, args.output_file, dates_dict)
if __name__ == "__main__":
    main()
    