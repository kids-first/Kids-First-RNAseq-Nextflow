#!/usr/bin/env python3
"""
Merge data from three CSV files based on genomic coordinates.
Each input CSV file contains columns: "coordinates", "circn", "strand", "gene", and a sample-specific score column.
The output CSV will contain merged data with scores from all three files for matching coordinates.

Usage: python script.py cx.csv dcc.csv fc.csv output.csv
"""

import csv
import re
import sys

def chr_key(chrom: str) -> int | float:
    """
    Convert chromosome string to a sortable key.

    Arguments:
        chrom: Chromosome string (e.g., "chr1", "chrX", "chrY", "chrM")
    Returns:
        Integer or float representing the chromosome number for sorting purposes.

    """
    match: re.Match = re.match(r'chr(\d+|X|Y|M)', chrom)
    if match:
        val = match.group(1)
        if val.isdigit():
            return int(val)
        elif val == 'X':
            return 23
        elif val == 'Y':
            return 24
        elif val == 'M':
            return 25
    return float('inf')

def coord_key(s: str) -> tuple:
    """
    Convert coordinate string to a sortable key.
    Arguments:
        s: Coordinate string in the format "chrN:start-end"
    Returns:
        Tuple containing chromosome key and start position as integer.
    """
    chrom, pos = s.split(':')
    start, end = pos.split('-')
    return (chr_key(chrom), int(start))

def read_csv(file_path: str) -> dict[str, dict[str, str | float]]:
    """
    Read CSV file and return a dictionary with coordinates as keys.
    Arguments:
        file_path: Path to the CSV file.
    Returns:
        Dictionary with coordinates as keys and a dictionary of relevant data (circ, strand, gene, sample_name, score) as values.
    """
    data = {}
    with open(file_path, mode='r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        sample_name = reader.fieldnames[-1]
        for row in reader:
            key = row['coordinates']
            data[key] = {k: row.get(k) for k in ["circn", "strand", "gene"]} | {"sample_name": sample_name} | {"score": row[sample_name]}
    return data

def write_csv(file_path: str, data: list[dict], fieldnames: list[str]) -> None:
    """
    Write data to a CSV file.
    Arguments:
        file_path: Path to the output CSV file.
        data: List of dictionaries containing the data to write.
        fieldnames: List of field names for the CSV header.
    Returns:
        None
    """
    with open(file_path, mode='w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

def main():
    if len(sys.argv) != 5:
        print(f"Usage: python {sys.argv[0]} cx.csv dcc.csv fc.csv output.csv")
        sys.exit(1)

    cx_file: str = sys.argv[1]
    dcc_file: str = sys.argv[2]
    fc_file: str = sys.argv[3]
    output_file: str = sys.argv[4]

    cx_data: dict = read_csv(cx_file)
    dcc_data: dict = read_csv(dcc_file)
    fc_data: dict = read_csv(fc_file)

    merged_data = []
    sorted_keys: list[str] = sorted(cx_data.keys(), key=coord_key)
    for key in sorted_keys:
        if key in dcc_data and key in fc_data:
            merged_data.append({
                "sample_name": cx_data[key]["sample_name"],
                "coordinates": key,
                "strand": cx_data[key]["strand"],
                "gene": cx_data[key]["gene"],
                "circn": cx_data[key]["circn"],
                "cx_score": cx_data[key]["score"],
                "dcc_score": dcc_data[key]["score"],
                "fc_score": fc_data[key]["score"]
            })
        else:
            print(f"Warning: Missing data for coordinates {key=} in one of the files.", file=sys.stderr)

    fieldnames: list[str] = list(merged_data[0].keys())
    write_csv(output_file, merged_data, fieldnames)

if __name__ == "__main__":
    main()

