#!/usr/bin/env python3
"""
DNA Variant Generator and Protein Sequence Translator
----------------------------------------------------

Author: Khan Lodoi <lodoi.khanburged@jichi.ac.jp>
Date: 10/29/2025
Version: 1.0.0

Description:
    This script generates DNA sequence variants by applying specified substitutions (A->G, T->C, or both) 
    within a given range of indices. It then translates these DNA variants into their corresponding protein sequences 
    based on the standard codon table. The results are saved in an output CSV file with both the variants and their 
    associated protein sequences.

Usage:
    python dna_to_protein_variants.py --input <input_csv_file> --output <output_csv_file> --substitution <substitution_type>

    Arguments:
        --input <input_csv_file>      : Path to the input CSV file containing DNA sequences and range columns.
                                        The CSV must contain at least two columns:
                                            - dna_sequence: The original DNA sequence.
                                            - range: The range (e.g., "5-10") indicating which positions to substitute.
        --output <output_csv_file>    : Path to the output CSV file where the variant DNA sequences and their protein 
                                        sequences will be written.
        --substitution <substitution_type> : The type of substitution to apply:
                                                - "AtoG"   : Substitute A with G within the specified range.
                                                - "TtoC"   : Substitute T with C within the specified range.
                                                - "both"   : Apply both A->G and T->C substitutions. Default is "both".

    Example:
        python dna_to_protein_variants.py --input input_file.csv --output output_file.csv --substitution AtoG

    Input CSV Format:
        The input CSV file must contain two columns: 'dna_sequence' and 'range'. Example:
        
        | dna_sequence       | range   |
        |--------------------|---------|
        | ATCGATTGAGCTGATCG  | 5-10    |
        | GCTAGCTAGCTAGCTA   | 2-8     |

    Output CSV Format:
        The output CSV file will contain the original DNA sequence, variant sequences (after substitution), 
        and their corresponding protein sequences. Example:
        
        | dna_sequence       | range   | variant_dna       | variant_protein_sequence |
        |--------------------|---------|-------------------|--------------------------|
        | ATCGATTGAGCTGATCG  | 5-10    | ATCGAGTGAGCTGATCG | M*...                    |
        |                    |         | ATCGGGTGAGCTGATCG | M*...                    |

    Requirements:
        No external libraries are required, this script uses Python's standard libraries only.

    Error Handling:
        - If the input CSV is missing required columns ('dna_sequence' or 'range'), the script will exit with an error message.
        - If a sequence range is invalid or in an incorrect format, a warning will be issued and the row will be skipped.
        - If an invalid substitution type is provided, the script will terminate with a usage message.

    Note:
        - The script assumes the input CSV is well-formed and contains the required columns.
        - The input sequences should be uppercase, and the range should be in the format "start-end".
        - Invalid ranges or sequences may cause errors, which will be printed to stderr.

"""
import csv
import argparse
import sys
import itertools

# Standard codon table
CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
}

def dna_to_protein(dna_seq: str) -> str:
    """Translate a DNA sequence into a protein sequence."""
    dna_seq = dna_seq.replace(' ', '').upper()
    protein_seq = ''
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        codon_clean = codon.replace('-', 'A')
        protein_seq += CODON_TABLE.get(codon_clean, 'X')
    return protein_seq

def generate_variants(dna_seq: str, start: int, end: int, substitution: str):
    """
    Generate all combinations based on the selected substitution(s):
    - A -> G (if 'AtoG' selected)
    - T -> C (if 'TtoC' selected)
    """
    seq_list = list(dna_seq)
    
    # Choose substitutions based on the flag
    if substitution == 'AtoG':
        # Find indices where A is located within the specified range
        indices = [i for i in range(start, end) if seq_list[i].upper() == 'A']
        replacements = ['A', 'G']
    elif substitution == 'TtoC':
        # Find indices where T is located within the specified range
        indices = [i for i in range(start, end) if seq_list[i].upper() == 'T']
        replacements = ['T', 'C']
    else:
        # For 'both', handle A -> G and T -> C
        indices = [i for i in range(start, end) if seq_list[i].upper() in ['A', 'T']]
        replacements = ['A', 'G', 'T', 'C']

    # No changes in range, return the original sequence
    if not indices:
        return [dna_seq]

    variants = []
    # Generate all combinations of replacements
    for combo in itertools.product(replacements, repeat=len(indices)):
        seq_copy = seq_list[:]
        for idx, nucleotide in zip(indices, combo):
            seq_copy[idx] = nucleotide
        variants.append(''.join(seq_copy))
    return variants

def process_csv(input_file: str, output_file: str, substitution: str):
    """Generate specified variants (A->G or T->C or both) and protein sequences for each DNA sequence."""
    with open(input_file, newline='') as csvfile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(csvfile)
        fieldnames = reader.fieldnames + ['variant_dna', 'variant_protein_sequence']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        if 'dna_sequence' not in reader.fieldnames or 'range' not in reader.fieldnames:
            print("Error: CSV must have 'dna_sequence' and 'range' columns", file=sys.stderr)
            sys.exit(1)

        for row in reader:
            dna_seq = row['dna_sequence']
            range_str = row['range']
            
            # Parse the start and end from the 'range' column
            try:
                start, end = map(int, range_str.split('-'))
            except ValueError:
                print(f"Warning: Invalid range format '{range_str}' for DNA sequence {dna_seq}. Skipping.", file=sys.stderr)
                continue
            
            # Generate all variants based on the selected substitution
            variants = generate_variants(dna_seq, start, end, substitution)
            
            # First variant (in the same row as the original sequence)
            first_variant = variants[0]
            first_variant_protein = dna_to_protein(first_variant)
            
            # Write the original sequence and the first variant in the same row
            row_copy = row.copy()
            row_copy['variant_dna'] = first_variant  # First variant dna_sequence
            row_copy['variant_protein_sequence'] = first_variant_protein  # Protein for the first variant
            # Also include the original sequence in the same row
            row_copy['dna_sequence'] = dna_seq  # Original sequence in the same row
            writer.writerow(row_copy)

            # Subsequent variants (new rows with empty dna_sequence)
            for var in variants[1:]:
                row_copy = row.copy()
                row_copy['dna_sequence'] = ''  # Empty dna_sequence for variants
                row_copy['variant_dna'] = var  # Only variant dna_sequence here
                row_copy['variant_protein_sequence'] = dna_to_protein(var)  # Protein for the variant
                writer.writerow(row_copy)

            print(f"{len(variants)} variants generated for DNA: {dna_seq}")

def main():
    parser = argparse.ArgumentParser(
        description="Generate DNA variants with specified substitutions (A->G, T->C) and compute protein sequences."
    )
    parser.add_argument("--input", required=True, help="Input CSV with dna_sequence and range columns")
    parser.add_argument("--output", required=True, help="Output CSV with variants and variant protein sequences")
    parser.add_argument("--substitution", choices=['AtoG', 'TtoC', 'both'], default='both', help="Specify which substitution(s) to apply: 'AtoG', 'TtoC', or 'both' (default)")

    args = parser.parse_args()

    process_csv(args.input, args.output, args.substitution)

if __name__ == "__main__":
    main()
