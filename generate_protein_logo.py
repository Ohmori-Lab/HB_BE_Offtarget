#!/usr/bin/env python3
"""
Protein Sequence Logo Generator
-----------------------------------------------------

Author: Khan Lodoi <lodoi.khanburged@jichi.ac.jp>
Date: 10/29/2025
Version: 1.0.0

Description:
    This script generates a protein sequence logo from a CSV file containing protein sequences.
    The sequences are processed to produce a Position Specific Scoring Matrix (PSSM), which is then
    used to create a sequence logo, which is saved as a PNG image.

    A sequence logo is a graphical representation of sequence conservation at each position in a set of aligned sequences.
    The heights of the letters in the logo represent the relative frequency of each amino acid at that position.

Usage:
    python generate_protein_logo.py --input <input_csv_file> --output <output_png_file> --column <column_name>

    Arguments:
        --input <input_csv_file>      : Path to the input CSV file containing protein sequences.
                                        The file should contain a column with protein sequences (default 'variant_protein_sequence').
        --output <output_png_file>    : Path to the output PNG file where the sequence logo will be saved.
                                        The output file must have a '.png' extension.
        --column <column_name>        : The name of the column in the CSV that contains the protein sequences (default is 'variant_protein_sequence').

    Example:
        python generate_protein_logo.py --input sequences.csv --output logo.png --column variant_protein_sequence

    Input CSV Format:
        The input CSV file must have a column containing protein sequences. Example:

        | variant_protein_sequence |
        |--------------------------|
        | MKTIIALSYIFCLVFA         |
        | MKTIIALSYIFCLVFA         |
        | MKTIIALSYIFCLVFA         |

    Output:
        The script generates a sequence logo in PNG format, visualizing the conservation of amino acids across the sequences at each position.

    Requirements:
        - pandas
        - numpy
        - matplotlib
        - logomaker

    Installation (if not already installed):
        pip install pandas numpy matplotlib logomaker

    Error Handling:
        - The script will raise an error if the input file does not exist, if the required column is missing, 
          or if invalid characters are found in the sequences.
        - The output file must have a '.png' extension.
        - Sequences should only contain valid amino acids (ACDEFGHIKLMNPQRSTVWY) and optionally '*' for gaps.
        - Sequences must be of the same length.

"""
import os
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import argparse
from typing import List

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Constants for Amino Acids and mapping
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
AA_DICT = {aa: idx for idx, aa in enumerate(AMINO_ACIDS)}

# Error Messages
ERROR_FILE_NOT_FOUND = "The file at {} does not exist."
ERROR_MISSING_COLUMN = "CSV file must contain a '{}' column."
ERROR_NO_SEQUENCES = "No sequences found in the CSV file."
ERROR_INVALID_OUTPUT_FORMAT = "Output file must have a '.png' extension."
ERROR_INVALID_CHARACTER = "Invalid character found in sequence: {}"


def read_sequences_from_csv(file_path: str, column_name: str) -> List[str]:
    """
    Reads protein sequences from a CSV file and returns them as a list.

    Args:
        file_path (str): Path to the input CSV file.
        column_name (str): The name of the column containing protein sequences.

    Returns:
        List[str]: List of protein sequences.
    
    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the required column is missing or if no sequences are found.
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(ERROR_FILE_NOT_FOUND.format(file_path))
    
    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        logger.error(f"Error reading CSV file: {e}")
        raise

    if column_name not in df.columns:
        raise ValueError(ERROR_MISSING_COLUMN.format(column_name))
    
    sequences = df[column_name].dropna().tolist()
    
    if not sequences:
        raise ValueError(ERROR_NO_SEQUENCES)
    
    logger.info(f"Successfully read {len(sequences)} sequences from {file_path}")
    return sequences


def validate_sequences(sequences: List[str]) -> None:
    """
    Validates that all sequences consist only of valid amino acids.

    Args:
        sequences (List[str]): List of protein sequences.
    
    Raises:
        ValueError: If any sequence contains invalid characters.
    """
    for seq in sequences:
        for aa in seq:
            if aa not in AMINO_ACIDS and aa != '*' and aa != '-':
                logger.error(f"Invalid character '{aa}' found in sequence: {seq}")
                raise ValueError(ERROR_INVALID_CHARACTER.format(aa))


def sequences_to_pssm(sequences: List[str]) -> pd.DataFrame:
    """
    Converts a list of protein sequences into a normalized Position Specific Scoring Matrix (PSSM).

    Args:
        sequences (List[str]): List of protein sequences.

    Returns:
        pd.DataFrame: Normalized PSSM matrix.
    
    Raises:
        ValueError: If sequences are of different lengths.
    """
    sequence_length = len(sequences[0])

    # Ensure all sequences are the same length
    if any(len(seq) != sequence_length for seq in sequences):
        raise ValueError("All sequences must be of the same length.")

    pssm = np.zeros((len(sequences), sequence_length, len(AMINO_ACIDS)))
    all_stars = np.ones(sequence_length)  # Track positions with '*'

    for i, seq in enumerate(sequences):
        for j, aa in enumerate(seq):
            if aa == '*':
                continue  # Skip '*' in the PSSM matrix
            if aa in AA_DICT:
                pssm[i, j, AA_DICT[aa]] = 1
                all_stars[j] = 0  # Mark as not all '*' for this position

    # Handle positions with all '*' as NaN
    pssm[:, all_stars == 1, :] = np.nan
    
    pssm_sum = np.sum(pssm, axis=0) + 1e-10  # Prevent division by zero
    pssm_norm = pssm_sum / np.sum(pssm_sum, axis=1, keepdims=True)
    pssm_norm = np.nan_to_num(pssm_norm)  # Replace NaNs with zeros

    pssm_df = pd.DataFrame(pssm_norm, columns=list(AMINO_ACIDS))
    
    logger.info("PSSM matrix generated successfully.")
    return pssm_df


def generate_sequence_logo(sequences: List[str], output_file: str) -> None:
    """
    Generates and saves a sequence logo based on the given protein sequences.

    Args:
        sequences (List[str]): List of protein sequences.
        output_file (str): Path to save the generated logo image.
    
    Raises:
        ValueError: If the output file is not a PNG file.
    """
    if not output_file.endswith('.png'):
        raise ValueError(ERROR_INVALID_OUTPUT_FORMAT)
    
    # Generate the PSSM dataframe
    pssm_df = sequences_to_pssm(sequences)
    
    # Adjust index to start from 1 (instead of the default 0)
    pssm_df.index = pssm_df.index + 1  # Shift index to start at 1
    
    # Generate the logo
    logo = logomaker.Logo(pssm_df)
    
    # Save the logo to file
    logo.fig.savefig(output_file, format='png')
    logger.info(f"Logo saved as {output_file}")


def parse_arguments() -> argparse.Namespace:
    """
    Parses command-line arguments and validates them.
    
    Returns:
        argparse.Namespace: Parsed arguments.
    
    Raises:
        ValueError: If the required arguments are missing or invalid.
    """
    parser = argparse.ArgumentParser(description="Generate a protein sequence logo from a CSV file.")
    
    parser.add_argument('--input', type=str, required=True, help="Path to the input CSV file containing protein sequences.")
    parser.add_argument('--output', type=str, required=True, help="Path to the output PNG file where the logo will be saved.")
    parser.add_argument('--column', type=str, default='variant_protein_sequence', help="Name of the column containing protein sequences.")
    
    args = parser.parse_args()

    # Validate input file
    if not os.path.isfile(args.input):
        raise ValueError(ERROR_FILE_NOT_FOUND.format(args.input))
    
    return args


def main() -> None:
    """
    Main function that orchestrates reading sequences, generating PSSM, and saving the sequence logo.
    """
    try:
        # Parse arguments
        args = parse_arguments()
        
        # Read protein sequences
        sequences = read_sequences_from_csv(args.input, args.column)
        
        # Validate sequences
        validate_sequences(sequences)
        
        # Generate and save sequence logo
        generate_sequence_logo(sequences, args.output)

    except Exception as e:
        logger.error(f"Error occurred: {e}")
        raise


if __name__ == "__main__":
    main()
