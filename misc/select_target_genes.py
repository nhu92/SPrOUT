import argparse
from Bio import SeqIO

def read_numbers_from_file(number_file):
    # Read the list of numbers from the provided file
    with open(number_file, 'r') as f:
        numbers = [line.strip() for line in f]
    return numbers

def filter_fasta_by_numbers(fasta_file, numbers, output_file):
    # Open the output file to write the selected sequences
    with open(output_file, 'w') as out_f:
        # Parse the fasta file
        for record in SeqIO.parse(fasta_file, 'fasta'):
            # Check if any number from the list is in the sequence name
            for number in numbers:
                if f"-{number}" in record.id:
                    # If a match is found, write the sequence to the output file
                    SeqIO.write(record, out_f, 'fasta')

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Filter sequences in a FASTA file by matching numbers in sequence names.")
    
    # Input files
    parser.add_argument("number_file", type=str, help="Path to the file containing list of numbers")
    parser.add_argument("fasta_file", type=str, help="Path to the input FASTA file")
    
    # Output file
    parser.add_argument("output_file", type=str, help="Path to the output FASTA file")

    # Parse arguments
    args = parser.parse_args()
    
    # Read numbers from the file
    numbers = read_numbers_from_file(args.number_file)
    
    # Filter the fasta sequences by the numbers
    filter_fasta_by_numbers(args.fasta_file, numbers, args.output_file)

    print(f"Filtered sequences saved to {args.output_file}")

if __name__ == "__main__":
    main()
