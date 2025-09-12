import argparse
from Bio import SeqIO

def extract_sequences(input_fasta, output_fasta, names_list):
    with open(input_fasta, "r") as input_handle, open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if any(name in record.id for name in names_list):
                SeqIO.write(record, output_handle, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract sequences from a FASTA file.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("output_fasta", help="Output FASTA file")
    parser.add_argument("names_file", help="File containing names to extract")
    
    args = parser.parse_args()

    with open(args.names_file, "r") as f:
        names_list = [line.strip() for line in f]

    extract_sequences(args.input_fasta, args.output_fasta, names_list)
