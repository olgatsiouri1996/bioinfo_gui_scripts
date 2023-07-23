import os
import glob
import textwrap
from gooey import *
from pyfaidx import Fasta
# Function to split the multi-FASTA file into smaller parts
def split_multi_fasta(input_file, number, output_directory):
    # Import fasta file
    features = Fasta(input_file)
    # Get the filename without extension
    file_name = os.path.splitext(os.path.basename(input_file))[0]
    # Split list of keys (sequence identifiers)
    keyslist = list(features.keys())
    split_lists = [keyslist[x:x+number] for x in range(0, len(keyslist), number)]
    # Extract many single-fasta files
    for count, lis in enumerate(split_lists, 1):
        output_file = os.path.join(output_directory, f"{file_name}_part{count}.fasta")
        with open(output_file, 'w') as output:
            for key in lis:
                output.write(f">{features[str(key)].long_name}\n")
                seq = str(features[str(key)][:].seq)
                wrapped_seq = textwrap.fill(seq, width=60)
                output.write(wrapped_seq + '\n')

# Function to split multiple multi-FASTA files in a directory into smaller parts
def split_multi_fasta_files_in_directory(directory, number, output_directory):
    fasta_files = glob.glob(os.path.join(directory, "*.fa")) + glob.glob(os.path.join(directory, "*.fasta")) + glob.glob(os.path.join(directory, "*.fna")) + glob.glob(os.path.join(directory, "*.ffn")) + glob.glob(os.path.join(directory, "*.faa")) + glob.glob(os.path.join(directory, "*.frn"))
    for fasta_file in fasta_files:
        split_multi_fasta(fasta_file, number, output_directory)

# Command-line interface using Gooey
@Gooey(program_name='Split Multi-FASTA', header_bg_color='#DCDCDC', terminal_font_color='#DCDCDC', terminal_panel_color='#DCDCDC',default_size=(610,610))
def main():
    ap = GooeyParser(description="Split one or many multi-FASTA to many")
    ap.add_argument("program", choices=["1", "2"], default="1", help="Program to choose: 1) split 1 multi-fasta file to many, 2) split many multi-fasta files to many")
    ap.add_argument("number", type=int, help="Number of fasta records per output fasta file (you can put any number you want as it makes sure the remaining fasta records will be written to a separate file as well)")
    ap.add_argument("output_directory", type=str, widget='DirChooser', help="Output directory to save the split fasta files")

    # Arguments for Program 1
    ap.add_argument("-in", "--input", required='program' in ["1"], widget='FileChooser', help="Input multi-fasta file")

    # Arguments for Program 2
    ap.add_argument("-dir", "--directory", required='program' in ["2"], type=str, widget='DirChooser', help="Directory to search for multi-fasta files with extensions: .fasta, .fna, .ffn, .faa, .frn, .fa files")

    args = vars(ap.parse_args())

    if args['program'] == "1":
        split_multi_fasta(args['input'], args['number'], args['output_directory'])
    else:
        split_multi_fasta_files_in_directory(args['directory'], args['number'], args['output_directory'])

if __name__ == '__main__':
    main()
