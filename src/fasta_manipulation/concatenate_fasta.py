# python3
import glob
import os
from gooey import Gooey, GooeyParser
# input parameters
@Gooey(required_cols=2, program_name="", header_bg_color= "#DCDCDC", terminal_font_color= "#DCDCDC", terminal_panel_color= "#DCDCDC", default_size=(800,330))
def main():
    # Create the GooeyParser
    parser = GooeyParser(description="Concatenate multiple FASTA files")
    # Add arguments for directory path, output file, and line width
    parser.add_argument("-dir","--input directory", widget="DirChooser", required=True, help="Choose the directory containing: .fasta, .fna, .ffn, .faa, .frn, .fa files")
    parser.add_argument("-out","--output file", widget="FileSaver", required=True, help="Choose the output file path")
    # Parse the arguments
    args = vars(parser.parse_args())
    # Find all the FASTA files in the specified directory
    fasta_files = glob.glob(os.path.join(args['input directory'], '*.f*'))
    # Open the output file in write mode
    with open(args['output file'], 'w') as outfile:
        for fasta_file in fasta_files:
            # Read the content of each FASTA file
            with open(fasta_file, 'r') as infile:
                outfile.write(infile.read())

if __name__ == "__main__":
    main()
