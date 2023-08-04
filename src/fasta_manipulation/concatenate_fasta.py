# python3
import glob
import os
from gooey import Gooey, GooeyParser
# input parameters
@Gooey(required_cols=2, program_name="", header_bg_color= "#DCDCDC", terminal_font_color= "#DCDCDC", terminal_panel_color= "#DCDCDC",default_size=(610,610))
def main():
    # Create the GooeyParser
    parser = GooeyParser(description="Concatenate multiple FASTA files")
    # Add arguments for directory path, output file, and line width
    parser.add_argument("input directory", widget="DirChooser", help="Choose the directory containing the files you want to concatenate")
    parser.add_argument("output file", widget="FileSaver", help="Choose the output file path")
    parser.add_argument("pattern type", default="search for all fasta extensions", choices=["search for all fasta extensions","prefix","intermediate","suffix"], help="Import pattern type of the inport files you want to concatenate")
    parser.add_argument("-pat","--pattern", required=False, help="Import pattern of the inport files you want to concatenate")
    # Parse the arguments
    args = vars(parser.parse_args())
    # Find all the FASTA files in the specified directory
    # Choose pattern type
    if args["pattern type"] == "search for all fasta extensions":
        fasta_files = glob.glob(os.path.join(args['input directory'], "*.fa")) + glob.glob(os.path.join(args['input directory'], "*.fasta")) + glob.glob(os.path.join(args['input directory'], "*.fna")) + glob.glob(os.path.join(args['input directory'], "*.ffn")) + glob.glob(os.path.join(args['input directory'], "*.faa")) + glob.glob(os.path.join(args['input directory'], "*.frn"))
    elif args["pattern type"] == "prefix":
        fasta_files = glob.glob(os.path.join(args['input directory'], args["pattern"]+"*"))
    elif args["pattern type"] == "intermediate":
        fasta_files = glob.glob(os.path.join(args['input directory'], "*"+args["pattern"]+"*"))
    else:
        fasta_files = glob.glob(os.path.join(args['input directory'], "*"+args["pattern"]))
    # Open the output file in write mode
    with open(args['output file'], 'w') as outfile:
        for fasta_file in fasta_files:
            # Read the content of each FASTA file
            with open(fasta_file, 'r') as infile:
                outfile.write(infile.read())

if __name__ == "__main__":
    main()
