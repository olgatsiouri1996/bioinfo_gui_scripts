import os
import glob
from gooey import *
from pyfaidx import Fasta
import textwrap

# Input parameters
@Gooey(required_cols=0, program_name='fasta formatter', header_bg_color='#DCDCDC', terminal_font_color='#DCDCDC', terminal_panel_color='#DCDCDC', default_size=(610, 660))
def main():
    ap = GooeyParser(description="changes the number of characters per FASTA sequence line in 1 or many FASTA files")
    ap.add_argument("program", type=str, default="one input/output fasta file", choices=["one input/output fasta file", "many input/output fasta files"], help="program to choose")
    ap.add_argument("width", type=int, default=60, choices=[60, 70, 80, 100, 120], help="number of characters per line")
    ap.add_argument("-in", "--input", required=False, widget="FileChooser", help="input fasta file")
    ap.add_argument("-out", "--output", required=False, widget="FileSaver", help="output fasta file")
    ap.add_argument("-indir", "--input directory", required=False, type=str, widget='DirChooser', help="directory to search for fasta files files with extensions: .fasta, .fna, .ffn, .faa, .frn, .fa")
    ap.add_argument("-outdir", "--output directory", required=False, type=str, widget='DirChooser', help="directory to save fasta files")
    args = vars(ap.parse_args())

    # Choose program
    if args['program'] == "one input/output fasta file":
        # Export to a new fasta file
        with open(args['output'], 'w') as output_file:
            fasta = Fasta(args['input'])
            for record in fasta:
                formatted_sequence = textwrap.fill(str(record), width=args['width'])
                output_file.write(f">{record.long_name}\n")
                output_file.write(formatted_sequence + "\n")

    else:
        # Import each fasta file from the working directory
        for filename in glob.glob(os.path.join(args['input directory'], "*.fa")) + glob.glob(os.path.join(args['input directory'], "*.fasta")) + glob.glob(os.path.join(args['input directory'], "*.fna")) + glob.glob(os.path.join(args['input directory'], "*.ffn")) + glob.glob(os.path.join(args['input directory'], "*.faa")) + glob.glob(os.path.join(args['input directory'], "*.frn")):
            output_filename = f"{os.path.splitext(os.path.basename(filename))[0]}_w{args['width']}.fasta"
            with open(os.path.join(args['output directory'], output_filename), 'w') as output_file:
                fasta = Fasta(filename)
                for record in fasta:
                    formatted_sequence = textwrap.fill(str(record), width=args['width'])
                    output_file.write(f">{record.long_name}\n")
                    output_file.write(formatted_sequence + "\n")

if __name__ == '__main__':
    main()
