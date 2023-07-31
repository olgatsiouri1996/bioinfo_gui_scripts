import os
import glob
from gooey import *
from Bio import SeqIO
import pandas as pd
import textwrap

# input parameters
@Gooey(required_cols=0, program_name='fasta formatter', header_bg_color='#DCDCDC', terminal_font_color='#DCDCDC', terminal_panel_color='#DCDCDC',default_size=(610,645))
def main():
    ap = GooeyParser(description="changes the width of sequences line in 1 or many FASTA files")
    ap.add_argument("program", type=str, default="one input/output fasta file", choices=["one input/output fasta file","many input/output fasta files",".txt file with fasta file names and width for each file"], help="program to choose")
    ap.add_argument("-in", "--input", required=False, widget="FileChooser", help="input fasta file")
    ap.add_argument("-txt", "--txt", required=False, widget="FileChooser", help="input txt file with 2 columns: filename without extension and width")
    ap.add_argument("-out", "--output", required=False, widget="FileSaver", help="output fasta file")
    ap.add_argument("-indir", "--input directory", required=False, type=str, widget='DirChooser', help="directory to search for fasta files")
    ap.add_argument("-outdir", "--output directory", required=False, type=str, widget='DirChooser', help="directory to save fasta files")
    ap.add_argument("-width", "--width", required=False, type=int, default=80, choices=[60,70,80,100,120], help="number of characters per line")
    args = vars(ap.parse_args())
    # create function to retrieve fasta description
    def retrieve_description(desc):
        try:
            description = str(desc).split(" ",1)[1]
        except:
            description = ""
        return description
    # choose program
    if args['program'] == "one input/output fasta file":
        # export to a new fasta file
        with open(args['output'], 'w') as output_file:
            for record in SeqIO.parse(args['input'], 'fasta'):
                formatted_sequence = textwrap.fill(str(record.seq), width=args['width'])
                output_file.write(f">{record.id} {retrieve_description(record.description)}\n")
                output_file.write(formatted_sequence + "\n")
                
    elif args['program'] == "many input/output fasta files":
        # import each fasta file from the working directory
        for filename in glob.glob(os.path.join(args['directory'], "*.fasta")):
            output_filename = f"{os.path.splitext(os.path.basename(filename))[0]}_w{args['width']}.fasta"
            with open(os.path.join(args['output directory'],output_filename), 'w') as output_file:
                for record in SeqIO.parse(filename, 'fasta'):
                    formatted_sequence = textwrap.fill(str(record.seq), width=args['width'])
                    output_file.write(f">{record.id} {retrieve_description(record.description)}\n")
                    output_file.write(formatted_sequence + "\n")
    
    else:
        df = pd.read_csv(args['txt'], header=None, sep="\t")
        headers = df.iloc[:, 0].values.tolist()
        widths = df.iloc[:, 1].values.tolist()
        for (a, b) in zip(headers, widths):
            input_filename = f"{a}.fasta"
            output_filename = f"{a}_w{b}.fasta"
            with open(os.path.join(args['output directory'],output_filename), 'w') as output_file:
                for record in SeqIO.parse(input_filename, 'fasta'):
                    formatted_sequence = textwrap.fill(str(record.seq), width=int(b))
                    output_file.write(f">{record.id} {retrieve_description(record.description)}\n")
                    output_file.write(formatted_sequence + "\n")


if __name__ == '__main__':
    main()
