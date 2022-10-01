# python3
import os
from gooey import *
import warnings
from Bio import BiopythonWarning
from Bio import SeqIO
# input parameters
@Gooey(required_cols=0, program_name='subset pdb to fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="converts a pdb file into a fasta file based on the chain, start and end locations")
    ap.add_argument("-in", "--input", required=False, widget='FileChooser', help="input pdb file(all names allowed except from \"out.fasta\") with HEADER section e.g the 1tii structure from Protein Data Bank")
    ap.add_argument("-dir", "--directory", required=False, type=str, widget='DirChooser', help="directory to search for pdb files")
    ap.add_argument("-chain", "--chain", required=False, default='A', help="chain from pdb file to select. Default is A")
    ap.add_argument("-start", "--start", required=False, default=1, type=int, help="amino acid in chain to start writing the fasta file")
    ap.add_argument("-end", "--end", required=False, type=int, help="amino acid in chain to end writing the fasta file")
    ap.add_argument("-type", "--type", required=False,default=1, type=int, widget='Dropdown', choices=[1,2],  help="type of input to choose: 1) 1 pdb file, 2) many pdb files.")
    ap.add_argument("-pro", "--program", required=False,default=1, type=int, widget='Dropdown', choices=[1,2], help="program to choose: 1) add both start and end location 2) the end location will be that of the latest amino acid in the chain.")
    args = vars(ap.parse_args())
    # main
    # ignore warnings
    warnings.simplefilter('ignore', BiopythonWarning)
    # create function to trim + convert to fasta
    def trim_to_fasta(fi):
        count = SeqIO.convert(fi, "pdb-atom", "out.fasta", "fasta")
        # insert fasta file
        for record in SeqIO.parse("out.fasta", "fasta"):
            if ''.join([":",args['chain']]) in record.id: # select sequence by chain
            # choose program
                if args['program'] == 1:
                    aa_end = args['end']
                else:
                    aa_end = len(record.seq)
            # trim input fasta sequence
                trimmed = record[int(args['start'] -1):aa_end]
        # export to fasta and 
        SeqIO.write(trimmed, "".join([str(fi).split(".")[0],"_",str(args['chain']),"_",str(args['start']),"_",str(aa_end),".fasta"]), "fasta")
        # remove intermediate file
        os.system("rm out.fasta")
        return
    # select between importing 1 or many pdb files
    if args['type'] == 1:
        trim_to_fasta(args['input'])
    else:
         # import each pdb file from the working directory
        for filename in sorted(os.listdir(os.chdir(args['directory']))):
            if filename.endswith(".pdb"):
                trim_to_fasta(filename)

if __name__ == '__main__':
    main()