# python3
import os
from gooey import *
from Bio import SeqIO
import sys
# input parameters
@Gooey(required_cols=0, program_name= 'fasta formatter', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="changes the width of sequences line in 1 or many FASTA files")
    ap.add_argument("-in", "--input", required=False, widget="FileChooser", help="input fasta file")
    ap.add_argument("-out", "--output", required=False, widget="FileSaver", help="output fasta file")
    ap.add_argument("-dir", "--directory", required=False, type=str, widget='DirChooser', help="directory to search for fasta files")
    ap.add_argument("-width", "--width", required=False, type=int, default=80, help="number of characters per line. Default 80")
    ap.add_argument("-pro", "--program", required=False, type=int, default=1, help="program to choose. 1) one input/output fasta file, 2) many input/output fasta files. Default is 1")
    args = vars(ap.parse_args())
    # main
    # create function to split the input sequence based on a specific number of characters
    def split_every_width(s,w): return [s[i:i+w] for i in range(0,len(s),w)]
    # choose program
    if args['program'] == 1:
        # export to a new fasta file
        sys.stdout = open(args['output'], 'a')
        for record in SeqIO.parse(args['input'],'fasta'):
                print(">"+record.id)
                print('\n'.join(split_every_width(str(record.seq), args['width']))) # add characters in new line after the number of characters surpasses the input width 
        sys.stdout.close()
    else:
        # import each fasta file from the working directory
        for filename in sorted(os.listdir(os.chdir(args['directory']))):
            if filename.endswith(".fa") or filename.endswith(".fasta"):
                # export to new fasta files with the user imported width value
                sys.stdout = open(''.join([filename.split(".")[0],"_","w",str(args['width']),".fasta"]), 'a')
                for record in SeqIO.parse(filename,'fasta'):
                        print(">"+record.id)
                        print('\n'.join(split_every_width(str(record.seq), args['width']))) # add characters in new line after the number of characters surpasses the input width 
                sys.stdout.close()

if __name__ == '__main__':
    main()
