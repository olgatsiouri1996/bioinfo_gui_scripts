# python3
import os
from gooey import *
from Bio import SeqIO
# input parameters
@Gooey(required_cols=1, program_name='trim multiple single-fasta files', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-start", "--start", required=False, default=1, type=int, help="region to start writing the fasta file")
    ap.add_argument("-stop", "--stop", required=False, type=int, help="region to stop writing the fasta file(it can be both a positive and  a negative number)")
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser', help="directory to search for fasta files")
    ap.add_argument("-pro", "--program", required=False,default=1, type=int, help="program to choose 1) add both start and stop location 2) the stop location with be that of the sequence length. Default is 1")
    args = vars(ap.parse_args())
# main
# create function to trim fasta records
    def fastatrim(fastarec,fastaseq):
        # choose program
        if args['program'] == 1:
            # fix the index for start parameter
            if args['start'] > 0:
                seq_start = args['start'] -1
            else:
                print("-start parameter must be a positive integer")
                exit(1)
            # add end parameter
            seq_end = args['stop']
        else:
            # fix the index for start parameter
            if args['start'] > 0:
                seq_start = args['start'] -1
            else:
                print("-start parameter must be a positive integer")
                exit(1)
            # add end parameter according to program 2
            args['stop'] = len(fastaseq)
            seq_end = args['stop']
        # subset each fasta record
        return fastarec[seq_start:seq_end]
# import each fasta file from the working directory
    for filename in sorted(os.listdir(os.chdir(args['directory']))):
        if filename.endswith(".fa") or filename.endswith(".fasta"):
            # read each file, trim and create SeqRecord to export
            record = SeqIO.read(filename, "fasta")
            sequence = fastatrim(record,record.seq)
            # export to fasta
            SeqIO.write(sequence, "".join([filename.split(".")[0],"_","trimmed",".fasta"]), "fasta")

if __name__ == '__main__':
    main()