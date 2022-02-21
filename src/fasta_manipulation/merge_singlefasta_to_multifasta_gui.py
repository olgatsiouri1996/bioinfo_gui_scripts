# python3
from gooey import *
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# input parameters
@Gooey(required_cols=2, program_name= 'merge single-fasta to multi-fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-dir", "--directory", required=True, widget='DirChooser', help="input directory with fasta files")
    ap.add_argument("-fa", "--multifasta", required=True, widget='FileSaver', help="output multi-fasta file")
    args = vars(ap.parse_args())
# main
# creat list
    records = []
# import each fasta file from the working directory
    for filename in sorted(os.listdir(os.chdir(args['directory']))):
        if filename.endswith(".fa") or filename.endswith(".fasta"):
            record = SeqIO.read(filename, "fasta")
            records.append(record)
# export all SeqRecords to a multi-fasta file
    SeqIO.write(records,args['multifasta'], "fasta")

if __name__ == '__main__':
    main()
