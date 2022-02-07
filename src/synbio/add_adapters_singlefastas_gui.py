# python3
import os
from gooey import *
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# imput parameters
@Gooey(required_cols=3, program_name='add adapters to single-fasta files', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-up", "--upstream", required=True, widget='FileChooser',  help="upstream adapter(fasta format)")
    ap.add_argument("-down", "--downstream", required=True, widget='FileChooser',  help="downstream adapter(fasta format)")
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser', help="directory to search for fasta files")
    args = vars(ap.parse_args())
    # main
    # first adapter
    x = SeqIO.read(args['upstream'], "fasta")
    # second adapter
    z = SeqIO.read(args['downstream'], "fasta")
    # sequence
    # import each fasta file from the working directory
    for filename in sorted(os.listdir(os.chdir(args['directory']))):
        if filename.endswith(".fa") or filename.endswith(".fasta"):
            y = SeqIO.read(filename, "fasta")
            # merge
            seqad = str(x.seq) + str(y.seq) + str(z.seq)
                # create SeqRecord
            record = SeqRecord(Seq(seqad),id=y.id,description="")
            # export to fasta
            SeqIO.write(record, "".join([filename.split(".")[0],"_","with","_","ad",".fasta"]), "fasta")

if __name__ == '__main__':
    main()
