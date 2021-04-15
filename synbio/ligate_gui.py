# python3
from gooey import *
from Bio import SeqIO
import sys
# imput parameters
@Gooey(required_cols=4, program_name='digital ligation', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="ligate linear vector with DNA insert")
    ap.add_argument("-vr", "--vector", required=True, widget='FileChooser', help="linear vector(fasta format)")
    ap.add_argument("-in", "--insert", required=True, widget='FileChooser', help="sequence to insert in the vector(fasta format)")
    ap.add_argument("-id", "--seqid", required=True, help="fasta header of the output file")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
    args = vars(ap.parse_args())
# main 
# linear vector
    for record in SeqIO.parse(args['vector'], "fasta"):
        x = str(record.seq)
# DNA insert
    for record in SeqIO.parse(args['insert'], "fasta"):
        y = str(record.seq)
# merge
    seqad = x + y
# output
    sys.stdout = open(args['output'], 'a')
    print(">"+args['seqid'], seqad, sep='\n')
    sys.stdout.close()

if __name__ == '__main__':
    main()
