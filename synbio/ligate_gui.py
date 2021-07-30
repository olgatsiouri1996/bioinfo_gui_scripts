# python3
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# imput parameters
@Gooey(required_cols=3, program_name='digital ligation', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="ligate linear vector with DNA insert")
    ap.add_argument("-vr", "--vector", required=True, widget='FileChooser', help="linear vector(fasta format)")
    ap.add_argument("-in", "--insert", required=True, widget='FileChooser', help="sequence/s to insert in the vector(single or multi fasta file)")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
    args = vars(ap.parse_args())
# main 
# linear vector
    for record in SeqIO.parse(args['vector'], "fasta"):
        x = str(record.seq)
# DNA insert
    records = [] # setup an empty list
    for record in SeqIO.parse(args['insert'], "fasta"):
        y = str(record.seq)
# merge
        seqad = x + y
        # add this record to the list
        records.append(SeqRecord(Seq(seqad),id=record.id,description=""))
# export to fasta
    SeqIO.write(records, args['output'], "fasta")

if __name__ == '__main__':
    main()
