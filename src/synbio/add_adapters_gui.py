# python3
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# imput parameters
@Gooey(required_cols=4, program_name='add adapters to sequence', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-up", "--upstream", required=True, widget='FileChooser', help="upstream adapter(fasta format)")
    ap.add_argument("-fa", "--fasta", required=True, widget='FileChooser', help="input single or multi fasta file with the sequence/s to add adapters")
    ap.add_argument("-down", "--downstream", required=True, widget='FileChooser', help="downstream adapter(fasta format)")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
    args = vars(ap.parse_args())
# main
# first adapter
    for record in SeqIO.parse(args['upstream'], "fasta"):
        x = str(record.seq)
# second adapter
    for record in SeqIO.parse(args['downstream'], "fasta"):
        z = str(record.seq)
# sequence
    records = [] # setup an empty list
    for record in SeqIO.parse(args['fasta'], "fasta"):
        y = str(record.seq)
# merge
        seqad = x + y + z
        # add this record to the list
        records.append(SeqRecord(Seq(seqad),id=record.id,description=""))
# export to fasta
    SeqIO.write(records, args['output'], "fasta")

if __name__ == '__main__':
    main()
