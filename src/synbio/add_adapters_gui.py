# python3
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# imput parameters
@Gooey(required_cols=2, program_name='add adapters', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description = "add a left or right adapter or both in a single or multi-fasta file")
    ap.add_argument("-i", "--input", required=True, widget='FileChooser',  help="input single or multi fasta file")
    ap.add_argument("-l", "--left", required=False, default="", type=str, help="adapter to the left of the sequence")
    ap.add_argument("-r", "--right", required=False, default="", type=str, help="adapter to the right of the sequence")
    ap.add_argument("-o", "--output", required=True, widget='FileSaver',  help="output fasta file")
    args = vars(ap.parse_args())
    # main
    records = [] # setup an empty list
    for record in SeqIO.parse(args['input'], "fasta"):
        # add this record to the list
        records.append(SeqRecord(Seq(''.join([args['left'],str(record.seq),args['right']])),id=record.id,description=record.description))
    # export to fasta
    SeqIO.write(records, args['output'], "fasta")

if __name__ == '__main__':
    main()