# python3
from gooey import *
from Bio import SeqIO
# imput parameters
@Gooey(required_cols=2, program_name='genbank to fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-gb", "--genbank", required=True, widget='FileChooser', help="input genbank file")
    ap.add_argument("-fa","--fasta", required=True, widget='FileSaver', help="output fasta file")
    args = vars(ap.parse_args())
# main
    count = SeqIO.convert(args['genbank'], "genbank", args['fasta'], "fasta")

if __name__ == '__main__':
    main()

