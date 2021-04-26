# python3
from gooey import *
from Bio import SeqIO
import sys
# imput parameters
@Gooey(required_cols=2, program_name='merge multifasta to singlefasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-mfa", "--multifasta", required=True, widget='FileChooser', help="input multi-fasta file to merge its sequences")
    ap.add_argument("-id", "--seqid", required=True, help="fasta header of the output file")
    ap.add_argument("-sfa", "--singlefasta", required=True, widget='FileSaver', help="output single-fasta file")
    args = vars(ap.parse_args())
# main
    sequences = []  # setup an empty list
    for record in SeqIO.parse(args['multifasta'], "fasta"):
        sequences.append(record.seq)
# output
    sys.stdout = open(args['singlefasta'], 'a')
    print(">"+args['seqid'], ''.join(map(str,sequences)), sep='\n')
    sys.stdout.close()

if __name__ == '__main__':
     main()

        
      
      
