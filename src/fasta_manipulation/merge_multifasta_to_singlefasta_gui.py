# python3
from gooey import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# imput parameters
@Gooey(required_cols=3, program_name='merge multifasta to singlefasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-mfa", "--multifasta", required=True, widget='FileChooser', help="input multi-fasta file to merge its sequences")
    ap.add_argument("-id", "--seqid", required=True, help="fasta header of the output file")
    ap.add_argument("-sp", "--spacer", required=False, type=str, default="", help="nucleotides or aminoacids to add between the merged fasta sequences. Default: no sequence to add")
    ap.add_argument("-sfa", "--singlefasta", required=True, widget='FileSaver', help="output single-fasta file")
    args = vars(ap.parse_args())
# main
    sequences = []  # setup an empty list
    for record in SeqIO.parse(args['multifasta'], "fasta"):
        sequences.append(record.seq)
# merge all sequences at 1 and create SeqRecord object
    merged_seqs = SeqRecord(Seq(args['spacer'].join(map(str,sequences))),id=args['seqid'],description="")
# export to fasta
    SeqIO.write(merged_seqs, args['singlefasta'], "fasta")

if __name__ == '__main__':
    main()
