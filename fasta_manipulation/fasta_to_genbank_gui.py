# python3
from gooey import *
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
# imput parameters
@Gooey(required_cols=2, program_name='genbank to fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-fa","--fasta", required=True, widget='FileChooser', help="input fasta file")
    ap.add_argument("-gb", "--genbank", required=True, widget='FileSaver', help="output genbank file")
    args = vars(ap.parse_args())
# main
    input_handle = open(args['fasta'], "rU")
    output_handle = open(args['genbank'], "w")
# import fasta
    sequences = list(SeqIO.parse(input_handle, "fasta"))

# asign generic_dna or generic_protein
    for seq in sequences:
        seq.seq.alphabet = generic_dna
# output
    count = SeqIO.write(sequences, output_handle, "genbank")

    output_handle.close()
    input_handle.close()

if __name__ == '__main__':
        main()
        
