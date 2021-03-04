# python3
from gooey import *
from Bio import SeqIO
from synbiopython.codon import table, taxonomy_utils, utils
import sys
# imput parameters
@Gooey(required_cols=3, program_name='codon optimization', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input multi or single fasta file")
    ap.add_argument("-taxid", "--taxid", required=True, help="taxonomy id to retrieve the codon table for optimization")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output single or multi fasta file")
    args = vars(ap.parse_args())
# main
    name = taxonomy_utils.get_organism_name(args['taxid'])
    name_table = table.get_table(name)
# output
    sys.stdout = open(args['output'], 'a')
    for record in SeqIO.parse(args['input'], "fasta"):
        print(">"+record.id,utils.optimise(name_table, str(record.seq)), sep='\n')
    sys.stdout.close()

if __name__ == '__main__':
    main()
