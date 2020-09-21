# python3
from gooey import *
from Bio import SeqIO
# input parameters
@Gooey(required_cols=2, program_name='tab to fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
  ap = GooeyParser(description="converts a tabular file with identifier and sequence, into a fasta file")
  ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input tab seperated file")
  ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")                                                                  
  args = vars(ap.parse_args())
  # main
  count = SeqIO.convert(args['input'], "tab", args['output'], "fasta")

if __name__ == '__main__':
    main()