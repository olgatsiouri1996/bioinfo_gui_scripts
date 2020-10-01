# python3
from gooey import *
from Bio import SeqIO
# input parameters
@Gooey(required_cols=5, program_name='subset multifasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
  ap = GooeyParser()
  ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
  ap.add_argument("-start", "--start", required=True, help="region to start writing the fasta file(starts from 0)")
  ap.add_argument("-stop", "--stop", required=True, help="region to stop writing the fasta file(it can also be a negative number to remove nucleotides from the end of the sequence)")
  ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
  args = vars(ap.parse_args())
# main
  sequences = []  # setup an empty list
  for record in SeqIO.parse(args['input'], "fasta"):
        # add this record to the list
      sequences.append(record[int(args['start']):int(args['stop'])])

  SeqIO.write(sequences, args['output'], "fasta")

if __name__ == '__main__':
    main()