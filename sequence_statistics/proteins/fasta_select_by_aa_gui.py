# python3
from gooey import *
from Bio import SeqIO
# input parameters
@Gooey(required_cols=5, program_name=' extract sequences by %aa content', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
  ap = GooeyParser(description="selects sequences from a fasta file based on a range of %aminoacid content")
  ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
  ap.add_argument("-max", "--max", required=True, help="max threshold of aa content, type = float")
  ap.add_argument("-min", "--min", required=True, help="min threshold of aa content, type = float")
  ap.add_argument("-aa", "--aa", required=True, help="aa to search the content for")
  ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
  args = vars(ap.parse_args())
# create aa_content function
  def aa_content(seq):
    return round((seq.count(args['aa']) / len(seq)) * 100, 2)
# main
  sequences = []  # setup an empty list
  for record in SeqIO.parse(args['input'], "fasta"):
    if float(args['min']) < aa_content(record.seq) < float(args['max']):
        # add this record to the list
        sequences.append(record)

  SeqIO.write(sequences, args['output'], "fasta")

if __name__ == '__main__':
    main()