# python3
from gooey import *
from Bio import SeqIO
# imput parameters
@Gooey(required_cols=5, program_name=' fasta subset by length', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
  ap = GooeyParser(description="filters a fasta file by keeping sequences under a rage of length")
  ap.add_argument("-in", "--input", required=True,widget='FileChooser', help="input fasta file")
  ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
  ap.add_argument("-max", "--max", required=True, help="max number of nucleotide/protein length")
  ap.add_argument("-min", "--min", required=True, help="min number of nucleotide/protein length")
  ap.add_argument("-headers", "--headers", required=True, widget='FileSaver', help="file to save the output fasta headers")
  args = vars(ap.parse_args())

# main
  sequences = []  # setup an empty list
  for record in SeqIO.parse(args['input'], "fasta"):
      if int(args['min']) < len(record.seq) < int(args['max']):
        # add this record to the list
         sequences.append(record)


  SeqIO.write(sequences, args['output'], "fasta")
# retrieve headers only
  headers = []  # setup an empty list
  for record in SeqIO.parse(args['input'], "fasta"):
      if int(args['min']) < len(record.seq) < int(args['max']):
        # add this record to the list
        headers.append(record.id)
  with open(args['headers'], 'w') as filehandle:
      for listitem in headers:
          filehandle.write('%s\n' % listitem)

if __name__ == '__main__':
    main()