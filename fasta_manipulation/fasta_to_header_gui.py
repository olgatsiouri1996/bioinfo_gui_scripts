# python3
from gooey import *
from Bio import SeqIO
# input arguments
@Gooey(required_cols=2, program_name='fasta to header', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
  ap = GooeyParser(description="exctract the headers from a multifasta file")
  ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
  ap.add_argument("-headers", "--headers", required=True, widget='FileSaver', help="file to save the output fasta headers")
  args = vars(ap.parse_args())
  # main
  headers = []
  for record in SeqIO.parse(args['input'], "fasta"):
      headers.append(record.id)
  with open(args['headers'], 'w') as filehandle:
      for listitem in headers:
          filehandle.write('%s\n' % listitem)

if __name__ == '__main__':
    main()