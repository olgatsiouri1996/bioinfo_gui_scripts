# python3
from gooey import *
from Bio import SeqIO
import pandas as pd
# input parameters
@Gooey(required_cols=2, program_name= 'count length from fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
  ap = GooeyParser(description="extract the id and length of sequence from a fasta or multifasta file")
  ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
  ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output txt file")
  args = vars(ap.parse_args())
# main
  headers = []
  lengths = [] # setup  empty lists
  for record in SeqIO.parse(args['input'], "fasta"):
        # add this record to the lists
     lengths.append(len(record.seq))
     headers.append(record.id)
# create data frame
  df = pd.DataFrame()
  df['id'] = headers
  df['length'] = lengths
# export
  with open(args['output'], 'a') as f:
    f.write(
        df.to_csv(header = True, index = False, sep= "\t", doublequote= False, line_terminator= '\n')
    )

if __name__ == '__main__':
    main()
