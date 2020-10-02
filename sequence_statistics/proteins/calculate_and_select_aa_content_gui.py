# python3
from gooey import *
from Bio import SeqIO
import pandas as pd
# input parameters
@Gooey(required_cols=5, program_name='calculate and select %aa content', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
  ap = GooeyParser(description="calculate a specific %aminoacid content into a selected range and output it in a tabular file")
  ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
  ap.add_argument("-max", "--max", required=True, help="max threshold of aa content, type = float")
  ap.add_argument("-min", "--min", required=True, help="min threshold of aa content, type = float")
  ap.add_argument("-aa", "--aa", required=True, help="amino acid to search the content for")
  ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output txt file")
  args = vars(ap.parse_args())
# create aa_content function
  def aa_content(seq):
    return round((seq.count(args['aa']) / len(seq)) * 100, 2)
# main
  content = []
  headers = []  # setup empty lists
  for record in SeqIO.parse(args['input'], "fasta"):
    if float(args['min']) < aa_content(record.seq) < float(args['max']):
        # add this record to the list
        headers.append(record.id)
        content.append(aa_content(record.seq))
# create data frame
  df = pd.DataFrame()
  df['protein_id'] = headers
  df[args['aa']] = content
# export
  with open(args['output_file'], 'a') as f:
    f.write(
        df.to_csv(header = True, index = False, sep= "\t")
    )

if __name__ == '__main__':
    main()
