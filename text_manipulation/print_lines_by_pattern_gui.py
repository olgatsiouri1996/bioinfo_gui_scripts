# python3
from gooey import *
import pandas as pd
# input parameters
@Gooey(required_cols=4, program_name='print lines by pattern', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
  ap = GooeyParser(description= "print lines that match a user specified pattern")
  ap.add_argument("-in", "--input", required=True, widget='FileChooser', help=" input txt/tsv file")
  ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output txt/tsv file")
  ap.add_argument("-p", "--pattern", required=True, help="specify the pattern to print the lines that have it")
  ap.add_argument("-c", "--column", required=True, help="specify the column number to search for the pattern(starts from 0)")
  args = vars(ap.parse_args())
# main
  df = pd.read_csv(args['input'], sep= "\t", encoding='latin-1')
  df = df.fillna("")
  bool2 = df.iloc[:, int(args['column'])].str.contains(args['pattern']) 
  df = df[bool2]
# export
  with open(args['output'], 'a') as f:
    f.write(
        df.to_csv(header = True, index = False, sep= "\t", doublequote= False, line_terminator= '\n')
    )

if __name__ == '__main__':
    main()
