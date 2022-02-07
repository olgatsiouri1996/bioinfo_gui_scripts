# python3
from gooey import *
import pandas as pd
# input parameters
@Gooey(required_cols=4, program_name='intersect tables', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
  ap = GooeyParser(description= "exports the common lines of  2 txt files")
  ap.add_argument("-tab1", "--table1", required=True, widget='FileChooser', help=" input txt file")
  ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output txt file")
  ap.add_argument("-tab2", "--table2", required=True, widget='FileChooser', help="txt file to match")
  ap.add_argument("-header", "--header", required=True, help="txt header to select data for merge")
  args = vars(ap.parse_args())
# main
  df1 = pd.read_csv(args['table1'], sep= "\t")
  df2 = pd.read_csv(args['table2'], sep= "\t")
  df_merge_col = pd.merge(df1, df2, on= args['header'], how='inner')
# export
  with open(args['output'], 'a') as f:
    f.write(
        df_merge_col.to_csv(header = True, index = False, sep= "\t", doublequote= False, line_terminator= '\n')
    )

if __name__ == '__main__':
    main()
