# python3
from gooey import *
import pandas as pd
# input parameters
@Gooey(required_cols=2, program_name='intesect csv', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
  ap = GooeyParser(description= "exports the common lines of  2 csv files")
  ap.add_argument("-tab1", "--table1", required=True, widget='FileChooser', help=" input csv file")
  ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output csv file")
  ap.add_argument("-tab2", "--table2", required=True, widget='FileChooser', help="csv file to match")
  ap.add_argument("-header", "--header", required=True, help="csv header to select data for merge")
  args = vars(ap.parse_args())
# main
  df1 = pd.read_csv(args['table1'])
  df2 = pd.read_csv(args['table2'])
  df_merge_col = pd.merge(df1, df2, on= args['header'], how='inner')
# export
  with open(args['output'], 'a') as f:
    f.write(
        df_merge_col.to_csv(header = True, index = False, doublequote= False)
    )

if __name__ == '__main__':
    main()