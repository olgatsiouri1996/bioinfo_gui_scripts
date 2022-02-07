# python3
from gooey import *
import sys
import COG
from more_itertools import flatten
# input parameters
@Gooey(required_cols=2, program_name='convert COGs to functional category descriptions', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description= "converts COG letters from a 1-column txt file to their functional category descriptions in tabular format")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input txt file with each COG in each line")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output tab seperated txt file with each COGs and their functional category descriptions")
    args = vars(ap.parse_args())
    # main
    sys.stdout = open(args['output'], 'a')    
    with open(args['input'], 'r') as f:
        for line in f:
            spl = line.split()
            func_tuple = COG.cat_from_letter(spl[0])
            func_tuple_flat = list(flatten(func_tuple))
            annot = list(flatten(func_tuple_flat))
            print('\t'.join(annot))
    sys.stdout.close()

if __name__ == '__main__':
    main()
    