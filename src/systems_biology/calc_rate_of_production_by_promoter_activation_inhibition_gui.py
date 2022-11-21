# python3
from gooey import *
import sys
# input parameters
@Gooey(required_cols=5, program_name='calculate the rate of protein production by promoter activation or inhibition', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description='calculate the rate of protein production by promoter activation or inhibition based on the input function')
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="1-column input txt file with the concentration of the transcriptional activator or repressor")
    ap.add_argument("-type", "--type", required=False, type=str, default='activator', choices=['activator', 'repressor'], help="type of transcriptional regulator")
    ap.add_argument("-beta", "--beta", required=True, type=float, help="maximal promoter activity")
    ap.add_argument("-kappa", "--kappa", required=True, type=float, help="activation or repression coefficient")
    ap.add_argument("-hill", "--hill", required=True, type=int, help="hill coefficient")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="2-column tab-seperated output txt file with the rate of protein production and the input column")
    args = vars(ap.parse_args())
    # inport variables
    b = args['beta']
    k = args['kappa']
    n = args['hill']
    # create function to calc the production rate by activation
    def activation_func(s):
        return (b*(s**n))/(k**n + s**n)
    # create function to calc the production rate by repression
    def repression_func(s):
        return (b*(k**n))/(k**n + s**n)
    # convert 1-column txt file to list
    tfs = (float(line.rstrip()) for line in open(args['input']))
    # choose an output based on the calculation type
    sys.stdout = open(args['output'],'w')
    if args['type'] == 'activator':
        print('rate_of_protein_production','activator_concentration',sep='\t',end='\n')
        for tf in tfs:
            print("%.3f" % round(activation_func(tf),3),"%.3f" % round(tf,3),sep='\t',end='\n')
    else:
        print('rate_of_protein_production','repressor_concentration',sep='\t',end='\n')
        for tf in tfs:
            print("%.3f" % round(repression_func(tf),3),"%.3f" % round(tf,3),sep='\t',end='\n')
    
    sys.stdout.close()

if __name__ == '__main__':
    main()