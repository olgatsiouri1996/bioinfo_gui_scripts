# python3
from gooey import *
import sys
# input parameters
@Gooey(required_cols=2, program_name='calculate reccesive allele propability', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description='calculate the probability of an individual carrying a reccesive allele, using the proportion of homozygous reccesive individuals in a given population')
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="1-column input txt file with the proportion of individuals homozygous for the reccesive allele")
    ap.add_argument("-sex", "--sex", required=False, type=str, default='no', choices=['no', 'yes'], help="linked for the X chromosome?")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="2-column tab-seperated output txt file with the probability of containing a reccesive allele and the input column")
    args = vars(ap.parse_args())
    # convert 1-column txt file to list
    with open(args['input'],'r') as f:
        A = f.readlines()
    A = [float(x.strip()) for x in A]
    # calculate the probability of an individual containing a reccesive allele and output both columns
    sys.stdout = open(args['output'],'w')
    print('carries_reccesive','proportion_of_homozygous_reccesive',sep='\t',end='\n')
    if args['sex'] == 'no':
        for i in A:
            print("%.3f" % round(2*i**0.5-i,3),"%.3f" % round(i,3),sep='\t',end='\n')
    else:
        for i in A:
            print("%.3f" % round(2*i*(1-i),3),"%.3f" % round(i,3),sep='\t',end='\n')

    sys.stdout.close()

if __name__ == '__main__':
    main()