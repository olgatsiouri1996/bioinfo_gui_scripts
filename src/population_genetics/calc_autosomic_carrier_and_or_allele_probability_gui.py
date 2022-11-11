# python3
from gooey import *
import sys
# input parameters
@Gooey(required_cols=2, program_name='calculate the autosomic allele and/or carrier probability', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description='calculate the probability of an individual carrying at least 1 autosomic allele copy and/or being heterozygous, using the proportion of homozygous individuals in a given population')
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="1-column input txt file with the proportion of individuals homozygous for the specific allele")
    ap.add_argument("-type", "--type", required=False, type=str, default='allele', choices=['allele', 'carrier', 'both'], help="type of calculation")
    ap.add_argument("-all", "--allele", required=False, type=str, default='this', choices=['this', 'other'], help="allele to calculate the probability of for each calculation type")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="2 or 3-column tab-seperated output txt file with the probability of containing an allele and/or being heterozygous and the input column(or the proportion of homozygous individuals for the other allele)")
    args = vars(ap.parse_args())
    # convert 1-column txt file to list
    with open(args['input'],'r') as f:
        A = f.readlines()
    A = (float(x.strip()) for x in A)
    if args['allele'] == 'this':
        pass
    else:
        A = ((1-(y**0.5))**2 for y in A)
    # choose an output based on the calculation type
    sys.stdout = open(args['output'],'w')
    if args['type'] == 'allele':
        print('probability_of_allele','proportion_of_homozygous_individuals',sep='\t',end='\n')
        for i in A:
            print("%.3f" % round(2*i**0.5-i,3),"%.3f" % round(i,3),sep='\t',end='\n')
    elif args['type'] == 'carrier':
        print('probability_of_carriers','proportion_of_homozygous_individuals',sep='\t',end='\n')
        for i in A:
            p = i**0.5
            print("%.3f" % round(2*p*(1-p),3),"%.3f" % round(i,3),sep='\t',end='\n')
    else:
        print('probability_of_allele','probability_of_carriers','proportion_of_homozygous_individuals',sep='\t',end='\n')
        for i in A:
            p = i**0.5
            print("%.3f" % round(2*i**0.5-i,3),"%.3f" % round(2*p*(1-p),3),"%.3f" % round(i,3),sep='\t',end='\n')

    sys.stdout.close()

if __name__ == '__main__':
    main()