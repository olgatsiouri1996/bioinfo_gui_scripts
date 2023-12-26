# python3
from gooey import *
# input parameters
@Gooey(required_cols=2, program_name='calculate the autosomic allele and/or carrier probability', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description='calculate the probability of an individual carrying at least 1 autosomic allele copy and/or being heterozygous, using the probability of homozygous individuals in a given diploid population')
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="1-column input txt file with the probability of individuals homozygous for the specific allele")
    ap.add_argument("-type", "--type", required=False, type=str, default='allele', choices=['allele', 'carrier', 'both'], help="type of calculation")
    ap.add_argument("-all", "--allele", required=False, type=str, default='this', choices=['this', 'other'], help="allele to calculate the probability of for each calculation type")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="3 or 4-column tab-seperated output txt file with the allele frequency, the probability of containing an allele and/or being heterozygous and the input column(for this or the other allele)")
    args = vars(ap.parse_args())
    # convert 1-column txt file to list
    probabilities = (float(line.rstrip()) for line in open(args['input']))
    # Adjust probabilities based on allele choice
    if args.allele == 'other':
        probabilities = ((1 - (y ** 0.5)) ** 2 for y in probabilities)

    # Choose an output based on the calculation type
    with open(args.output, 'w') as output_file:
        if args.type == 'allele':
            output_file.write('probability_of_containing_allele\tprobability_of_homozygous_individuals\tallele_frequency\n')
            for i in probabilities:
                output_file.write(f"{round(2 * i**0.5 - i, 3)}\t{round(i, 3)}\t{round(i**0.5, 3)}\n")
        elif args.type == 'carrier':
            output_file.write('probability_of_carriers\tprobability_of_homozygous_individuals\tallele_frequency\n')
            for i in probabilities:
                p = i**0.5
                output_file.write(f"{round(2 * p * (1 - p), 3)}\t{round(i, 3)}\t{round(p, 3)}\n")
        else:
            output_file.write('probability_of_containing_allele\tprobability_of_carriers\tprobability_of_homozygous_individuals\tallele_frequency\n')
            for i in probabilities:
                p = i**0.5
                output_file.write(f"{round(2 * i**0.5 - i, 3)}\t{round(2 * p * (1 - p), 3)}\t{round(i, 3)}\t{round(p, 3)}\n")

if __name__ == '__main__':
    main()