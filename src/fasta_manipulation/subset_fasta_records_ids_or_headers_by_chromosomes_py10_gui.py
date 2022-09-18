# python3
from gooey import *
import sys
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=2, program_name= 'extract fasta records, ids or headers by chromosomes/scaffolds or contings', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input multi-fasta file(if it has records with the same identifier only the 1st record is picked)")
    ap.add_argument("-ma", "--match", required=True, widget='FileChooser', help="1-column txt file with words to search in each fasta identifier e.g 5g in Lj5g3v1117110.1.1(no regular expression, differentiates between capital or non capital letters)")
    ap.add_argument("-out", "--output", required=False, widget='FileSaver', help="output multi-fasta file. File can be appended")
    ap.add_argument("-headers", "--headers", required=False, widget='FileSaver', help="1 or 2-column tab seperated txt file to save the output fasta identifiers or full fasta headers respectively. File can be appended")
    ap.add_argument("-pro", "--program", required=False, type=int, default=1, widget='Dropdown', choices=[1,2,3], help="Program to choose: 1) collect fasta records with headers that match the pattern 2) collect only fasta identifiers, 3) collect full fasta headers(id and description are tab seperated).")
    args = vars(ap.parse_args())
    # main
    # import the txt file with headers you want to extract the sequence from the input fasta
    with open(args['match'], 'r') as f:
        chroms = f.readlines()
    chroms = [x.strip() for x in chroms]
    # create function to split the input sequence based on a specific number of characters(60)
    def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
    # import fasta file
    try:
        features = Fasta(args['input'])
    except ValueError:
        features = Fasta(args['input'], duplicate_action="first")
    # choose program
    program = args['program']
    match program:
        case 1:
            # iterate input headers to extract sequences and export as multi-fasta
            sys.stdout = open(args['output'], 'a')
            for chrom in chroms:
                for key in features.keys():
                    if str(chrom) in key:
                        print(''.join([">",features[key].long_name]).replace('\r',''))
                        print('\n'.join(split_every_60(features[key][:].seq)))
            sys.stdout.close()
        case 2:
            # export to 1-column txt file
            with open(args['headers'], 'a') as filehandle:
                for chrom in chroms:
                    for key in features.keys():
                        if str(chrom) in key:
                            filehandle.write('%s\n' % key)
        case 3:
            # export to 1-column txt file
            with open(args['headers'], 'a') as filehandle:
                for chrom in chroms:
                    for key in features.keys():
                        if str(chrom) in key:
                            try:
                                description = str(features[key].long_name).split(' ',1)[1]
                            except IndexError:
                                description = ""

                            filehandle.write('%s\n' % '\t'.join([key,description]))

if __name__ == '__main__':
    main()
