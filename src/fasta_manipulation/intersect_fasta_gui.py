# python3
from gooey import *
from pyfaidx import Fasta
# imput parameters
@Gooey(required_cols=3, program_name='intersect fasta',default_size=(1000, 530), fin_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-fa1", "--fasta1", required=True, widget='FileChooser', help="input multi fasta file")
    ap.add_argument("-fa2", "--fasta2", required=True, widget='FileChooser', help="input multi fasta file")
    ap.add_argument("-pro", "--program", type=str, default='export the sequences of fasta1 that have the same identifiers as fasta2', choices=['export the sequences of fasta1 that have the same identifiers as fasta2','export the sequences of fasta1 that do not have the same identifiers as fasta2'], required=False, help="program to choose")
    ap.add_argument("-type", "--type", type=str, default='1 multi-fasta file', choices=['1 multi-fasta file','2-column tab seperated txt file with id and seq as columns','2-column tab seperated txt file with id and description as columns','3-column tab seperated txt file with id, description and seq as columns'], required=False, help="output type to export to")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver',  help="output multi fasta or 2 or 3-column tab seperated txt file")
    args = vars(ap.parse_args())
    # main
    # helper function to wrap fasta sequence to 60 characters per line
    def wrap_fasta_seq(seq):
        return '\n'.join([seq[i:i+60] for i in range(0, len(seq), 60)])
    # create fasta index
    features1 = Fasta(args['fasta1'])
    features2 = Fasta(args['fasta2'])
    ## choose program
    if args['program'] == 'export the sequences of fasta1 that have the same identifiers as fasta2':
        # find common ids of the 2 files
        final = (set(features1.keys()).intersection(features2.keys()))
    else:
        final = (set(features1.keys()).difference(features2.keys()))
    # choose export type
    # choose program
    type = args['type']
    match type:
        case '1 multi-fasta file':
            # export to fasta
            with open(args['output'], 'w') as f:
                for fin in final:
                    f.write(f'>{str(features1[str(fin)].long_name).rstrip()}\n{wrap_fasta_seq(features1[str(fin)][:].seq)}\n')
        case '2-column tab seperated txt file with id and seq as columns':
            with  open(args['output'], 'w') as f:
                    f.write(f'{"id"}\t{"seq"}\n')
                    for fin in final:
                        f.write(f'{str(fin)}\t{features1[str(fin)][:].seq}\n')
        case '2-column tab seperated txt file with id and description as columns':
            with  open(args['output'], 'w') as f:
                    f.write(f'{"id"}\t{"description"}\n')
                    for fin in final:
                        try:
                            f.write(f'{str(fin)}\t{str(str(features1[str(fin)].long_name).rstrip()).split(" ",1)[1]}\n')
                        except IndexError:
                            f.write(f'{str(fin)}\t{""}\n')
        case '3-column tab seperated txt file with id, description and seq as columns':
            with  open(args['output'], 'w') as f:
                    f.write(f'{"id"}\t{"description"}\t{"seq"}\n')
                    for fin in final:
                        try:
                            f.write(f'{str(fin)}\t{str(str(features1[str(fin)].long_name).rstrip()).split(" ",1)[1]}\t{features1[str(fin)][:].seq}\n')
                        except IndexError:
                            f.write(f'{str(fin)}\t{""}\t{features1[str(fin)][:].seq}\n')

if __name__ == '__main__':
    main()