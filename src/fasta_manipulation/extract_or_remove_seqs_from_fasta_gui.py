# python3
from gooey import *
import os
from pyfaidx import Fasta
import textwrap
# input parameters
@Gooey(required_cols=2, program_name= 'extract or remove sequences from fasta',default_size=(1030, 530), header_bg_color= '#F5F5F5', body_bg_color='#F5F5F5', terminal_font_color= '#F5F5F5', terminal_panel_color= '#F5F5F5')
def main():
    ap = GooeyParser(description="use a txt file with fasta headers to extract or remove sequences from fasta file")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input multi-fasta file")
    ap.add_argument("-ids", "--ids", required=True, widget='FileChooser', help="1 or multiple column tab seperated file with no column names, with fasta identifiers in the 1st column to retrieve the output fasta sequences")
    ap.add_argument("-act", "--action",type=str, default='extract', required=False, choices=['extract','remove'], widget='Dropdown', help="choose to extract or remove sequences")
    ap.add_argument("-ty", "--type",type=str, default='1 multi-fasta file', required=False, choices=['1 multi-fasta file','many single-fasta files','2-column txt file(id,seq)','3-column txt file(id,description,seq)','3-column txt file(id,full_fasta_header,seq)','2-column txt file(id,description)','1-column txt file(id)'], widget='Dropdown', help="output type")
    ap.add_argument("-out", "--output", required=False, widget='FileSaver', help="output multi-fasta or a 1-column or a 2-column or a 3-column txt file")
    ap.add_argument("-dir", "--directory", required=False, type=str, widget='DirChooser',  help="output directory to save the single-fasta files")
    args = vars(ap.parse_args())
    # main
    # import the txt file with headers you want to extract the sequence from the input fasta
    headers = (str(line.rstrip()).split()[0] for line in open(args['ids']))
    # import fasta file
    features = Fasta(args['input'])
    # choose to remove sequences
    if 'remove' in args['action']:
        # collect ids of import fasta file
        keyslist = (features.keys())
        # remove those that exist in the headers generator
        final_keys = (set(keyslist).difference(headers))
    else:
        final_keys = headers
    # choose program
    program = args['type']
    match program:
        case '1 multi-fasta file':
            # export to 1 multi-fasta
            with  open(args['output'], 'w') as f:
                for key in final_keys:
                    f.write(f'>{str(features[str(key)].long_name).rstrip()}\n{textwrap.fill(features[str(key)][:].seq, width=60)}\n')                
        case 'many single-fasta files':
            # extract many sigle-fasta files
            os.chdir(args['directory'])
            for key in final_keys:
                with open(str(key)+".fasta", 'w') as f:
                    f.write(f'>{str(features[str(key)].long_name).rstrip()}\n{textwrap.fill(features[str(key)][:].seq, width=60)}\n')
                # extract a multi-fasta file
        case '2-column txt file(id,seq)':                        
            # iterate input headers to extract sequences and export as 2 column txt
            with  open(args['output'], 'w') as f:
                f.write(f'{"id"}\t{"seq"}\n')
                for key in final_keys:
                    f.write(f'{str(key)}\t{features[str(key)][:].seq}\n')
        case '3-column txt file(id,description,seq)':
            # iterate input headers to extract sequences and export as 3 column txt
            with  open(args['output'], 'w') as f:
                f.write(f'{"id"}\t{"description"}\t{"seq"}\n')
                for key in final_keys:
                    try:
                        f.write(f'{str(key)}\t{str(str(features[str(key)].long_name).rstrip()).split(" ",1)[1]}\t{features[str(key)][:].seq}\n')
                    except IndexError:
                        f.write(f'{str(key)}\t{""}\t{features[str(key)][:].seq}\n')
        case '3-column txt file(id,full_fasta_header,seq)':           
            # iterate input headers to extract sequences and export as 3 column txt
            with  open(args['output'], 'w') as f:
                f.write(f'{"id"}\t{"full_fasta_header"}\t{"seq"}\n')
                for key in final_keys:
                    f.write(f'{str(key)}\t{str(features[str(key)].long_name).rstrip()}\t{features[str(key)][:].seq}\n')
        case '2-column txt file(id,description)':
            # iterate input headers to extract sequences and export as 2 column txt
            with  open(args['output'], 'w') as f:
                f.write(f'{"id"}\t{"description"}\n')
                for key in final_keys:
                    try:
                        f.write(f'{str(key)}\t{str(str(features[str(key)].long_name).rstrip()).split(" ",1)[1]}\n')
                    except IndexError:
                        f.write(f'{str(key)}\t{""}\n')
        case '1-column txt file(id)':
            # iterate input headers to extract sequences and export as 1-column txt
            with  open(args['output'], 'w') as f:
                f.write(f'{"id"}\n')
                for key in final_keys:
                    f.write(f'{str(key)}\n')
                                        
if __name__ == '__main__':
    main()
            