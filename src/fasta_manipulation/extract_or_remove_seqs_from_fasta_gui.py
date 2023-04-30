# python3
from gooey import *
import os
from pyfaidx import Fasta
# input parameters
@Gooey(required_cols=2, program_name= 'extract or remove sequences from fasta',default_size=(870, 530), header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="use a txt file with fasta headers to extract or remove sequences from fasta file")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input multi-fasta file")
    ap.add_argument("-ids", "--ids", required=True, widget='FileChooser', help="file with fasta headers to retrieve the output fasta sequences")
    ap.add_argument("-pro", "--program",type=str, default='extract sequences from a multi-fasta file', required=False, choices=['extract sequences from a multi-fasta file','extract many single-fasta files','remove sequences from a multi-fasta file','remove sequences and export to many single-fasta files','extract sequences to a 2-column txt file','remove sequences and export the rest to a 2-column txt file','extract sequences to a 3-column txt file','remove sequences and export the rest to a 3-colmn txt file'], widget='Dropdown', help="program to choose")
    ap.add_argument("-out", "--output", required=False, widget='FileSaver', help="output multi-fasta or a 2-column or a 3-column txt file with id, seq or id, description, seq as headers respectively")
    ap.add_argument("-dir", "--directory", required=False, type=str, widget='DirChooser',  help="output directory to save the single-fasta files.")
    args = vars(ap.parse_args())
    # main
    # helper function to wrap fasta sequence to 60 characters per line
    def wrap_fasta_seq(seq):
        return '\n'.join([seq[i:i+60] for i in range(0, len(seq), 60)])
    # import the txt file with headers you want to extract the sequence from the input fasta
    headers = (line.rstrip() for line in open(args['ids']))
    # import fasta file
    features = Fasta(args['input'])
    # choose program
    program = args['program']
    match program:
        # extract a multi-fasta file
        case 'extract sequences from a multi-fasta file':           
            # iterate input headers to extract sequences and export as multi-fasta
            with  open(args['output'], 'w') as f:
                for header in headers:
                    f.write(f'>{str(features[str(header)].long_name).rstrip()}\n{wrap_fasta_seq(features[str(header)][:].seq)}\n')
        case 'extract many single-fasta files':
            # extract many single fasta files
            os.chdir(args['directory'])
            for header in headers:
                with open(str(header)+".fasta", 'w') as f:
                    f.write(f'>{str(features[str(header)].long_name).rstrip()}\n{wrap_fasta_seq(features[str(header)][:].seq)}\n')
        case 'remove sequences from a multi-fasta file':
            # remove ids
            keyslist = (features.keys())
            final_keys = (set(keyslist).difference(headers))
            # export to 1 multi-fasta
            with  open(args['output'], 'w') as f:
                for key in final_keys:
                    f.write(f'>{str(features[str(key)].long_name).rstrip()}\n{wrap_fasta_seq(features[str(key)][:].seq)}\n')                
        case 'remove sequences and export to many single-fasta files':
            # remove ids
            keyslist = (features.keys())
            final_keys = (set(keyslist).difference(headers))
            # extract many sigle-fasta files
            os.chdir(args['directory'])
            for key in final_keys:
                with open(str(key)+".fasta", 'w') as f:
                    f.write(f'>{str(features[str(key)].long_name).rstrip()}\n{wrap_fasta_seq(features[str(key)][:].seq)}\n')
                # extract a multi-fasta file
        case 'extract sequences to a 2-column txt file':           
            # iterate input headers to extract sequences and export as 3 column txt
            with  open(args['output'], 'w') as f:
                f.write(f'{"id"}\t{"seq"}\n')
                for header in headers:
                    f.write(f'{str(header)}\t{features[str(header)][:].seq}\n')
        case 'remove sequences and export the rest to a 2-column txt file':
            # remove ids           
            keyslist = (features.keys())
            final_keys = (set(keyslist).difference(headers))            
            # iterate input headers to extract sequences and export as 3 column txt
            with  open(args['output'], 'w') as f:
                f.write(f'{"id"}\t{"seq"}\n')
                for key in final_keys:
                    f.write(f'{str(key)}\t{features[str(key)][:].seq}\n')
        case 'extract sequences to a 3-column txt file':           
            # iterate input headers to extract sequences and export as 3 column txt
            with  open(args['output'], 'w') as f:
                f.write(f'{"id"}\t{"description"}\t{"seq"}\n')
                for header in headers:
                    try:
                        f.write(f'{str(header)}\t{str(str(features[str(header)].long_name).rstrip()).split(" ",1)[1]}\t{features[str(header)][:].seq}\n')
                    except IndexError:
                        f.write(f'{str(header)}\t{""}\t{features[str(header)][:].seq}\n')
        case 'remove sequences and export the rest to a 3-column txt file':
            # remove ids           
            keyslist = (features.keys())
            final_keys = (set(keyslist).difference(headers))            
            # iterate input headers to extract sequences and export as 3 column txt
            with  open(args['output'], 'w') as f:
                f.write(f'{"id"}\t{"description"}\t{"seq"}\n')
                for key in final_keys:
                    try:
                        f.write(f'{str(key)}\t{str(str(features[str(key)].long_name).rstrip()).split(" ",1)[1]}\t{features[str(key)][:].seq}\n')
                    except IndexError:
                        f.write(f'{str(key)}\t{""}\t{features[str(key)][:].seq}\n')
                    
if __name__ == '__main__':
    main()
            