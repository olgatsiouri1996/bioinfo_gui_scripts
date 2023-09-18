# python 3
from gooey import *
from pyfaidx import Fasta
import textwrap
# input parameters
@Gooey(required_cols=2, program_name='deduplicate fasta by id',default_size=(610,440), header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="keep only the 1st of the duplicated sequences(that have the same identifier and sequence)")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input multi-fasta file")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output multi-fasta or txt file")
    ap.add_argument("-type", "--output type", required=False,default='multi-fasta',choices=['multi-fasta','3-column txt(1d,description,seq)','1-column txt(id)'], help="output type")
    args = vars(ap.parse_args())
    # create fasta index
    features = Fasta(args['input'], duplicate_action="first")
    # select output type
    if args['output type']=="multi-fasta":
        # create fasta format function
        def make_fasta_format(key):
            header = str(features[key].long_name).rstrip()
            seq= textwrap.fill(features[key][:].seq,width=60)
            return f">{header}\n{seq}"
        # select sequences
        fasta_records = map(make_fasta_format,features.keys())
        # export to fasta
        with open(args['output'], 'w') as f:
            f.write('\n'.join(fasta_records))
        
    elif args['output type']== "3-column txt(1d,description,seq)":
        # make 3 column txt format
        def make_txt_format(key):
            header = str(features[key].long_name).rstrip()
            try:
                description = header.split(' ',1)[1]
            except IndexError:
                description = ''
            seq= features[key][:].seq
            return f"{key}\t{description}\t{seq}"
        # select sequences
        fasta_records = map(make_txt_format,features.keys())
        # export to 3-column txt
        with open(args['output'], 'w') as f:
            f.write(f'{"id"}\t{"description"}\t{"seq"}\n')
            f.write('\n'.join(fasta_records))

    else:
        # export to 1-column ids txt
        with open(args['output'], 'w') as f:
            f.write('\n'.join(features.keys()))

if __name__ == '__main__':
    main()