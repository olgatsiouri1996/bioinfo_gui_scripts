# python3
from gooey import *
import os
# input parameters
@Gooey(required_cols=2, program_name= 'merge singlefasta to multifasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-dir", "--directory", required=True, widget='DirChooser', help="input directory with fasta files")
    ap.add_argument("-fa", "--multifasta", required=True, widget='FileSaver', help="output multi-fasta file")
    args = vars(ap.parse_args())
# main
    DIR = args['directory']
    oh = open( args['multifasta'], 'w')
    for f in os.listdir(DIR):
        fh = open(os.path.join(DIR, f))
        for line in fh:
            oh.write(line)
        fh.close()
    oh.close()

if __name__ == '__main__':
        main()


