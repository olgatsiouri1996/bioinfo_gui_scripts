# python3
from gooey import *
import os
# imput parameters
@Gooey(required_cols=0, program_name='change file extension', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=False, widget='FileChooser', help="input file to change the extension")
    ap.add_argument("-ext", "--extension", required=False, type=str, default='.fasta', help="extension to change into")
    ap.add_argument("-num", "--number", required=False, type=str, default='one', widget='Dropdown', choices=['one','many'], help="number of  files to change extensions")
    ap.add_argument("-dir", "--directory", required=False, type=str, widget='DirChooser', help="directory with files to change extension")
    args = vars(ap.parse_args())
    # main
    # function to change the file extension
    def change_extension(fi):
        base = os.path.splitext(fi)[0]
        os.rename(fi, base + args['extension'])
    # choose number of input files
    if args['number'] == 'one':
        change_extension(args['input'])
    else:
        for filename in sorted(os.listdir(os.chdir(args['directory']))):
            change_extension(filename)

if __name__ == '__main__':
    main()
