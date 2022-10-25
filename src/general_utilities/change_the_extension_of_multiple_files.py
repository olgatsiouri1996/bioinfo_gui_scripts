# python3
from gooey import *
import os
# imput parameters
@Gooey(required_cols=1, program_name='change the extension of multiple files', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser', help="directory with files to change extension")
    ap.add_argument("-ext", "--extension", required=False, type=str, default='.fasta', help="extension to change into")
    args = vars(ap.parse_args())
    # main
    # change the file extension
    for filename in sorted(os.listdir(os.chdir(args['directory']))):
        base = os.path.splitext(filename)[0]
        os.rename(filename, base + args['extension'])

if __name__ == '__main__':
    main()
