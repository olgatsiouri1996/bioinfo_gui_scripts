# python3
from gooey import *
import os
import shutil
# imput parameters
@Gooey(required_cols=3, program_name='copy some files to folder', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input txt file with filenames to copy")
    ap.add_argument("-ext", "--extension", required=False, type=str, default='.fasta', help="extension of the input files")
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser', help="directory with files to copy")
    ap.add_argument("-out", "--output", required=True, type=str, widget='DirChooser', help="directory to save copied files")
    args = vars(ap.parse_args())
    # main
    # add filenames to copy
    with open(args['input'], 'r') as f:
        input_files = f.readlines()
    input_files = [x.strip() + args['extension'] for x in input_files]
    # go to the directory where the files you want to copy are located
    os.chdir(args['directory'])
    # copy each file
    for file in input_files:
        shutil.copy(file,args['output'])

if __name__ == '__main__':
    main()
