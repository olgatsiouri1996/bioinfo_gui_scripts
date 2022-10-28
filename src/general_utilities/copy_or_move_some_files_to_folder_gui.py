# python3
from gooey import *
import os
import shutil
# imput parameters
@Gooey(required_cols=3, program_name='copy or move some files to folder', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input 1-column txt file with filenames to copy. if number is set to many then each filename should be followed by an extension e.g. file.txt")
    ap.add_argument("-ext", "--extension", required=False, type=str, default='.fasta', help="extension of the input files")
    ap.add_argument("-act", "--action", required=False, type=str, default='copy', choices=['copy','move'], help="copy or move selected files")
    ap.add_argument("-num", "--number", required=False, type=str, default='one', choices=['one','many'], help="number of extensions selected files have") 
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser', help="directory with files to copy or move")
    ap.add_argument("-out", "--output", required=True, type=str, widget='DirChooser', help="directory to save copied or moved files")
    args = vars(ap.parse_args())
    # main
    # add filenames to copy
    with open(args['input'], 'r') as f:
        input_files = f.readlines()
    # select number of extensions
    if args['number']=='one':
        input_files = [x.strip() + args['extension'] for x in input_files]
    else:
        input_files = [x.strip() for x in input_files]
    # go to the directory where the files you want to copy are located
    os.chdir(args['directory'])
    # choose to copy or move files
    if args['action'] == 'copy':
        # copy each file
        for file in input_files:
            shutil.copy(file,args['output'])
    else:
        # copy each file
        for file in input_files:
            shutil.move(file,args['output']) 

if __name__ == '__main__':
    main()
