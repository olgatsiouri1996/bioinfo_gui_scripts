# python3
from gooey import *
import os
import glob
import shutil
# imput parameters
@Gooey(required_cols=2, program_name='copy or move some files to folder', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser', help="directory with files to copy or move")
    ap.add_argument("-num", "--number", required=False, type=str, default='one', choices=['one','many'], help="number of extensions selected files have") 
    ap.add_argument("-in", "--input", required=False, widget='FileChooser', help="input 1-column txt file with filenames to copy. Each filename should be followed by an extension e.g. file.txt")
    ap.add_argument("-pat", "--pattern", required=False, type=str, default='.fasta', help="prefix or suffix of the input files")
    ap.add_argument("-type", "--type", required=False, type=str, default='suffix', choices=['suffix','prefix'], help="select the files based on suffix or prefix")
    ap.add_argument("-act", "--action", required=False, type=str, default='copy', choices=['copy','move'], help="copy or move selected files")
    ap.add_argument("-out", "--output", required=True, type=str, widget='DirChooser', help="directory to save copied or moved files")
    args = vars(ap.parse_args())
    # main
    # add filenames to copy
    with open(args['input'], 'r') as f:
        input_files = f.readlines()
    # select number of extensions
    if args['number']=='one':
        if args['type']=='suffix':
            input_files = glob.glob(os.path.join(args['directory'], '*'+args['pattern']))
        else:
            input_files = glob.glob(os.path.join(args['directory'], args['pattern']+'*'))
    else:
        input_files = (x.rstrip() for x in input_files)
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
