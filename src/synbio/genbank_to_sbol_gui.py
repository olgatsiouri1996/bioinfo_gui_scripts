# python3
from gooey import *
import synbiopython
import synbiopython.genbabel as stdgen
# imput parameters
@Gooey(required_cols=2, program_name='genbank to sbol', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-gb", "--genbank", required=True, widget='FileChooser', help="input genbank file")
    ap.add_argument("-sbol", "--sbol", required=True, widget='FileSaver', help="outpul sbol file")
    args = vars(ap.parse_args())
# main
    stdconv = stdgen.GenSBOLconv()
    uri_Prefix_igb = 'http://synbiohub.org/public/igem'
    stdconv.run_sbolvalidator(args['genbank'],'SBOL2', uri_Prefix_igb, outputfile = args['sbol'])

if __name__ == '__main__':
    main()
