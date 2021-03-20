# python3
from gooey import *
import synbiopython
import synbiopython.genbabel as stdgen
# imput parameters
@Gooey(required_cols=3, program_name='circut plotter', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=True, help="input circut components")
    ap.add_argument("-reg", "--regulations", required=True, help="regulations of circut components(e.g activation)")
    ap.add_argument("-plot", "--plot", required=True, widget='FileSaver', help="file to save the circut plot")
    args = vars(ap.parse_args())
# main
    simplot = stdgen.SimpleDNAplot()
    Input = args['input']
    Regulations = args['regulations']
    maxdnalength, figure = simplot.plot_circuit(Input, Regulations, args['plot'])

if __name__ == '__main__':
    main()
