# python3
from gooey import *
from biopandas.pdb import PandasPdb
import matplotlib
import matplotlib.pyplot as plt
# input parameters
@Gooey(required_cols=5, program_name='plot B factor by chain', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="creates a line plot with B factor values from a selected chain")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input pdb file")
    ap.add_argument("-chain", "--chain", required=True, help="chain from pdb file to select")
    ap.add_argument("-col", "--col", required=True, help="colour of line in plot(you can use hex colour codes)")
    ap.add_argument("-plot", "--plot", required=False, widget='FileSaver', help="export plot to file")
    ap.add_argument("-dpi", "--dpi", default= 300, required=False, help="dpi of exported plot")
    ap.add_argument("-type", "--type", required=True, help="type of plot file(e.g png etc)")
    args = vars(ap.parse_args())
    # main
    ppdb = PandasPdb()
    ppdb.read_pdb(args['input'])
    ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == args['chain']]
    # export figure
    ppdb.df['ATOM']['b_factor'].plot(kind='line', color= args['col'])
    plt.title('B-Factors Along the Amino Acid Chain')
    plt.xlabel('Residue Number')
    plt.ylabel('B-factor in $A^2$')
    matplotlib.pyplot.savefig(args['plot'], dpi=int(args['dpi']), format=args['type'])

if __name__ == '__main__':
    main()