# python3
from gooey import *
from biopandas.pdb import PandasPdb
import matplotlib
import matplotlib.pyplot as plt

# input parameters
@Gooey(required_cols=3, program_name='plot B factor by chain', header_bg_color='#DCDCDC', terminal_font_color='#DCDCDC', terminal_panel_color='#DCDCDC',default_size=(610,750))
def main():
    ap = GooeyParser(description="creates a line plot or histogram with B factor values from a selected chain")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input pdb file")
    ap.add_argument("-chain", "--chain", required=False,default='A', help="chain from pdb file to select")
    ap.add_argument("-col", "--col", required=True, help="colour of line/histogram in plot (you can use hex colour codes)")
    ap.add_argument("-plot", "--plot", required=True, widget='FileSaver', help="export plot to file")
    ap.add_argument("-dpi", "--dpi", default=300, required=False, help="dpi of exported plot")
    ap.add_argument("-type", "--type", required=False, choices=['pdf', 'svg'], default='pdf', help="type of plot file")
    ap.add_argument("-start", "--start", required=False, type=int, default=1, help="start residue position in the chain")
    ap.add_argument("-end", "--end", required=False, type=int, help="end residue position in the chain (optional)")
    ap.add_argument("-plot_type", "--plot_type", required=False, choices=['line', 'histogram'], default='line', help="type of plot (line or histogram)")
    ap.add_argument("-width", "--width", required=False, default=8, type=float, help="width of the output figure in inches")
    ap.add_argument("-height", "--height", required=False, default=6, type=float, help="height of the output figure in inches")
    args = vars(ap.parse_args())

    # main
    ppdb = PandasPdb()
    ppdb.read_pdb(args['input'])

    # Filter by chain and residue positions
    chain_id = args['chain']
    start_residue = args['start']
    end_residue = args['end']

    if end_residue is None:
        # If end residue is not provided, set it to the residue number of the last residue in the chain
        chain_residue_numbers = ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id'] == chain_id]['residue_number']
        end_residue = max(chain_residue_numbers)

    ppdb.df['ATOM'] = ppdb.df['ATOM'][(ppdb.df['ATOM']['chain_id'] == chain_id) &
                                      (ppdb.df['ATOM']['residue_number'] >= start_residue) &
                                      (ppdb.df['ATOM']['residue_number'] <= end_residue)]

    # Check if any data exists after filtering
    if ppdb.df['ATOM'].empty:
        print("No data available for the specified chain and residue positions.")
        return

    # Set figure size
    plt.figure(figsize=(args['width'], args['height']))

    if args['plot_type'] == 'line':
        # Plot B-factor values as a line plot
        ppdb.df['ATOM']['b_factor'].plot(kind='line', color=args['col'])
    elif args['plot_type'] == 'histogram':
        # Plot B-factor values as a histogram
        ppdb.df['ATOM']['b_factor'].plot(kind='hist', color=args['col'], bins=20, alpha=0.7)

    plt.title('B-Factors Along the Amino Acid Chain')
    plt.xlabel('Residue Number')
    plt.ylabel('B-factor in $A^2$')

    # Export figure
    plt.savefig(args['plot'], dpi=int(args['dpi']), format=args['type'])

if __name__ == '__main__':
    main()
