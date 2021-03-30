# python3
from gooey import *
from pymut import *
# input parameters
@Gooey(required_cols=4, program_name='mutate aminoacid in structure using a rotamer library', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input pdb file")
    ap.add_argument("-lib", "--library", required=True, widget='FileChooser', help="input rotamer library file")
    ap.add_argument("-chain", "--chain", required=False, default= 'A', type=str,  help="input genbank file")
    ap.add_argument("-res", "--residue", required=True, type=int,  help="number of residue to mutate")
    ap.add_argument("-mut", "--mutation", required=True,  help="resulted aminoacid after the mutation")
    ap.add_argument("-type", "--type", required=False, default= 'first',  type=str,  help="mutation type(first: Select the most likely rotamer based on probability in the library. random: Select a rotamer randomly based on the probability in the library. best: Select the best rotamer based on VdW energy)")
    ap.add_argument("-out","--output", required=True, widget='FileSaver', help="output pdb file")
    args = vars(ap.parse_args())
# main
    rotamer_lib = load_rotamers(args['library'].format(DATA_DIR))
    pdb_obj = PDB(args['input'])
    pdb_obj.mutate(chain=args['chain'], mutate_to=args['mutation'], res_num=args['residue'], mutation_type=args['type'], rotamer_lib=rotamer_lib)
    pdb_obj.dump(args['output'])

if __name__ == '__main__':
    main()

