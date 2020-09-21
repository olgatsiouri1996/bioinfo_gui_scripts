# python3
from gooey import *
import Bio
from Bio.PDB import *
# input parameters
@Gooey(required_cols=5, program_name='subset pdb', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
  ap = GooeyParser(description="subsets a pdb file by selecting the chain and residues from it")
  ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input pdb file")
  ap.add_argument("-chain", "--chain", required=True, help="chain from pdb file to select")
  ap.add_argument("-start", "--start", required=True, help="amino acid in chain to start writing the pdb file")
  ap.add_argument("-end", "--end", required=True, help="amino acid in chain to end writing the pdb file")
  ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output pdb file")
  args = vars(ap.parse_args())
  #main
  parser = PDBParser()
  s = parser.get_structure("name", args['input'])
  Bio.PDB.Dice.extract(s, args['chain'], int(args['start']), int(args['end']), args['output'])

if __name__ == '__main__':
    main()