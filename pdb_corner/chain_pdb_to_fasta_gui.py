# python3
from gooey import *
from Bio.PDB import *
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# input parameters
@Gooey(required_cols=2, program_name='pdb chain to fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="subsets a pdb file by selecting the model and chain and convert the output to fasta format")
    ap.add_argument("-pdb", "--pdb", required=True, widget='FileChooser', help="input pdb file")
    ap.add_argument("-model", "--model",default=0, required=False, help="model from pdb file to select(integer). Default is 0(1 model only)")
    ap.add_argument("-chain", "--chain", required=True, help="chain from pdb file to select")
    args = vars(ap.parse_args())
# main
# select chain
    parser = PDBParser()
    s = parser.get_structure("name", args['pdb'])
    fill = s[int(args['model'])][args['chain']]
# retrieve the pdb id of the input file
    filename = os.path.split(args['pdb'])[1]
    pdb_id = filename.split(".")[0]
# export to fasta
    ppb = PPBuilder()
    for pp in ppb.build_peptides(fill):
        record = SeqRecord(Seq(str(pp.get_sequence())),id="".join([str(pdb_id),"_",str(args['chain'])]),description="")
        SeqIO.write(record, "".join([str(pdb_id),"_",str(args['chain']),".fasta"]), "fasta")
     
if __name__ == '__main__':
    main()
