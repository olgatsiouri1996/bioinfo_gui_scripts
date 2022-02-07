# python3
import os
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
# imput parameters
@Gooey(required_cols=1, program_name='ligate single-fasta files with vectors', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="ligate vectors in genbank format with annotations, with inserts in single-fasta files")
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser',  help="directory to search for vectors in genbank format and inserts in single-fasta files")
    args = vars(ap.parse_args())
    # main 
    # linear vectors
    # import each genbank file from the working directory
    for filename in sorted(os.listdir(os.chdir(args['directory']))):
        if filename.endswith(".gb") or filename.endswith(".gbk"): 
            plasmid = SeqIO.read(filename, "genbank")
            x = str(plasmid.seq)
            # retrieve the name of the input file
            gb_file = filename.split(".")[0]
            # DNA insert
            # import each fasta file from the working directory
            for filename in sorted(os.listdir(os.chdir(args['directory']))):
                if filename.endswith(".fa") or filename.endswith(".fasta"):
                    record = SeqIO.read(filename, "fasta")
                    y = str(record.seq)
                    # merge
                    seqad = x + y
                    # add this record to the list
                    ligated = SeqRecord(Seq(seqad),id='_'.join([record.id,gb_file]),description="",annotations={"molecule_type":"DNA","topology":"circular"})
                    ligated.features = plasmid.features
                    # export to genbank
                    SeqIO.write(ligated,"".join([filename.split(".")[0],"_",gb_file,".gb"]), "genbank")

if __name__ == '__main__':
    main()