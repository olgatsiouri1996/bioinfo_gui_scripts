# python3
import os
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqFeature
# imput parameters
@Gooey(required_cols=2, program_name='ligate single-fastas with vector', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="ligate vector with inserts in single-fasta files")
    ap.add_argument("-vr", "--vector", required=True,  help="vector in genbank format")
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser', help="directory to search for fasta files")
    args = vars(ap.parse_args())
    # main 
    # linear vector
    plasmid = SeqIO.read(args['vector'], "genbank")
    x = str(plasmid.seq)
    # DNA insert
    # import each fasta file from the working directory
    for filename in sorted(os.listdir(os.chdir(args['directory']))):
        if filename.endswith(".fa") or filename.endswith(".fasta"):
            record = SeqIO.read(filename, "fasta")
            y = str(record.seq)
            # merge
            seqad = x + y
            # add this record to the list
            ligated = SeqRecord(Seq(seqad,generic_dna),id='_'.join([record.id,args['vector']]),description="",annotations={"topology":"circular"})
            ligated.features = plasmid.features
            # export to genbank
            SeqIO.write(ligated,"".join([filename.split(".")[0],"_",args['vector'].split(".")[0],".gb"]), "genbank")

if __name__ == '__main__':
    main()
