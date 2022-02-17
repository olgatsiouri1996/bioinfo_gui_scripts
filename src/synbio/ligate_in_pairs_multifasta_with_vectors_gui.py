# python3
import os
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
import  pandas as pd
# imput parameters
@Gooey(required_cols=3, program_name='ligate in pairs a multi-fasta files with vectors', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="ligate in pairs vectors in genbank format with annotations, with inserts in a multi-fasta file")
    ap.add_argument("-txt", "--txt", required=False, widget='FileChooser', help="input 1-column tab-seperated txt file with genbank filename in each row(with extensions .fa, .fasta )")
    ap.add_argument("-fasta", "--fasta", required=False, widget='FileChooser', help="multi-fasta file to inport")
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser',  help="directory to search for vectors in genbank format ")
    args = vars(ap.parse_args())
# main
# set working directory
    os.chdir(args['directory'])
# inport txt file and convert each column to list
    df_txt = pd.read_csv(args['txt'], header=None, sep="\t")
    gb_list = df_txt.iloc[:,0].values.tolist()
# create lists
    gb_seqs = []
    fasta_seqs = []
    gb_features = []
    fasta_list = []
# linear vectors
# import each genbank file from the list
    for i in gb_list:
        plasmid = SeqIO.read(i, "genbank")
        gb_seqs.append(str(plasmid.seq))
        gb_features.append(plasmid.features)
# DNA insert
# import each fasta file from the list
    for record in SeqIO.parse(args['fasta'], "fasta"):
        fasta_seqs.append(str(record.seq))
        fasta_list.append(record.id)
# iterate all below lists in pairs
    for (a,b,c,d,e) in zip(gb_list,fasta_list,gb_seqs,gb_features,fasta_seqs):
        # merge
        seqad = str(c) + str(e)
        # add this record to the list
        ligated = SeqRecord(Seq(seqad),id='_'.join([str(a).split(".")[0], str(b)]),description="",annotations={"molecule_type":"DNA","topology":"circular"})
        ligated.features = d
        # export to genbank
        SeqIO.write(ligated,"".join([str(a).split(".")[0],"_",str(b),".gb"]), "genbank")

if __name__ == '__main__':
    main()