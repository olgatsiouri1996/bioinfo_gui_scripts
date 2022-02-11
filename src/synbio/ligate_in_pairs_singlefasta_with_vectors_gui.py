# python3
import os
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
import pandas as pd
# imput parameters
@Gooey(required_cols=2, program_name='ligate in pairs single-fasta files with vectors', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="ligate in pairs vectors in genbank format with annotations, with inserts in single-fasta files")
    ap.add_argument("-txt", "--txt", required=True, type=str, widget='FileChooser', help="input tab-seperated txt file with fasta and genbank filenames in each row(with extensions .gb, .gbk, .fa, .fasta and column names genbank and fasta respectively)")
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser',  help="directory to search for vectors in genbank format and inserts in single-fasta files")
    args = vars(ap.parse_args())
# main
# set working directory
    os.chdir(args['directory'])
# inport txt file and convert each column to list
    df_txt = pd.read_csv(args['txt'], sep="\t")
    gb_list = df_txt.iloc[:,0].values.tolist()
    fasta_list = df_txt.iloc[:,1].values.tolist()
# create lists
    gb_seqs = []
    fasta_seqs = []
    gb_features = []
# linear vectors
# import each genbank file from the list
    for i in gb_list:
        plasmid = SeqIO.read(i, "genbank")
        gb_seqs.append(str(plasmid.seq))
        gb_features.append(plasmid.features)
# DNA insert
# import each fasta file from the list
    for i in fasta_list:
        record = SeqIO.read(i, "fasta")
        fasta_seqs.append(str(record.seq))
# iterate all below lists in pairs
    for (a,b,c,d,e) in zip(gb_list,fasta_list,gb_seqs,gb_features,fasta_seqs):
    # merge
        seqad = str(c) + str(e)
    # add this record to the list
        ligated = SeqRecord(Seq(seqad),id='_'.join([str(a).split(".")[0], str(b).split(".")[0]]),description="",annotations={"molecule_type":"DNA","topology":"circular"})
        ligated.features = d
    # export to genbank
        SeqIO.write(ligated,"".join([str(a).split(".")[0],"_",str(b).split(".")[0],".gb"]), "genbank")

if __name__ == '__main__':
    main()