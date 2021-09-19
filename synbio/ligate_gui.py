# python3
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqFeature
# imput parameters
@Gooey(required_cols=3, program_name='digital ligation', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="ligate vector with insert")
    ap.add_argument("-vr", "--vector", required=True, widget='FileChooser', help="linear vector in genbank format")
    ap.add_argument("-in", "--insert", required=True, widget='FileChooser', help="sequence to insert in the vector in fasta format")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output genbank file with circular sequence")
    args = vars(ap.parse_args())
# main 
# linear vector
    plasmid = SeqIO.read(args['vector'], "genbank")
    x = str(plasmid.seq)
# DNA insert
    record = SeqIO.read(args['insert'], "fasta")
    y = str(record.seq)
# merge
    seqad = x + y
# add this record to the list
    ligated = SeqRecord(Seq(seqad,generic_dna),id='_'.join([record.id,plasmid.id]),description="",annotations={"topology":"circular"})
    ligated.features = plasmid.features
# export to genbank
    SeqIO.write(ligated,args['output'], "genbank")

if __name__ == '__main__':
    main()
