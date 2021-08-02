# python3
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
# imput parameters
@Gooey(required_cols=3, program_name='digital ligation', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="ligate linear vector with DNA insert")
    ap.add_argument("-vr", "--vector", required=True, widget='FileChooser', help="linear vector/s (single or multi fasta file)")
    ap.add_argument("-in", "--insert", required=True, widget='FileChooser', help="sequence/s to insert in the vector(single or multi fasta file)")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output genbank file with circular sequence/s")
    args = vars(ap.parse_args())
# main 
# linear vector
    output_handle = open(args['output'], "w")
    records = [] # setup an empty list
    for plasmid in SeqIO.parse(args['vector'], "fasta"):
        x = str(plasmid.seq)
        # DNA insert
        for record in SeqIO.parse(args['insert'], "fasta"):
            y = str(record.seq)
        # merge
            seqad = y + x
            # add this record to the list
            records.append(SeqRecord(Seq(seqad,generic_dna),id='_'.join([record.id,plasmid.id]),description="",annotations={"topology":"circular"}))
    # export to fasta
    SeqIO.write(records,output_handle, "genbank")
    output_handle.close()

if __name__ == '__main__':
    main()
