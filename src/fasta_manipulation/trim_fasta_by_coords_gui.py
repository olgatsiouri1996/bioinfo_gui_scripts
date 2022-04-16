# python3
from gooey import *
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import  pandas as pd
# input parameters
@Gooey(required_cols=1, program_name='trim a multi-fasta file or multiple single-fasta files based on a txt file with start and end coordinates', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-in", "--input", required=False, widget='FileChooser', help="input fasta file")
    ap.add_argument("-coords", "--coordinates", required=True, widget='FileChooser', help="input 4-column tab-seperated txt file with id, start, end positions and strand(+, -) respectively in each row")
    ap.add_argument("-type", "--type", required=False,default=1, type=int, help="type of fasta to import 1) 1 multi-fasta file 2)  many single-fasta files")
    ap.add_argument("-fastadir", "--fastadir", required=False, type=str, widget='DirChooser', help="directory to search for fasta files")
    ap.add_argument("-outdir", "--outdir", required=False, type=str, widget='DirChooser', help="directory to save output fasta files")
    ap.add_argument("-out", "--output", required=False, widget='FileSaver', help="output multi-fasta file")
    args = vars(ap.parse_args())
    # main
    # inport txt file and convert each column to list
    df_txt = pd.read_csv(args['coordinates'], header=None, sep="\t")
    ids = df_txt.iloc[:,0].values.tolist()
    seq_start = df_txt.iloc[:,1].values.tolist()
    seq_start[:] = [i - 1 for i in seq_start]
    seq_end = df_txt.iloc[:,2].values.tolist()
    seq_strand = df_txt.iloc[:,3].values.tolist()
    # create function to choose strand from
    def choose_strand(seq,strand):
        if strand == "+":
            rec=str(seq)
        else:
            rec=str(Seq(seq).complement())

        return rec
    # create function to rename strand
    def rename_strand(strand):
        if strand == "+":
            return "plus"
        else:
            return "minus"
    # setup empty lists
    records = []
    trimmed_records = []
    # choose fasta type to import
    if args['type'] == 1:    
        # iterate for each record
        for i in ids:
            for record in SeqIO.parse(args['input'], "fasta"):
                if i == record.id:
                    records.append(record)
        # iterate all below lists in pairs
        for (a, b, c, d) in zip(records, seq_start, seq_end, seq_strand):
            trimmed_records.append(SeqRecord(Seq(choose_strand(str(a.seq),str(d))[int(b):int(c)]), id='_'.join([str(a.id),str(b + 1),str(c),rename_strand(str(d))]), description=""))
        # export to fasta
        SeqIO.write(trimmed_records, args['output'], "fasta")
    else:
        # import each fasta file from the working directory
        os.chdir(args['fastadir'])
        for i in ids:
            # read each file
            record = SeqIO.read(''.join([i,".fasta"]), "fasta")
            # add this record to the lists
            records.append(record)
        # select directory to save the output files
        os.chdir(args['outdir'])
        # iterate all below lists in pairs
        for (a, b, c, d) in zip(records, seq_start, seq_end, seq_strand):
            trimmed_record = SeqRecord(Seq(choose_strand(str(a.seq),str(d))[int(b):int(c)]), id='_'.join([str(a.id),str(b + 1),str(c),rename_strand(str(d))]), description="")
            # export to fasta
            SeqIO.write(trimmed_record, "".join([trimmed_record.id,".fasta"]), "fasta")

if __name__ == '__main__':
    main()