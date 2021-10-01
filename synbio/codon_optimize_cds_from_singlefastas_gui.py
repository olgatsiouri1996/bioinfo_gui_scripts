# python3
import os
from gooey import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dnachisel import *
# imput parameters
@Gooey(required_cols=2, program_name='codon optimize cds from single-fastas', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-org","--organism", required=True, help="organism to input(use either the names of the genomes avaliable on dnachisel or use the taxid of the organisms in http://www.kazusa.or.jp/codon/)")
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser', help="directory to search for fasta files")
    args = vars(ap.parse_args())
    # main
    # import each fasta file from the working directory
    for filename in sorted(os.listdir(os.chdir(args['directory']))):
        if filename.endswith(".fa") or filename.endswith(".fasta"):
            record = SeqIO.read(filename, "fasta")
            problem = DnaOptimizationProblem(sequence=str(record.seq),
            constraints=[EnforceTranslation()],
            objectives=[CodonOptimize(species= args['organism'])])
            problem.optimize()
            # create seqRecord with optimized sequence
            optimized_seq=SeqRecord(Seq(problem.sequence),id="".join([record.id,"_",args['organism']]),description="")
            # export to fasta
            SeqIO.write(optimized_seq, "".join([filename.split(".")[0],"_",args['organism'],".fasta"]), "fasta")

if __name__ == '__main__':
    main()