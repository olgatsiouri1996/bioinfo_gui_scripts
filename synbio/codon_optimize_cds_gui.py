# python3
from gooey import *
from Bio import SeqIO
from dnachisel import *
# imput parameters
@Gooey(required_cols=3, program_name='codon optimize cds', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser()
    ap.add_argument("-fa", "--fasta", required=True, widget='FileChooser', help="input fasta file")
    ap.add_argument("-org","--organism", required=True, help="organism to input(use either the names of the genomes avaliable on dnachisel or use the taxid of the organisms that exist in http://www.kazusa.or.jp/codon/)")
    ap.add_argument("-opt","--optimized", required=True, widget='FileSaver', help="optimized genbank file")
    args = vars(ap.parse_args())
# main
    for record in SeqIO.parse(args['fasta'], "fasta"):
        problem = DnaOptimizationProblem(sequence=str(record.seq),
        objectives=[CodonOptimize(species= args['organism'])])
        problem.optimize()
    final_record = problem.to_record(args['optimized'])

if __name__ == '__main__':
    main()

