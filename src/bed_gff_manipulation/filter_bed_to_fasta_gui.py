from gooey import *
from pyfaidx import Fasta
import pandas as pd
import warnings
import textwrap

# input parameters
@Gooey(
    required_cols=3,
    program_name='filter bed to fasta',
    header_bg_color='#DCDCDC',
    terminal_font_color='#DCDCDC',
    terminal_panel_color='#DCDCDC',
    default_size=(610,580)
)
def main():
    ap = GooeyParser()
    ap.add_argument("-bed", "--bed", required=True, widget='FileChooser',
                    help="input bed file(made with bedops's gff2bed function, every feature in the '.gff' or '.gff3' file should have an 'ID' tag in the 'attributes' column)")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
    ap.add_argument("-fea", "--feature", required=False, default='gene', choices=['gene','CDS','exon','intron','promoter'], type=str, help="specify the feature to collect sequences for")
    ap.add_argument("-pro", "--program", required=False, default='filter the bed file by feature before retrieving sequences', type=str, choices=['filter the bed file by feature before retrieving sequences','do not filter'], help="program to choose")
    args = vars(ap.parse_args())

    # ignore warnings
    warnings.filterwarnings('ignore')

    # Define the process_row function
    def process_row(row, output_file, features):
        chrom, start, end, ids, strand = row

        if str(strand) == "+":
            header = f">{ids} {chrom}:{int(start) + 1}-{end}".replace('\r', '')
            seq = features[str(chrom)][int(start):end].seq
        else:
            header = f">{ids} {chrom}:{int(start) + 1}-{end} reverse complement".replace('\r', '')
            seq = features[str(chrom)][int(start):end].reverse.complement.seq

        output_file.write(header + '\n')
        output_file.write(textwrap.fill(seq, width=60) + '\n')

    # import bed with no headers specified
    df = pd.read_csv(args['bed'], sep="\t", header=None)

    # choose program
    if args['program'] == 'filter the bed file by feature before retrieving sequences':
        # select the rows containing the feature
        bool2 = df.iloc[:, 7].str.contains(args['feature'])
        df = df[bool2]

    # import fasta file
    features = Fasta(args['input'])

    # Process the DataFrame and write results to the output file
    with open(args['output'], 'w') as output_file:
        for _, row in df.iloc[:, [0, 1, 2, 3, 5]].iterrows():
            process_row(row, output_file, features)

if __name__ == '__main__':
    main()
