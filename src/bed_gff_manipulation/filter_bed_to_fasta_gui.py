from gooey import *
from pyfaidx import Fasta
import pandas as pd
import warnings
import textwrap

# input parameters
@Gooey(
    program_name='filter bed to fasta or txt',
    header_bg_color='#DCDCDC',
    terminal_font_color='#DCDCDC',
    terminal_panel_color='#DCDCDC',
    default_size=(610,580)
)
def main():
    ap = GooeyParser()
    ap.add_argument("-bed", "--bed", required=True, widget='FileChooser',
                    help="input bed file with 0-based start coordinates and columns 1-4, 6 and 8 as: chromosome/scaffold, start coordinate, end coordinate, name, strand and feature type")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
    ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output file (choose .fasta or .txt)")
    ap.add_argument("-fea", "--feature", required=False, default='gene', choices=['gene','CDS','exon','intron','promoter'], type=str, help="specify the feature to collect sequences for")
    ap.add_argument("-strandedness", "--strandedness", required=False, default='no', choices=['yes', 'no'], type=str, help="reverse complement sequences on '-' strand")
    ap.add_argument("-program", "--program", required=False, default='filter', choices=['filter', 'no_filter'], type=str, help="choose whether to filter input bed file by feature")
    args = vars(ap.parse_args())

    # ignore warnings
    warnings.filterwarnings('ignore')

    # Define the process_row function
    def process_row(row, fasta, output_format, strandedness):
        chrom, start, end, ids, strand, feature = row

        if args['program'] == 'filter' and feature != args['feature']:
            return None  # Skip this row and continue processing other rows

        if str(strand) == "-" and strandedness == 'yes':
            header = f"{ids} {chrom}:{int(start) + 1}-{end} ({strand}) reverse complement".replace('\r', '')
            seq = fasta[str(chrom)][int(start):end].reverse.complement.seq

        if str(strand) == "-" and strandedness == 'no' or str(strand) == "+" and strandedness == 'no' or str(strand) == "+" and strandedness == 'yes':
            header = f"{ids} {chrom}:{int(start) + 1}-{end} ({strand})".replace('\r', '')    
            seq = fasta[str(chrom)][int(start):end].seq

        if output_format == 'fasta':
            header = f">{header}"
            return f"{header}\n{textwrap.fill(seq, width=60)}\n"
        elif output_format == 'txt':
            return f"{ids}\t{chrom}\t{int(start) + 1}\t{end}\t{strand}\t{seq}"

    # Read the BED file into a Pandas DataFrame, selecting columns 0 to 5 and column 7
    df = pd.read_csv(args['bed'], sep="\t", header=None, usecols=[0, 1, 2, 3, 5, 7])

    # Import the FASTA file
    fasta = Fasta(args['input'])

    # Determine the output format based on the chosen output file extension
    output_format = 'fasta' if args['output'].endswith('.fasta') else 'txt'

    # Process the DataFrame and write results to the output file
    with open(args['output'], 'w') as output_file:
        if output_format == 'fasta':
            for _, row in df.iterrows():
                result = process_row(row, fasta, output_format, args['strandedness'])
                if result is not None:
                    output_file.write(result)
        elif output_format == 'txt':
            header = "id\tchrom\tstart\tend\tstrand\tsequence"
            output_file.write(header + '\n')
            for _, row in df.iterrows():
                result = process_row(row, fasta, output_format, args['strandedness'])
                if result is not None:
                    output_file.write(result + '\n')

if __name__ == '__main__':
    main()
