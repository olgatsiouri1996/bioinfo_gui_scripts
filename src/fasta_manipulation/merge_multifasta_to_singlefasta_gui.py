# python3
from gooey import *
from pyfaidx import Fasta
import textwrap
# imput parameters
@Gooey(required_cols=3, program_name="merge multi-fasta to single-fasta",default_size=(800, 530), header_bg_color= "#DCDCDC", terminal_font_color= "#DCDCDC", terminal_panel_color= "#DCDCDC")
def main():
    ap = GooeyParser()
    ap.add_argument("-mfa", "--multi fasta", required=True, widget="FileChooser", help="input multi-fasta file to merge its sequences")
    ap.add_argument("-type", "--sequence type", required=False, type=str, default="none",choices=["none","DNA", "aa"], help="sequence type of unknown DNA/AAs to add between the merged fasta sequences. Default: no sequence to add")
    ap.add_argument("-num", "--sequence number", required=False, type=int, help="number of nucleotides or aminoacids to add between the merged fasta sequences. If sequence type is none then ignore this option")
    ap.add_argument("-id", "--fasta header", required=True, type=str, help="fasta header of the output file")
    ap.add_argument("-wi", "--fasta width", required=False, type=int, default=80, choices=[60,70,80,100,120], help="number of fasta sequence characters per line")
    ap.add_argument("-sfa", "--single fasta", required=True, widget="FileSaver", help="output single-fasta file")
    args = vars(ap.parse_args())
    # main
    # index fasta file
    features = Fasta(args["multi fasta"],as_raw=True)
    # store all sequences to a list
    sequences = [features[key][:] for key in features.keys()]
    # select sequence type
    if args["sequence type"]=="none":
        spacer = ""
    elif args["sequence type"]=="DNA":
        spacer = "N"*args["sequence number"]
    else:
        spacer = "X"*args["sequence number"]
    # merge all sequences at 1 
    merged_seqs = spacer.join(sequences)
    # export to fasta
    with open(args["single fasta"], "w") as f:
        f.write(f'>{args["fasta header"]}\n{textwrap.fill(merged_seqs, width=args["fasta width"])}')

if __name__ == "__main__":
    main()

