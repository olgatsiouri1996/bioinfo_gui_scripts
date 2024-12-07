# python3
from gooey import *
from pyfaidx import Fasta
import textwrap
# imput parameters
@Gooey(required_cols=3, program_name="merge multi-fasta to single-fasta",default_size=(800, 530), header_bg_color= "#DCDCDC", terminal_font_color= "#DCDCDC", terminal_panel_color= "#DCDCDC")
def main():
    ap = GooeyParser()
    ap.add_argument("-mfa", "--Multi FASTA", required=True, widget="FileChooser", help="Input multi-fasta file to merge its sequences")
    ap.add_argument("-type", "--Sequence type", required=False, type=str, default="none",choices=["none","DNA", "aa"], help="Sequence type of unknown DNA/AAs to add between the merged fasta sequences. Default: no sequence to add")
    ap.add_argument("-num", "--Sequence number", required=False, default=10, type=int, help="Number of nucleotides or aminoacids to add between the merged fasta sequences. If Sequence type is none then ignore this option")
    ap.add_argument("-id", "--FASTA header", required=True, type=str, help="FASTA header of the output file")
    ap.add_argument("-wi", "--FASTA width", required=False, type=int, default=80, choices=[60,70,80,100,120], help="Number of fasta sequence characters per line")
    ap.add_argument("-sfa", "--Single FASTA", required=True, widget="FileSaver", help="Output single-fasta file")
    args = vars(ap.parse_args())
    # main
    # index fasta file
    features = Fasta(args["Multi FASTA"], as_raw=True)
    # store all sequences to a list
    sequences = map(lambda key: features[key][:], features.keys())
    # select Sequence type
    if args["Sequence type"]=="none":
        spacer = ""
    elif args["Sequence type"]=="DNA":
        spacer = "N"*args["Sequence number"]
    else:
        spacer = "X"*args["Sequence number"]
    # merge all sequences at 1 
    merged_seqs = spacer.join(sequences)
    # export to fasta
    with open(args["Single FASTA"], "w") as f:
        f.write(f'>{args["FASTA header"]}\n{textwrap.fill(merged_seqs, width=args["FASTA width"])}')

if __name__ == "__main__":
    main()

