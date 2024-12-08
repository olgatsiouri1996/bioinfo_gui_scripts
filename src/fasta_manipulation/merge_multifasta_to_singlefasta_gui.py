# python3
from gooey import *
from pyfaidx import Fasta
import textwrap
# imput parameters
@Gooey(required_cols=3, program_name="merge multi-fasta to single-fasta",default_size=(800, 530), header_bg_color= "#DCDCDC", terminal_font_color= "#DCDCDC", terminal_panel_color= "#DCDCDC")
def main():
    ap = GooeyParser()
    ap.add_argument("-mfa", "--Multi FASTA", required=True, widget="FileChooser", help="Input multi-FASTA file to merge its sequences")
    ap.add_argument("-type", "--Sequence type", required=False, type=str, default="DNA",choices=["DNA", "aa"], help="Sequence type")
    ap.add_argument("-num", "--Unknown number", required=False, default=0, type=int, widget="IntegerField", help="Number of NNs/XXs to add between the merged fasta sequences")
    ap.add_argument("-id", "--FASTA header", required=True, type=str, help="FASTA header of the output file")
    ap.add_argument("-wi", "--FASTA width", required=False, type=int, default=80, choices=[60,70,80,100,120], help="Number of fasta sequence characters per line")
    ap.add_argument("-sfa", "--Single FASTA", required=True, widget="FileSaver", gooey_options = {'wildcard': "FASTA files |*.fasta"}, help="Output single-fasta file")
    args = vars(ap.parse_args())
    # main
    # index fasta file
    features = Fasta(args["Multi FASTA"], as_raw=True)
    # store all sequences to a list
    sequences = map(lambda key: features[key][:], features.keys())
    # select Sequence type
    if args["Unknown number"]==0:
        spacer = ""
    else:
        if args["Sequence type"]=="DNA":
            spacer = "N"*args["Unknown number"]
        else:
            spacer = "X"*args["Unknown number"]
    # merge all sequences at 1 
    merged_seq = spacer.join(sequences)
    # wrap sequence
    wrapped_seq = textwrap.fill(merged_seq, width=args["FASTA width"])
    # create fasta format
    fasta_formatted = f'>{args["FASTA header"]}\n{wrapped_seq}'
    # export to fasta
    with open(args["Single FASTA"], "w") as f:
        f.write(fasta_formatted)

if __name__ == "__main__":
    main()

