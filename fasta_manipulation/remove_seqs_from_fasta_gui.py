# python3
from gooey import *
from Bio import SeqIO
# input parameters
@Gooey(required_cols=3, program_name= 'remove sequences from fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
  ap = GooeyParser(description="use a txt file with fasta headers to remove sequences from a fasta file")
  ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input fasta file")
  ap.add_argument("-out", "--output", required=True, widget='FileSaver', help="output fasta file")
  ap.add_argument("-headers", "--headers", required=True, widget='FileChooser', help="file with fasta headers to exlude from the output fasta sequences")
  args = vars(ap.parse_args())
# main
  wanted = set()
  with open(args['headers']) as f:
    for line in f:
        line = line.strip()
        if line != "":
            wanted.add(line)

  fasta_sequences = SeqIO.parse(open(args['input']),'fasta')
  with open(args['output'], "w") as f:
    for seq in fasta_sequences:
        if seq.id not in wanted:
            SeqIO.write([seq], f, "fasta")

if __name__ == '__main__':
    main()
