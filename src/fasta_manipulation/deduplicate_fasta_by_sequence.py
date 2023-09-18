import textwrap
from pyfaidx import Fasta
from gooey import Gooey, GooeyParser

@Gooey(program_name='Deduplicate FASTA by Sequence', default_size=(610, 440), header_bg_color='#DCDCDC', terminal_font_color='#DCDCDC', terminal_panel_color='#DCDCDC')
def main():
    parser = GooeyParser(description="Remove duplicates based on sequence from a FASTA file")
    parser.add_argument('input file', widget="FileChooser", help="Input FASTA file")
    parser.add_argument('output file', widget="FileSaver", help="Output file name")
    parser.add_argument('output format', widget="Dropdown", choices=["FASTA", "Tab Separated"], default="FASTA", help="Output format")

    args = vars(parser.parse_args())
    
    sequences = Fasta(args['input file'])

    if args['output format'] == "FASTA":
        seen_sequences = set()

    with open(args['output file'], 'w') as out_f:
        if args['output format'] == "FASTA":
            for record in sequences:
                sequence = str(record)
                if sequence not in seen_sequences:
                    seen_sequences.add(sequence)
                    wrapped_sequence = textwrap.fill(sequence, width=60)  # Wrap the sequence
                    out_f.write(f'>{record.long_name}\n{wrapped_sequence}\n')
        elif args['output format'] == "Tab Separated":
            out_f.write("id\tdescription\tsequence\n")
            seen_sequences = set()
            for record in sequences:
                sequence = str(record)
                if sequence not in seen_sequences:
                    seen_sequences.add(sequence)
                    # Extract the description from the long_name, or use an empty string if no whitespace is found
                    description = record.long_name.split(None, 1)[1] if len(record.long_name.split(None, 1)) > 1 else ""
                    out_f.write(f'{record.name}\t{description}\t{sequence}\n')

if __name__ == "__main__":
    main()
