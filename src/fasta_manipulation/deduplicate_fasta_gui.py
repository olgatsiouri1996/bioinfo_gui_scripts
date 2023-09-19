import textwrap
from pyfaidx import Fasta
from gooey import Gooey, GooeyParser

@Gooey(program_name='Deduplicate FASTA by Identifier and/or Sequence', default_size=(840, 440), header_bg_color='#DCDCDC', terminal_font_color='#DCDCDC', terminal_panel_color='#DCDCDC')
def main():
    parser = GooeyParser(description="Remove duplicates from a FASTA file, by keeping the first occurrence based on Identifier and/or Sequence")
    parser.add_argument('input file', widget="FileChooser", help="Input FASTA file")
    parser.add_argument('output file', widget="FileSaver", help="Output file name")
    parser.add_argument('deduplication option', widget="Dropdown", choices=["Identifier", "Sequence", "Both"], default="Identifier", help="Deduplication Option")
    parser.add_argument('output option', widget="Dropdown", choices=["multi-fasta", "3-Column txt with identifier, description and sequence in each row", "1-Column txt with identifier in each row"], default="multi-fasta", help="Output Option")

    args = vars(parser.parse_args())
    
    if args['deduplication option'] == "Sequence":
        sequences = Fasta(args['input file'])
    else:
        sequences = Fasta(args['input file'], duplicate_action="first")
    
    def write_record(out_f, record):
        if args['output option'] == "multi-fasta":
            wrapped_sequence = textwrap.fill(str(record), width=60)
            out_f.write(f'>{record.long_name}\n{wrapped_sequence}\n')
        elif args['output option'] == "3-Column txt with identifier, description and sequence in each row":
            description = str(record.long_name).split(None, 1)[1] if len(str(record.long_name).split(None, 1)) > 1 else ""
            out_f.write(f'{record.name}\t{description}\t{str(record)}\n')
        else:
            out_f.write(f'{record.name}\n')

    seen_records = set()
    with open(args['output file'], 'w') as out_f:

        if args['output option'] == "3-Column txt with identifier, description and sequence in each row":
            out_f.write("id\tdescription\tsequence\n")

        for record in sequences:
            sequence = str(record)

            if args['deduplication option'] == "Identifier":
                write_record(out_f, record)
            else:
                if sequence not in seen_records:
                    seen_records.add(sequence)
                    write_record(out_f, record)

if __name__ == "__main__":
    main()
