import pandas as pd
import os
import warnings
from gooey import Gooey, GooeyParser

# Custom function to split values by semicolon and return the first part
def split_and_get_first(x):
    if isinstance(x, str):
        parts = x.split(';')
        if len(parts) > 0:
            return parts[0]
    return x

# Custom function to split values by equal sign and return the second part (if present)
def split_and_get_second(x):
    if isinstance(x, str):
        parts = x.split('=')
        if len(parts) > 1:
            return parts[1]
    return x

@Gooey(program_name="GFF3 to BED Converter",default_size=(610,370))
def main():
    parser = GooeyParser(description="Convert GFF3 to BED format")

    parser.add_argument("input file", widget="FileChooser", help="Select the GFF3 input file", gooey_options={"wildcard": "*.gff3;*.gff"})
    parser.add_argument("output dir", widget="DirChooser", help="Select the output directory")

    args = vars(parser.parse_args())

    # Extract the filename without the extension
    filename = os.path.splitext(os.path.basename(args['input file']))[0]

    # Suppress pandas warnings
    warnings.filterwarnings("ignore")

    # Read the GFF3 file into a DataFrame
    gff3 = pd.read_csv(args['input file'], sep='\t', header=None, comment='#')

    # Convert the first column (named '8') to strings
    gff3[8] = gff3[8].astype(str)

    # Select and reorder specific columns
    dat = gff3[[0, 3, 4, 8, 5, 6, 1, 2, 7]].copy()

    # Rename the first '8' column to '8.1' and the second '8' column to '8.2'
    dat = dat.rename(columns={8: '8.1'})
    gff3 = gff3.rename(columns={8: '8.2'})

    # Apply the custom function to split values in the '8.1' column by semicolon and get the first part
    dat['8.1'] = dat['8.1'].apply(split_and_get_first)
    dat['8.1'] = dat['8.1'].apply(split_and_get_second)

    # Append the '8.2' column to the final 'dat' DataFrame
    dat['8.2'] = gff3['8.2']

    # Subtract 1 from the V4 column
    dat[3] = dat[3] - 1

    # Construct the output file path
    output_file = os.path.join(args['output dir'], filename + '.bed')

    # Write the modified DataFrame to the BED file in the selected output directory
    dat.to_csv(output_file, sep='\t', header=False, index=False)

if __name__ == "__main__":
    main()
