# python3
import os
import sys
from gooey import *
from biopandas.pdb import PandasPdb
# input parameters
@Gooey(required_cols=2, program_name='subset many pdb to single or multi fasta', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
    ap = GooeyParser(description="converts each pdb files into single fasta files or makes a multi-fasta file based on the chain, start and end locations")
    ap.add_argument("-in", "--input", required=True, widget='FileChooser', help="input 4-column txt file with pdb filename(no extension),chain,start and end coordinates")
    ap.add_argument("-mfa", "--multifasta", required=False, widget='FileSaver', help="output multi-fasta file")
    ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser', help="directory to search for pdb files(the directory can contain many filetypes)")
    ap.add_argument("-out", "--output", required=False, type=str, widget='DirChooser', help="directory to save the single fasta files(the directory can contain many filetypes)")
    ap.add_argument("-type", "--type", required=False,default=1, type=int, widget='Dropdown', choices=[1,2],  help="type of output to choose: 1) 1 multi-fasta file, 2) many single-fasta files.")
    args = vars(ap.parse_args())
    # main
    # create function to split the input sequence based on a specific number of characters(60)
    def split_every_60(s): return [str(s)[i:i+60] for i in range(0,len(str(s)),60)]
    # convert a 4-column txt file to 4 generators 1 per column
    pdbname = (str(line.rstrip()).split()[0] for line in open(args['input']))
    chain = (str(line.rstrip()).split()[1] for line in open(args['input']))
    start = (int(str(line.rstrip()).split()[2]) for line in open(args['input']))
    end = (int(str(line.rstrip()).split()[3]) for line in open(args['input']))
    # create function to trim + convert to fasta
    def trim_to_fasta(fi,ch,st,en):
        # insert pdb file
        ppdb = PandasPdb()
        ppdb.read_pdb(os.path.join(args['directory'],''.join([fi,'.pdb'])))
        # convert 3 letters aa to 1
        one = ppdb.amino3to1()
        # select chain and convert to string
        sequence = ''.join(one.loc[one['chain_id']==ch,'residue_name'])
        # subset by location
        prot = sequence[int(st -1):en]
        # remove dataframes and lists
        del ppdb; del one
        return prot
    # select between exporting 1 or many fasta files
    if args['type'] == 1:
        sys.stdout = open(args['multifasta'],'a')
        for (a,b,c,d) in zip(pdbname,chain,start,end):
            print(''.join([">",a,"_",b,"_",str(c),"_",str(d)]).replace('\r',''))
            print('\n'.join(split_every_60(trim_to_fasta(a,b,c,d))))
        sys.stdout.close()
    else:
        os.chdir(args['output'])
        for (a,b,c,d) in zip(pdbname,chain,start,end):
            sys.stdout = open(''.join([a,"_",b,"_",str(c),"_",str(d),".fasta"]), 'a')
            print(''.join([">",a,"_",b,"_",str(c),"_",str(d)]).replace('\r',''))
            print('\n'.join(split_every_60(trim_to_fasta(a,b,c,d))))
            sys.stdout.close()            

         
if __name__ == '__main__':
    main()