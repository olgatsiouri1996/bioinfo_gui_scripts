# python3
from gooey import *
import os
from dnachisel import *
# main
@Gooey(required_cols=1, program_name='sequence optimize plasmids', header_bg_color= '#DCDCDC', terminal_font_color= '#DCDCDC', terminal_panel_color= '#DCDCDC')
def main():
	ap = GooeyParser()
	ap.add_argument("-dir", "--directory", required=True, type=str, widget='DirChooser', help="directory with genbank files")
	args = vars(ap.parse_args())
# import each genbank file from the working directory
	for filename in sorted(os.listdir(os.chdir(args['directory']))):
		if filename.endswith(".gb") or filename.endswith(".gbk"):
# optimize the sequence from each file
			problem = DnaOptimizationProblem.from_record(filename)
			problem.resolve_constraints()
			problem.optimize()
# create biopython type sequence record from the optimized sequence and export to genbank
			problem.to_record(filepath= "".join([filename.split(".")[0],"_opt",".gb"]), with_sequence_edits=False)

if __name__ == "__main__":
	main()
