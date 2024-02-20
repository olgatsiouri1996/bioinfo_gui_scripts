import glob
import pandas as pd
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from Bio.PDB import *
import warnings
from Bio import BiopythonWarning

# Ignore Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)

def parse_pdb(filename):
    try:
        parser = PDBParser()
        s = parser.get_structure("name", filename)
        fill = s[0]
        dssp = DSSP(fill, filename, dssp='mkdssp')
        df = pd.DataFrame(dssp)
        df = df.loc[:, 2]
        struct_list = df.values.tolist()
        df1 = pd.DataFrame()
        df1['struct_list'] = struct_list
        df1 = df1['struct_list'].value_counts()
        df1 = round((df1 / df1.sum(axis=0)) * 100, 2)
        return df1
    except Exception as e:
        messagebox.showerror("Error", str(e))
        return None

def calculate_secondary_structure(input_dir, output_file):
    try:
        file_list = sorted(glob.glob(input_dir + '/*.pdb'))
        stats = pd.concat(map(parse_pdb, file_list), axis=1)
        if stats.empty:
            raise ValueError("No valid PDB files found in the directory.")
        stats.columns = list(map(lambda filename: filename.split("/")[-1].split(".pdb")[0], file_list))
        stats = stats.fillna(0)
        stats_t = stats.T
            
        with open(output_file, 'w') as f:
            f.write(
                stats_t.to_csv(header=True, index=True, sep="\t", doublequote=False, lineterminator='\n')
            )
        
        messagebox.showinfo("Analysis Complete", "Secondary structure analysis is complete. Results are saved to {}".format(output_file))
    except Exception as e:
        messagebox.showerror("Error", str(e))

def browse_input_dir():
    try:
        input_dir = filedialog.askdirectory()
        entry_input_dir.delete(0, tk.END)
        entry_input_dir.insert(0, input_dir)
    except Exception as e:
        messagebox.showerror("Error", str(e))

def browse_output_file():
    try:
        output_file = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])
        entry_output_file.delete(0, tk.END)
        entry_output_file.insert(0, output_file)
    except Exception as e:
        messagebox.showerror("Error", str(e))

def run_analysis():
    try:
        input_dir = entry_input_dir.get()
        output_file = entry_output_file.get()
        if not input_dir:
            raise ValueError("Please select input directory.")
        if not output_file:
            raise ValueError("Please select output file.")
        calculate_secondary_structure(input_dir, output_file)
    except Exception as e:
        messagebox.showerror("Error", str(e))

# Create GUI
root = tk.Tk()
root.title("PDB Secondary Structure Analysis")

frame = tk.Frame(root)
frame.pack(padx=10, pady=10)

label_input_dir = tk.Label(frame, text="Input Directory:")
label_input_dir.grid(row=0, column=0, sticky="w")

entry_input_dir = tk.Entry(frame, width=40)
entry_input_dir.grid(row=0, column=1, padx=5, pady=5)

button_browse_input = tk.Button(frame, text="Browse", command=browse_input_dir)
button_browse_input.grid(row=0, column=2, padx=5, pady=5)

label_output_file = tk.Label(frame, text="Output File:")
label_output_file.grid(row=1, column=0, sticky="w")

entry_output_file = tk.Entry(frame, width=40)
entry_output_file.grid(row=1, column=1, padx=5, pady=5)

button_browse_output = tk.Button(frame, text="Browse", command=browse_output_file)
button_browse_output.grid(row=1, column=2, padx=5, pady=5)

button_run = tk.Button(frame, text="Run Analysis", command=run_analysis)
button_run.grid(row=2, column=1, pady=10)

root.mainloop()
