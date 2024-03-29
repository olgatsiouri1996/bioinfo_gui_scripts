import tkinter as tk
from tkinter import filedialog, messagebox
import numpy as np
import math

def calculate_probabilities(input_file, calculation_type, round_type, exponent_value, output_file):
    num = int(num_entry.get())
    
    # Load numeric data from the input file into a NumPy array
    if calculation_type =='factorial':
        dat = np.loadtxt(input_file,dtype=int)
    else:
        dat = np.loadtxt(input_file)

    if round_type == 'yes':
        forma = '%.{}f'.format(num)
    else:
        forma = '%s'
    
    try:
        match calculation_type:            
            case 'sin':
                res = np.sin(dat)
            case 'cos':
                res = np.cos(dat)
            case 'tan':
                res = np.tan(dat)
            case 'inverse sin':
                res = np.arcsin(dat)
            case 'inverse cos':
                res = np.arccos(dat)
            case 'inverse tan':
                res = np.arctan(dat)
            case 'sinh':
                res = np.sinh(dat)
            case 'cosh':
                res = np.cosh(dat)
            case 'tanh':
                res = np.tanh(dat)
            case 'inverse sinh':
                res = np.arcsinh(dat)
            case 'inverse cosh':
                res = np.arccosh(dat)
            case 'inverse tanh':
                res = np.arctanh(dat) 
            case 'ln':
                res = np.log(dat)
            case 'log2':
                res = np.log2(dat)
            case 'log10':
                res = np.log10(dat)
            case 'e^x':
                res = np.exp(dat)
            case 'square root':
                res = np.sqrt(dat)
            case 'cubic root':
                res = np.cbrt(dat)
            case 'x^n':
                res = np.power(dat, exponent_value)
            case 'x^1/n':
                res = np.power(dat, 1/exponent_value)
            case 'n^x':
                res = np.power(exponent_value, dat)
            case 'factorial':
                res =  np.array(list(map(math.factorial, dat)))
            case 'x+n':
                res =  dat + exponent_value
            case 'x*n':
                res =  dat*exponent_value
            case 'x-n':
                res = dat - exponent_value
            case 'n-x':
                res = exponent_value - dat
            case 'x/n':
                res = dat/exponent_value
            case 'n/x':
                res = exponent_value/dat
            case 'x+y':
                res =  dat[:, 0] + dat[:, 1]
            case 'x*y':
                res =  dat[:, 0]*dat[:, 1]
            case 'x-y':
                res = dat[:, 0] - dat[:, 1]
            case 'y-x':
                res = dat[:, 1] - dat[:, 0]
            case 'x/y':
                res = dat[:, 0]/dat[:, 1]
            case 'y/x':
                res = dat[:, 1]/dat[:, 0]
            case 'x^y':
                res = np.power(dat[:, 0], dat[:, 1])
            case 'y^x':
                res = np.power(dat[:, 1], dat[:, 0])
            case 'cumulative sum':
                res = np.cumsum(dat)
                        
        # Save results to the output file
        np.savetxt(output_file, res, delimiter='\t', fmt=forma)
        messagebox.showinfo("Success", "Calculation completed successfully")
    except Exception as e:
        messagebox.showerror("Error", f"Error writing output file: {str(e)}")

def open_file(entry):
    file_path = filedialog.askopenfilename()
    entry.delete(0, tk.END)
    entry.insert(tk.END, file_path)

def save_file(entry):
    file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])
    entry.delete(0, tk.END)
    entry.insert(tk.END, file_path)

def toggle_entries(calculation_type):
    if calculation_type in ['x^n','x^1/n','n^x','x+n','x*n','x-n','n-x','x/n','n/x']:
        exponent_label.config(state=tk.NORMAL)
        exponent_entry.config(state=tk.NORMAL)
    else:
        exponent_label.config(state=tk.DISABLED)
        exponent_entry.config(state=tk.DISABLED)

def toggle_round(round_type):
    if round_type == 'yes':
        num_label.config(state=tk.NORMAL)
        num_entry.config(state=tk.NORMAL)
    else:
        num_label.config(state=tk.DISABLED)
        num_entry.config(state=tk.DISABLED)

def main():
    root = tk.Tk()
    root.title('Multi Value Scientific Calculator')
    # Program description
    description_label = tk.Label(root)
    description_label.pack()

    # Input file
    input_frame = tk.Frame(root)
    input_label = tk.Label(input_frame, text='Input 1/2-Column txt File(no headers):')
    input_label.pack(side='left')
    input_entry = tk.Entry(input_frame)
    input_entry.pack(side='left')
    input_button = tk.Button(input_frame, text='Browse', command=lambda: open_file(input_entry))
    input_button.pack(side='left')
    input_frame.pack()

    # Calculation type
    calculation_type_frame = tk.Frame(root)
    calculation_type_label = tk.Label(calculation_type_frame, text='Calculation Type:')
    calculation_type_label.pack(side='left')
    calculation_type_var = tk.StringVar(root)
    calculation_type_var.set('sin')
    calculation_type_dropdown = tk.OptionMenu(calculation_type_frame, calculation_type_var, 'sin', 'cos', 'tan','inverse sin', 'inverse cos', 'inverse tan','sinh', 'cosh', 'tanh','inverse sinh', 'inverse cosh', 'inverse tanh','ln','log2','log10','e^x','square root','cubic root','x^n','x^1/n','n^x','factorial','x+n','x*n','x-n','n-x','x/n','n/x','x+y','x*y','x-y','y-x','x/y','y/x','x^y','y^x','cumulative sum', command=toggle_entries)
    calculation_type_dropdown.pack(side='left')
    calculation_type_frame.pack()

    round_type_frame = tk.Frame(root)
    round_type_label = tk.Label(round_type_frame, text='Round Data:')
    round_type_label.pack(side='left')
    round_type_var = tk.StringVar(root)
    round_type_var.set('yes')
    round_type_dropdown = tk.OptionMenu(round_type_frame, round_type_var, 'yes', 'no', command=toggle_round)
    round_type_dropdown.pack(side='left')
    round_type_frame.pack()

    # Exponent value for power or root
    exponent_frame = tk.Frame(root)
    global exponent_label
    exponent_label = tk.Label(exponent_frame, text='n:', state=tk.DISABLED)
    exponent_label.pack(side='left')
    global exponent_entry
    exponent_entry = tk.Entry(exponent_frame, state=tk.DISABLED)
    exponent_entry.pack(side='left')
    exponent_frame.pack()

    # Number of records for multi-fasta files
    global num_label
    num_label = tk.Label(root, text="Number of digits after the decimal point:",state=tk.NORMAL)
    num_label.pack()
    global num_entry
    num_entry = tk.Entry(root)
    num_entry.insert(0, "3")
    num_entry.pack()

    # Output file
    output_frame = tk.Frame(root)
    output_label = tk.Label(output_frame, text='Output 1-Column txt File(no headers):')
    output_label.pack(side='left')
    output_entry = tk.Entry(output_frame)
    output_entry.pack(side='left')
    output_button = tk.Button(output_frame, text='Browse', command=lambda: save_file(output_entry))
    output_button.pack(side='left')
    output_frame.pack()

    # Run button
    run_button = tk.Button(root, text='Run', command=lambda: calculate_probabilities(
        input_entry.get(),
        calculation_type_var.get(),
        round_type_var.get(),
        float(exponent_entry.get() or 0),
        output_entry.get()
    ))
    run_button.pack()

    root.mainloop()

if __name__ == '__main__':
    main()
