import tkinter as tk
from tkinter import filedialog, messagebox
import numpy as np

def calculate_probabilities(input_file, calculation_type, allele_choice, output_file):
    num = int(num_entry.get())
    
    # Load probabilities from the input file into a NumPy array
    probabilities = np.loadtxt(input_file)
    
    try:
        if allele_choice == 'other':
            probabilities = (1 - np.sqrt(probabilities)) ** 2
        
        if calculation_type == 'allele':
            res = np.column_stack((2 * np.sqrt(probabilities) - probabilities, probabilities, np.sqrt(probabilities)))
            header = 'probability_of_containing_allele\tprobability_of_homozygous_individuals\tallele_frequency\n'
        elif calculation_type == 'carrier':
            p = np.sqrt(probabilities)
            res = np.column_stack((2 * p * (1 - p), probabilities, p))
            header = 'probability_of_carriers\tprobability_of_homozygous_individuals\tallele_frequency\n'
        else:
            p = np.sqrt(probabilities)
            res = np.column_stack((2 * np.sqrt(probabilities) - probabilities, 2 * p * (1 - p), probabilities, p))
            header = 'probability_of_containing_allele\tprobability_of_carriers\tprobability_of_homozygous_individuals\tallele_frequency\n'

        # Save results to the output file
        np.savetxt(output_file, res, delimiter='\t', header=header, fmt='%.{}f'.format(num))
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

def main():
    root = tk.Tk()
    root.title('Genetic Probability Calculator')

    # Program description
    description_label = tk.Label(root, text='Calculate the probability of an individual carrying at least 1 autosomic allele copy and/or being heterozygous, using the probability of homozygous individuals in a given diploid population')
    description_label.pack()

    # Input file
    input_frame = tk.Frame(root)
    input_label = tk.Label(input_frame, text='Input File:')
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
    calculation_type_var.set('allele')
    calculation_type_dropdown = tk.OptionMenu(calculation_type_frame, calculation_type_var, 'allele', 'carrier', 'both')
    calculation_type_dropdown.pack(side='left')
    calculation_type_frame.pack()

    # Allele choice
    allele_frame = tk.Frame(root)
    allele_label = tk.Label(allele_frame, text='Allele Choice:')
    allele_label.pack(side='left')
    allele_var = tk.StringVar(root)
    allele_var.set('this')
    allele_dropdown = tk.OptionMenu(allele_frame, allele_var, 'this', 'other')
    allele_dropdown.pack(side='left')
    allele_frame.pack()

    # Number of records for multi-fasta files
    num_label = tk.Label(root, text="Number of digits after the decimal point:")
    num_label.pack()
    global num_entry
    num_entry = tk.Entry(root)
    num_entry.insert(0, "3")
    num_entry.pack()

    # Output file
    output_frame = tk.Frame(root)
    output_label = tk.Label(output_frame, text='Output File:')
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
        allele_var.get(),
        output_entry.get()
    ))
    run_button.pack()

    root.mainloop()

if __name__ == '__main__':
    main()
