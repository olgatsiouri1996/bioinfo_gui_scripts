import tkinter as tk
from tkinter import filedialog, messagebox

def calculate_u_production_rate(s, a1, b):
    return a1 / (1 + s**b)

def calculate_v_production_rate(s, a2, g):
    return a2 / (1 + s**g)

def process_input_file(input_filename, a1, a2, b, g, output_filename):
    try:
        with open(input_filename, 'r') as file:
            tf1, tf2 = zip(*(map(float, line.split()) for line in file))
    except Exception as e:
        messagebox.showerror("Error", f"Error reading input file: {str(e)}")
        return

    try:
        with open(output_filename, 'w') as output_file:
            output_file.write(f'rate_of_protein1_production\trate_of_protein2_production\tprotein1_concentration\tprotein2_concentration\n')

            for values in zip(
                map(lambda s: calculate_u_production_rate(s, a1, b), tf2),
                map(lambda s: calculate_v_production_rate(s, a2, g), tf1),
                tf1,
                tf2
            ):
                output_file.write('\t'.join(f"{x:.3f}" for x in values) + '\n')

        messagebox.showinfo("Success", "Calculation completed successfully. Output written to: " + output_filename)
    except Exception as e:
        messagebox.showerror("Error", f"Error writing output file: {str(e)}")

def open_file_dialog(entry_var):
    filename = filedialog.askopenfilename()
    entry_var.set(filename)

def save_file_dialog(entry_var):
    filename = filedialog.asksaveasfilename(defaultextension=".txt")
    entry_var.set(filename)

def run_calculation():
    input_filename = input_entry.get()
    output_filename = output_entry.get()

    if not input_filename or not output_filename:
        messagebox.showerror("Error", "Please provide both input and output file names.")
        return

    a1 = float(alpha1_entry.get())
    a2 = float(alpha2_entry.get())
    b = int(beta_entry.get())
    g = int(gamma_entry.get())

    process_input_file(input_filename, a1, a2, b, g, output_filename)

# GUI setup
root = tk.Tk()
root.title("Transcriptional Repressor Rate of Production Calculator")

# Input File
input_label = tk.Label(root, text="Input File")
input_label.pack()
input_entry = tk.Entry(root)
input_entry.pack()
input_button = tk.Button(root, text="Browse", command=lambda: input_entry.insert(0, filedialog.askopenfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])))
input_button.pack()

# Alpha1
alpha1_label = tk.Label(root, text="a1")
alpha1_label.pack()
alpha1_entry = tk.Entry(root)
alpha1_entry.pack()

# Alpha2
alpha2_label = tk.Label(root, text="a2")
alpha2_label.pack()
alpha2_entry = tk.Entry(root)
alpha2_entry.pack()

# Beta
beta_label = tk.Label(root, text="β")
beta_label.pack()
beta_entry = tk.Entry(root)
beta_entry.pack()

# Gamma
gamma_label = tk.Label(root, text="γ")
gamma_label.pack()
gamma_entry = tk.Entry(root)
gamma_entry.pack()

# Output File
output_label = tk.Label(root, text="Output File")
output_label.pack()
output_entry = tk.Entry(root)
output_entry.pack()
output_button = tk.Button(root, text="Browse", command=lambda: output_entry.insert(0, filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])))
output_button.pack()

# Run Button
run_button = tk.Button(root, text="Run Calculation", command=run_calculation)
run_button.pack()

root.mainloop()
