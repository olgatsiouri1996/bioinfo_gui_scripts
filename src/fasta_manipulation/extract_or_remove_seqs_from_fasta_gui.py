import os
import tkinter as tk
from tkinter import filedialog, messagebox
from pyfaidx import Fasta
import textwrap

def enable_output_file_widgets():
    output_entry.config(state=tk.NORMAL)
    directory_entry.config(state=tk.DISABLED)
    prefix_entry.config(state=tk.DISABLED)
    num_entry.config(state=tk.DISABLED)

def enable_output_directory_widgets():
    output_entry.config(state=tk.DISABLED)
    directory_entry.config(state=tk.NORMAL)
    prefix_entry.config(state=tk.DISABLED)
    num_entry.config(state=tk.DISABLED)

def enable_output_directory_prefix_widgets():
    output_entry.config(state=tk.DISABLED)
    directory_entry.config(state=tk.NORMAL)
    prefix_entry.config(state=tk.NORMAL)
    num_entry.config(state=tk.NORMAL)

def update_output_widgets(*args):
    selected_output_type = output_type_var.get()
    
    if selected_output_type in [
        '1 multi-fasta file', '2-column txt file(id,seq)',
        '3-column txt file(id,description,seq)',
        '3-column txt file(id,full_fasta_header,seq)',
        '2-column txt file(id,description)', '1-column txt file(id)'
    ]:
        enable_output_file_widgets()
    elif selected_output_type in ['many single-fasta files', 'many multi-fasta files 1 per group']:
        enable_output_directory_widgets()
    elif selected_output_type == 'many multi-fasta files by number of records':
        enable_output_directory_prefix_widgets()

def get_sequences():
    input_file = input_entry.get()
    ids_file = ids_entry.get()
    action = action_var.get()
    output_type = output_type_var.get()
    output_file = output_entry.get()
    output_directory = directory_entry.get()
    prefix = prefix_entry.get()
    num = num_entry.get()

    try:
        # import fasta file
        features = Fasta(input_file)
        # collect ids of import fasta file
        keyslist = (features.keys())
        # import the txt file with headers you want to extract the sequence from the input fasta
        txt = (str(line.rstrip()).split()[0] for line in open(ids_file))
        # merge txt file ids to import fasta file ids and keep the common ones
        headers=(set(txt).intersection(keyslist))
        # choose to remove sequences
        if action == 'remove':
            # remove those that exist in the headers generator
            final_keys = (set(keyslist).difference(headers))
        else:
            final_keys = headers
        
        match output_type:
            case '1 multi-fasta file':
                # export to 1 multi-fasta
                with  open(output_file, 'w') as f:
                    for key in final_keys:
                        f.write(f'>{str(features[key].long_name).rstrip()}\n{textwrap.fill(features[key][:].seq, width=60)}\n')                
            case 'many single-fasta files':
                # extract many single-fasta files
                os.chdir(output_directory)
                for key in final_keys:
                    with open(f'{key}.fasta', 'w') as f:
                        f.write(f'>{str(features[key].long_name).rstrip()}\n{textwrap.fill(features[key][:].seq, width=60)}\n')
            case 'many multi-fasta files 1 per group':
                # create a dictionary with identifier and group
                identifier_to_group = {}
                groups = (str(line.rstrip()).split()[1] for line in open(ids_file))
                for (identifier, group) in zip(headers, groups):
                    identifier_to_group[identifier] = group
                # populate group_to_sequences dictionary
                group_to_sequences = {}
                for identifier, group in identifier_to_group.items():
                    if group not in group_to_sequences:
                        group_to_sequences[group] = []
                    sequence = textwrap.fill(features[str(identifier)][:].seq, width=60)
                    description = features[str(identifier)].long_name
                    group_to_sequences[group].append((description, sequence))
                # write sequences to group-specific output files
                os.chdir(output_directory)
                for group, sequences in group_to_sequences.items():
                    output_file_name = f"{group}.fasta"
                    with open(output_file_name, "w") as output_file:
                        for description, sequence in sequences:
                            output_file.write(f">{str(description).rstrip()}\n{sequence}\n")
            case 'many multi-fasta files by number of records':
                split_size = int(num)
                split_lists = [list(final_keys)[i:i+split_size] for i in range(0, len(list(final_keys)), split_size)]
                for i, split_identifiers in enumerate(split_lists, 1):
                    output_file = os.path.join(output_directory, f"{prefix}_part{i}.fasta")
                    with open(output_file, 'w') as output:
                        for key in split_identifiers:
                            header = str(features[key].long_name).rstrip()
                            output.write(f">{header}\n")
                            seq = str(features[key][:].seq).rstrip()
                            wrapped_seq = textwrap.fill(seq, width=60)
                            output.write(wrapped_seq + '\n')
            case '2-column txt file(id,seq)':
                # Iterate input headers to extract sequences and export as a 2-column txt file
                with open(output_file, 'w') as f:
                    f.write("id\tseq\n")
                    for key in final_keys:
                        f.write(f"{key}\t{features[key][:].seq}\n")
            case '3-column txt file(id,description,seq)':
                # Iterate input headers to extract sequences and export as a 3-column txt file
                with open(output_file, 'w') as f:
                    f.write("id\tdescription\tseq\n")
                    for key in final_keys:
                        try:
                            f.write(f"{key}\t{str(features[key].long_name).rstrip().split(' ', 1)[1]}\t{features[key][:].seq}\n")
                        except IndexError:
                            f.write(f"{key}\t{''}\t{features[key][:].seq}\n")
            case '3-column txt file(id,full_fasta_header,seq)':
                # Iterate input headers to extract sequences and export as a 3-column txt file
                with open(output_file, 'w') as f:
                    f.write("id\tfull_fasta_header\tseq\n")
                    for key in final_keys:
                        f.write(f"{key}\t{str(features[key].long_name).rstrip()}\t{features[key][:].seq}\n")
            case '2-column txt file(id,description)':
                # Iterate input headers to extract sequences and export as a 2-column txt file
                with open(output_file, 'w') as f:
                    f.write("id\tdescription\n")
                    for key in final_keys:
                        try:
                            f.write(f"{key}\t{str(features[key].long_name).rstrip().split(' ', 1)[1]}\n")
                        except IndexError:
                            f.write(f"{key}\t{''}\n")
            case '1-column txt file(id)':
                # Iterate input headers to extract identifiers and export as a 1-column txt file
                with open(output_file, 'w') as f:
                    for key in final_keys:
                        f.write(f"{key}\n")

        # Notify when processing is finished
        messagebox.showinfo("Processing Complete", "The project is finished.")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {str(e)}")

# Create the main Tkinter window
root = tk.Tk()
root.title("Sequence Extractor or Remover")

# Input file
input_label = tk.Label(root, text="Input multi-fasta file")
input_label.pack()
input_entry = tk.Entry(root)
input_entry.pack()
input_button = tk.Button(root, text="Browse", command=lambda: input_entry.insert(0, filedialog.askopenfilename()))
input_button.pack()

# IDs file
ids_label = tk.Label(root, text="IDs file")
ids_label.pack()
ids_entry = tk.Entry(root)
ids_entry.pack()
ids_button = tk.Button(root, text="Browse", command=lambda: ids_entry.insert(0, filedialog.askopenfilename()))
ids_button.pack()

# Action
action_label = tk.Label(root, text="Action")
action_label.pack()
action_var = tk.StringVar()
action_var.set('extract')  # Default selection
action_options = ['extract', 'remove']
action_menu = tk.OptionMenu(root, action_var, *action_options)
action_menu.pack()

# Output type
output_type_label = tk.Label(root, text="Output type")
output_type_label.pack()
output_type_var = tk.StringVar()
output_type_var.set('1 multi-fasta file')  # Default selection
output_type_options = [
    '1 multi-fasta file', '2-column txt file(id,seq)', '3-column txt file(id,description,seq)',
    '3-column txt file(id,full_fasta_header,seq)', '2-column txt file(id,description)', '1-column txt file(id)',
    'many single-fasta files', 'many multi-fasta files 1 per group', 'many multi-fasta files by number of records'
]
output_type_menu = tk.OptionMenu(root, output_type_var, *output_type_options, command=update_output_widgets)
output_type_menu.pack()

# Output file
output_label = tk.Label(root, text="Output file")
output_label.pack()
output_entry = tk.Entry(root)
output_entry.pack()
output_extension_var = tk.StringVar()
output_extension_var.set('fasta')  # Default selection
output_extension_options = ['fasta', 'txt']
output_extension_menu = tk.OptionMenu(root, output_extension_var, *output_extension_options)
output_extension_menu.pack()
output_button = tk.Button(root, text="Browse", command=lambda: output_entry.insert(0, filedialog.asksaveasfilename(defaultextension=".fasta", filetypes=[("Fasta files", "*.fasta"), ("Text files", "*.txt")])))
output_button.pack()

# Output directory
directory_label = tk.Label(root, text="Output directory")
directory_label.pack()
directory_entry = tk.Entry(root)
directory_entry.pack()
directory_button = tk.Button(root, text="Browse", command=lambda: directory_entry.insert(0, filedialog.askdirectory()))
directory_button.pack()

# File Prefix
prefix_label = tk.Label(root, text="File Prefix")
prefix_label.pack()
prefix_entry = tk.Entry(root)
prefix_entry.pack()

# Number of records for multi-fasta files
num_label = tk.Label(root, text="Number of records per file")
num_label.pack()
num_entry = tk.Entry(root)
num_entry.pack()

# Extract button
extract_button = tk.Button(root, text="Get Sequences", command=get_sequences)
extract_button.pack()

# Initialize the output widgets based on the default value
update_output_widgets()

root.mainloop()
