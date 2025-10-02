"""
Use the A.I to design an app with a graphical user interface that is able to 
integrate your alg from assignment 1 and 2. Your app must have a button which
allows the user to choose a Fasta File.
Fasta Files contain a specific biological format. The output should be shown 
on the main window by using a TextBox obj or something similair.

Fasta Files have the following format:
1. the first line is an information line that shows the ID of teh seq. and other type of information
2. the starting from 2nd line we have the rows seq. which can be DNA, ARN or proteins that is split 
in 80 chars lines until the end of the file 

Use the AI to simmulate a fasta file for the input

"""

import tkinter as tk
from tkinter import filedialog, scrolledtext


def exercise_one_logic(sequence):
    """Counts occurrences of each letter in the sequence."""
    alphabet_dict = {}
    for s in sequence:
        alphabet_dict[s] = alphabet_dict.get(s, 0) + 1
    
    result = "Exercise 1 (Alphabet):\n"
    for letter in sorted(alphabet_dict.keys()):
        result += f"{letter}\n"
    return result


def exercise_two_logic(sequence):
    """Calculates relative frequencies of letters in the sequence."""
    alphabet_dict = {}
    for s in sequence:
        alphabet_dict[s] = alphabet_dict.get(s, 0) + 1

    total = sum(alphabet_dict.values())
    
    result = "Exercise 2 (Frequencies):\n"
    for letter in sorted(alphabet_dict.keys()):
        percentage = alphabet_dict[letter] / total * 100
        result += f"{letter}: {percentage:.2f}%\n"
    return result


def read_fasta_file(filepath):
    """Reads a FASTA file and returns only the sequence (ignoring header)."""
    with open(filepath, "r") as f:
        lines = f.readlines()
    
    # First line is header, rest is sequence
    sequence = "".join(line.strip() for line in lines[1:])
    return sequence


def load_fasta():
    filepath = filedialog.askopenfilename(
        title="Select a FASTA file",
        filetypes=(("FASTA files", "*.fasta *.fa"), ("All files", "*.*"))
    )
    if filepath:
        sequence = read_fasta_file(filepath)
        output = ""
        output += exercise_one_logic(sequence) + "\n"
        output += exercise_two_logic(sequence)
        text_box.delete(1.0, tk.END)
        text_box.insert(tk.END, f"FASTA Sequence:\n{sequence}\n\n{output}")


# GUI Setup
root = tk.Tk()
root.title("FASTA Analyzer")

btn = tk.Button(root, text="Load FASTA File", command=load_fasta)
btn.pack(pady=10)

text_box = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=60, height=20)
text_box.pack(padx=10, pady=10)

root.mainloop()
