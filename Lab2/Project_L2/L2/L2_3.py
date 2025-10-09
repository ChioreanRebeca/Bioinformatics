# Design an app using AI which contains GUI which allows the user to select a Fasta File. 
# The content of the file should be analized using a sliding window of 30 positions.
# The content for each sliding window should be used in order to extract the relative 
# frequences of the symbols found in the alphabet of the sequence.
# Thus your input would be the DNA sequence fom the fasta file and the output should be
#  the values of the relative frequencies of each symbol in the sequence.
# Translate in lines on a chart, thus your chart in the case of DNA should have 4 lines 
# which reflect the values found over the sequence.



import tkinter as tk
from tkinter import filedialog, scrolledtext
import matplotlib.pyplot as plt

WINDOW_SIZE = 30

def read_fasta_file(filepath):
    """Reads a FASTA file and returns only the sequence (ignoring header)."""
    with open(filepath, "r") as f:
        lines = f.readlines()
    sequence = "".join(line.strip() for line in lines[1:])
    return sequence

def sliding_window_frequencies(sequence):
    """Calculate relative frequencies per sliding window."""
    alphabet = sorted(set(sequence))
    freqs = {letter: [] for letter in alphabet}
    
    for i in range(len(sequence) - WINDOW_SIZE + 1):
        window = sequence[i:i+WINDOW_SIZE]
        total = len(window)
        count_dict = {}
        for s in window:
            count_dict[s] = count_dict.get(s, 0) + 1
        for letter in alphabet:
            freqs[letter].append(count_dict.get(letter, 0)/total)
    
    return freqs

def plot_frequencies(freqs):
    """Plot the relative frequencies for each symbol."""
    plt.figure(figsize=(10, 6))
    x = range(len(next(iter(freqs.values()))))
    for letter, values in freqs.items():
        plt.plot(x, values, label=letter)
    plt.xlabel("Window Start Position")
    plt.ylabel("Relative Frequency")
    plt.title("Sliding Window Symbol Frequencies")
    plt.legend()
    plt.show()

def count_di_tri(sequence):
    """Count dinucleotides and trinucleotides with brute force."""
    dinucleoids = {}
    trinucleoids = {}
    
    for i in range(len(sequence)-1):
        dinuc = sequence[i:i+2]
        dinucleoids[dinuc] = dinucleoids.get(dinuc, 0) + 1
    for i in range(len(sequence)-2):
        trinuc = sequence[i:i+3]
        trinucleoids[trinuc] = trinucleoids.get(trinuc, 0) + 1


    total_dn = 0
    total_tn = 0

    for d in dinucleoids.keys():
        total_dn += dinucleoids[d]

    for t in trinucleoids.keys():
        total_tn += trinucleoids[t]

    dinuc_perc = {d: round(dinucleoids[d]/total_dn*100, 2) for d in dinucleoids}
    trinuc_perc = {t: round(trinucleoids[t]/total_tn*100, 2) for t in trinucleoids}
    return dinuc_perc, trinuc_perc

def load_fasta():
    filepath = filedialog.askopenfilename(
        title="Select a FASTA file",
        filetypes=(("FASTA files", "*.fasta *.fa"), ("All files", "*.*"))
    )
    if filepath:
        sequence = read_fasta_file(filepath)
        dinuc, trinuc = count_di_tri(sequence)
        freqs = sliding_window_frequencies(sequence)
        
        # Display text output
        text_box.delete(1.0, tk.END)
        text_box.insert(tk.END, f"FASTA Sequence:\n{sequence}\n\n")
        text_box.insert(tk.END, "Dinucleotide counts:\n")
        for k, v in sorted(dinuc.items()):
            text_box.insert(tk.END, f"{k}: {v}\n")
        text_box.insert(tk.END, "\nTrinucleotide counts:\n")
        for k, v in sorted(trinuc.items()):
            text_box.insert(tk.END, f"{k}: {v}\n")
        
        # Plot chart
        plot_frequencies(freqs)

# GUI Setup
root = tk.Tk()
root.title("FASTA Analyzer")

btn = tk.Button(root, text="Load FASTA File", command=load_fasta)
btn.pack(pady=10)

text_box = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=80, height=25)
text_box.pack(padx=10, pady=10)

root.mainloop()
