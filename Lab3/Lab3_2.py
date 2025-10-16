"""
# Use the AI to design an application that uses the sliding window methodology 
# in order to scan a DNA sequence from a FASTA file, 
# and display the melting temperature along the sequence by using a chart. 
# THe chart should have 2 signals, one for each formula.
# Note: the sliding window should have 9 positions.

"""

import math
import tkinter as tk
from tkinter import filedialog, scrolledtext
import matplotlib.pyplot as plt


Na_plus = 0.001  # Fixed sodium concentration


def melting_temperature_simple(G: int, C: int, A: int, T: int):
    """Simple formula for melting temperature."""
    return 4 * (G + C) + 2 * (A + T)


def melting_temperature_alternative(CG_perc, length):
    """Alternative formula for melting temperature."""
    return 81.5 + 16.6 * math.log10(Na_plus) + 0.41 * (CG_perc) - 600 / length


def read_fasta_file(filepath):
    """Reads a FASTA file and returns only the DNA sequence (ignoring header)."""
    with open(filepath, "r") as f:
        lines = f.readlines()
    # Skip header lines starting with '>'
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    return sequence.upper()


def sliding_window(sequence, window_size=9):
    """Yields each window substring of length window_size."""
    for i in range(len(sequence) - window_size + 1):
        yield sequence[i:i + window_size]


def compute_tm_for_sequence(sequence):
    """Computes melting temperature arrays using sliding window."""
    window_size = 9
    tm_simple_values = []
    tm_alt_values = []
    positions = []

    for i, window in enumerate(sliding_window(sequence, window_size)):
        G = window.count('G')
        C = window.count('C')
        A = window.count('A')
        T = window.count('T')

        # Simple formula
        tm1 = melting_temperature_simple(G, C, A, T)
        tm_simple_values.append(tm1)

        # Alternative formula
        CG_perc = (C + G) / window_size * 100
        tm2 = melting_temperature_alternative(CG_perc, window_size)
        tm_alt_values.append(tm2)

        positions.append(i + 1)

    return positions, tm_simple_values, tm_alt_values


def plot_tm_chart(positions, tm_simple, tm_alt):
    """Displays a chart with both Tm signals."""
    plt.figure(figsize=(10, 5))
    plt.plot(positions, tm_simple, label="Tm Simple (4(G+C)+2(A+T))", color='blue')
    plt.plot(positions, tm_alt, label="Tm Alternative 81.5+16.6(log10([Na+]))+0.41*(%GC)-600/length", color='red')
    plt.title("DNA Melting Temperature along Sequence (Sliding Window = 9)")
    plt.xlabel("Position in Sequence")
    plt.ylabel("Melting Temperature (°C)")
    plt.legend()
    plt.grid(True)
    plt.show()


def analyze_fasta():
    filepath = filedialog.askopenfilename(
        title="Select a FASTA file",
        filetypes=(("FASTA files", "*.fasta *.fa"), ("All files", "*.*"))
    )

    if filepath:
        sequence = read_fasta_file(filepath)
        text_box.delete(1.0, tk.END)
        text_box.insert(tk.END, f"FASTA Sequence:\n{sequence}\n\n")

        if len(sequence) < 9:
            text_box.insert(tk.END, "Error: Sequence too short for 9-position sliding window.")
            return

        positions, tm_simple, tm_alt = compute_tm_for_sequence(sequence)

        # Display results in text box
        text_box.insert(tk.END, f"Computed {len(tm_simple)} windows.\n")
        text_box.insert(tk.END, "Showing first few Tm values:\n")
        for i in range(min(10, len(tm_simple))):
            text_box.insert(
                tk.END,
                f"Pos {positions[i]}: Tm_simple={tm_simple[i]:.2f} °C, Tm_alt={tm_alt[i]:.2f} °C\n"
            )

        # Show chart
        plot_tm_chart(positions, tm_simple, tm_alt)


# GUI Setup
root = tk.Tk()
root.title("DNA Melting Temperature Analyzer")

btn = tk.Button(root, text="Load FASTA File", command=analyze_fasta)
btn.pack(pady=10)

text_box = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=70, height=20)
text_box.pack(padx=10, pady=10)

root.mainloop()
