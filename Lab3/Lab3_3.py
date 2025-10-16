"""
#Show the minimum and maximum values over the two signals. Also, allow the user to set a threshold(a filter) 
# that is able to take into consideration only the value above the threshold. 
# These values above the threshold should be shown to the user on a second chart as horizontal bars.
# Thus, the chunks of the signal that are above the threshold are shown as a horizontal bar over the signal. 
# Whenever the signal is below the threshold, the chart should show empty space.
"""
import math
import tkinter as tk
from tkinter import filedialog, scrolledtext
import matplotlib.pyplot as plt
import numpy as np

Na_plus = 0.001  # Fixed sodium concentration


def melting_temperature_simple(G: int, C: int, A: int, T: int):
    """Simple formula for melting temperature."""
    return 4 * (G + C) + 2 * (A + T)


def read_fasta_file(filepath):
    """Reads a FASTA file and returns only the DNA sequence (ignoring header)."""
    with open(filepath, "r") as f:
        lines = f.readlines()
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    return sequence.upper()


def sliding_window(sequence, window_size=9):
    """Yields each window substring of length window_size."""
    for i in range(len(sequence) - window_size + 1):
        yield sequence[i:i + window_size]


def compute_tm_for_sequence(sequence):
    """Computes melting temperature array using sliding window."""
    window_size = 9
    tm_values = []
    positions = []

    for i, window in enumerate(sliding_window(sequence, window_size)):
        G = window.count('G')
        C = window.count('C')
        A = window.count('A')
        T = window.count('T')

        tm = melting_temperature_simple(G, C, A, T)
        tm_values.append(tm)
        positions.append(i + 1)

    return positions, tm_values


def plot_tm_chart(positions, tm_values, threshold=None):
    """Displays the signal and threshold filter visualization."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True, sharey=True)

    # === Top Chart: Full signal ===
    ax1.plot(positions, tm_values, label="Tm Simple (4(G+C)+2(A+T))", color='blue')
    ax1.set_title("DNA Melting Temperature along Sequence (Sliding Window = 9)")
    ax1.set_xlabel("Position in Sequence")
    ax1.set_ylabel("Melting Temperature (°C)")
    if threshold is not None:
        ax1.axhline(y=threshold, color='green', linestyle='--', label=f"Threshold = {threshold} °C")
    ax1.legend()
    ax1.grid(True)

    # === Bottom Chart: Highlight regions above threshold ===
    if threshold is not None:
        above = np.array(tm_values) > threshold
        start_idx = None

        for i in range(len(above)):
            if above[i] and start_idx is None:
                start_idx = positions[i]
            elif not above[i] and start_idx is not None:
                ax2.hlines(y=threshold, xmin=start_idx, xmax=positions[i - 1],
                           color='blue', linewidth=6)
                start_idx = None

        if start_idx is not None:
            ax2.hlines(y=threshold, xmin=start_idx, xmax=positions[-1],
                       color='blue', linewidth=6)

        ax2.axhline(y=threshold, color='green', linestyle='--', label=f"Threshold = {threshold} °C")

    ax2.set_title("Regions Above Threshold")
    ax2.set_xlabel("Position in Sequence")
    ax2.set_ylabel("Melting Temperature (°C)")
    ax2.legend()
    ax2.grid(True)

    # === Ensure same limits and tick divisions ===
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_ylim(ax1.get_ylim())
    ax2.set_yticks(ax1.get_yticks())
    ax2.set_xticks(ax1.get_xticks())

    plt.tight_layout()
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

        positions, tm_values = compute_tm_for_sequence(sequence)

        # Display stats
        min_tm, max_tm = min(tm_values), max(tm_values)
        text_box.insert(tk.END, f"Computed {len(tm_values)} windows.\n\n")
        text_box.insert(tk.END, f"Tm Simple: Min = {min_tm:.2f} °C, Max = {max_tm:.2f} °C\n\n")

        # Get threshold value
        try:
            threshold = float(threshold_entry.get())
        except ValueError:
            threshold = None

        if threshold is not None:
            text_box.insert(tk.END, f"Threshold set to: {threshold} °C\n\n")

        # Show chart
        plot_tm_chart(positions, tm_values, threshold)


# === GUI Setup ===
root = tk.Tk()
root.title("DNA Melting Temperature Analyzer")

frame = tk.Frame(root)
frame.pack(pady=10)

btn = tk.Button(frame, text="Load FASTA File", command=analyze_fasta)
btn.grid(row=0, column=0, padx=5)

tk.Label(frame, text="Threshold (°C):").grid(row=0, column=1)
threshold_entry = tk.Entry(frame, width=10)
threshold_entry.grid(row=0, column=2, padx=5)
threshold_entry.insert(0, "30")  # Default threshold

text_box = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=70, height=20)
text_box.pack(padx=10, pady=10)

root.mainloop()
