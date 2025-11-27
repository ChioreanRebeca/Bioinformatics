"""
Download 10 influenza virus variants from NCBI and 
analyze their genomes by using the application from
assigment 1.
1) Make an electrophoresys gel for each genome
2) eliminate all lines that are in common between the gel simulation
such that the differnces will be shown
3) Merge all elctrophoresys gel simulations(that show only the differences) 
in one general elctrophoresys gel.
"""


import os
import glob
import re
from collections import OrderedDict
import math
import matplotlib.pyplot as plt


FASTA_FOLDER = "influenza"   
FASTA_GLOB = os.path.join(FASTA_FOLDER, "*.fasta")

FASTA_GLOB_ALT = os.path.join(FASTA_FOLDER, "*.fa")

ENZYMES = [
    {"name": "EcoRI",  "site": "GAATTC", "cut": 1},
    {"name": "HaeIII", "site": "GGCC",   "cut": 2}, 
]

BIN_SIZE = 10


FIGSIZE = (10, 8)
LANE_WIDTH = 0.8
BAND_LINEWIDTH = 3.0
OUTPUT_IMAGE = "merged_diff_gel.png"



def read_fasta_first_sequence(path):

    seq_lines = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                
                continue
            seq_lines.append(line)
    sequence = "".join(seq_lines).upper()
    if not sequence:
        raise ValueError(f"No sequence found in {path}")
    return sequence


def find_cut_positions(sequence, site, cut_offset):

    pattern = re.compile("(?=" + re.escape(site) + ")")
    cuts = []
    for m in pattern.finditer(sequence):
        cut_pos = m.start() + cut_offset

        cut_pos = max(0, min(len(sequence), cut_pos))
        cuts.append(cut_pos)
    return cuts


def digest_sequence(sequence, enzymes):

    cuts = []
    for enz in enzymes:
        cuts.extend(find_cut_positions(sequence, enz["site"], enz["cut"]))
   
    cuts = sorted(set(cuts))
    
    fragments = []
    prev = 0
    for c in cuts:
        fragments.append(c - prev)
        prev = c
    fragments.append(len(sequence) - prev)
   
    fragments = [f for f in fragments if f > 0]
    fragments.sort(reverse=True)
    return fragments


def round_band_size(size, bin_size=BIN_SIZE):
    return int(round(size / bin_size) * bin_size)


def bands_to_rounded_set(bands, bin_size=BIN_SIZE):
    return set(round_band_size(b, bin_size) for b in bands)


def main():
   
    files = sorted(glob.glob(FASTA_GLOB) + glob.glob(FASTA_GLOB_ALT))
    if not files:
        print(f"No FASTA files found in {FASTA_FOLDER}. Looked for {FASTA_GLOB} and {FASTA_GLOB_ALT}.")
        return

    print(f"Found {len(files)} FASTA files. Processing...")

    sample_names = []
    sample_fragments = OrderedDict()   # sample -> list of fragments (raw sizes)
    sample_rounded_sets = OrderedDict()  # sample -> set of rounded band sizes

    for path in files:
        name = os.path.splitext(os.path.basename(path))[0]
        sample_names.append(name)
        seq = read_fasta_first_sequence(path)
        fragments = digest_sequence(seq, ENZYMES)
        sample_fragments[name] = fragments
        sample_rounded_sets[name] = bands_to_rounded_set(fragments, BIN_SIZE)
        print(f"Simulated gel bands for {name}: {fragments}")

    all_sets = list(sample_rounded_sets.values())
    if not all_sets:
        print("No band data found.")
        return

    common_bands = set.intersection(*all_sets)
    common_bands_sorted = sorted(common_bands, reverse=True)
    print("\nBands common to ALL samples (rounded):", common_bands_sorted)

    diff_only_per_sample = OrderedDict()
    for name in sample_names:
        rounded_set = sample_rounded_sets[name]
        diff_set = sorted(rounded_set - common_bands, reverse=True)
        raw_fragments = []
        for b in diff_set:
            found = None
            for raw in sample_fragments[name]:
                if round_band_size(raw, BIN_SIZE) == b:
                    found = raw
                    break
            if found is not None:
                raw_fragments.append(found)
            else:
                raw_fragments.append(b)
        diff_only_per_sample[name] = sorted(raw_fragments, reverse=True)
        print(f"Difference-only bands for {name}: {diff_only_per_sample[name]}")

    merged_diff_bands = sorted(set().union(*(bands_to_rounded_set(v, BIN_SIZE) for v in diff_only_per_sample.values())), reverse=True)
    print("\nMerged unique difference bands (rounded):", merged_diff_bands)
    print(f"Total unique difference bands across all samples: {len(merged_diff_bands)}")

    num_lanes = len(sample_names)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.set_facecolor("black")
    plt.title("Merged Difference-only Electrophoresis Gel", color="white", fontsize=14)

    max_len = 0
    for frags in sample_fragments.values():
        if frags:
            max_len = max(max_len, max(frags))
    if max_len == 0:
        max_len = 1  # avoid division by zero

    for i, name in enumerate(sample_names):
        lane_x = i 
        fragments = diff_only_per_sample[name]
       
        ax.text(lane_x, 1.02, name, color="white", ha="center", va="bottom", rotation=45, fontsize=8)

       
        if not fragments:
            ax.vlines(x=lane_x, ymin=0.02, ymax=0.98, colors="dimgray", linewidth=0.5, alpha=0.5)
            continue

       
        ax.vlines(x=lane_x, ymin=0.02, ymax=0.98, colors="dimgray", linewidth=0.8, alpha=0.7)

        for frag in fragments:
           
            y = 1.0 - (frag / max_len)  # in 0..1
            y = max(0.02, min(0.98, y))  # clamp to plotting area
          
            half_width = LANE_WIDTH / 2.0
            ax.hlines(y=y, xmin=lane_x - half_width, xmax=lane_x + half_width,
                      colors="white", linewidth=BAND_LINEWIDTH, alpha=0.9)
           
            ax.text(lane_x + half_width + 0.05, y, f"{int(frag)}", color="white", fontsize=7, va="center")

    ax.set_xlim(-1, num_lanes)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([])  # hide x ticks
    ax.set_yticks([])  # hide y ticks

    num_scale_ticks = 6
    tick_values = [int(max_len * i / (num_scale_ticks - 1)) for i in range(num_scale_ticks)]
    tick_ys = [1.0 - (tv / max_len) for tv in tick_values]
    for tv, ty in zip(tick_values, tick_ys):
        ax.text(-0.9, ty, f"{tv} bp", color="white", fontsize=8, ha="left", va="center")

    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()

    plt.savefig(OUTPUT_IMAGE, dpi=300, facecolor=fig.get_facecolor())
    print(f"\nSaved merged difference-only gel to '{OUTPUT_IMAGE}'")
    plt.show()


if __name__ == "__main__":
    main()