"""
1. Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology), between 1000 and 3000 nucleotides (letters).
2. Use 5 restriction enzymes (enzyme name, recognized sequence, cleavage site):

EcoRI     	
5'GAATTC
3'CTTAAG
     
5'---G     AATTC---3'
3'---CTTAA     G---5'

BamHI     	
5'GGATCC
3'CCTAGG
     
5'---G     GATCC---3'
3'---CCTAG     G---5'

HindIII     	
5'AAGCTT
3'TTCGAA
     
5'---A     AGCTT---3'
3'---TTCGA     A---5'

TaqI     	
5'TCGA
3'AGCT
     
5'---T   CGA---3'
3'---AGC   T---5'

HaeIII*     	
5'GGCC
3'CCGG
     
5'---GG  CC---3'
3'---CC  GG---5'

Input of the implementation:
1. The recognized sequence for each restriction enzyme.
2. A DNA sequence to be digested.

Output of the implementation:
1. Number of cleavages for a each restriction enzyme.
2. Cleavage positions and length of fragments.
3. A simulation of the electrophoresis gel based on the number of restriction enzymes used.
"""
import random
import re


def load_fasta(filename, min_len=1000, max_len=3000):

    with open(filename, "r") as f:
        data = f.read().splitlines()

    seq = "".join([line.strip() for line in data if not line.startswith(">")]).upper()

    if len(seq) < min_len:
        raise ValueError("Sequence too short!")
    
    target_len = random.randint(min_len, min(max_len, len(seq)))
    start = random.randint(0, len(seq) - target_len)
    end = start + target_len

    print(f"Loaded DNA length: {len(seq)}")
    print(f"Using subsequence [{start}:{end}] (length {target_len})\n")

    return seq[start:end]

def find_cuts(sequence, enzyme):
    site = enzyme["site"]
    cut_offset = enzyme["cut"]

    cuts = []

    for match in re.finditer(site, sequence):
        cut_position = match.start() + cut_offset
        cuts.append(cut_position)

    return cuts


def compute_fragments(sequence_length, cuts):
    if not cuts:
        return [sequence_length]  # no cutting

    cuts = sorted(cuts)
    fragments = []

    prev = 0
    for c in cuts:
        fragments.append(c - prev)
        prev = c

    fragments.append(sequence_length - prev)

    return sorted(fragments, reverse=True)


def simulate_gel(fragments_list, lane_labels):

    max_len = max(frag for frags in fragments_list for frag in frags)

    gel = "\n================ ELECTROPHORESIS GEL ================\n"
    gel += "(Bands lower = larger fragments)\n\n"

    width = 60  # gel width

    for name, fragments in zip(lane_labels, fragments_list):
        gel += f"Lane: {name}\n"
        for frag in sorted(fragments, reverse=True):
            band_width = max(1, int((frag / max_len) * width))
            gel += "[" + "=" * band_width + "] " + f"({frag} bp)\n"
        gel += "\n"

    return gel

import matplotlib.pyplot as plt
import numpy as np

def show_gel_window(fragments_list, lane_labels):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_facecolor("black")
    plt.title("Restriction Digest Gel Simulation", color="white")

    max_len = max(frag for frags in fragments_list for frag in frags)

    num_lanes = len(fragments_list)
    lane_width = 1

    for i, (fragments, name) in enumerate(zip(fragments_list, lane_labels)):
        x_center = i * lane_width

        for frag in fragments:
            y_pos = 1 - (frag / max_len)

            thickness = 0.015
            smear = np.random.uniform(-0.02, 0.02)

            ax.hlines(
                y=y_pos + smear,
                xmin=x_center - 0.4,
                xmax=x_center + 0.4,
                colors="white",
                linewidth=2.5
            )

        ax.text(x_center, 1.02, name, color="white", ha="center")

    ax.set_xlim(-1, num_lanes)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.show()



if __name__ == "__main__":
    enzymes = [
        {"name": "EcoRI",   "site": "GAATTC", "cut": 1},
        {"name": "BamHI",   "site": "GGATCC", "cut": 1},
        {"name": "HindIII", "site": "AAGCTT", "cut": 1},
        {"name": "TaqI",    "site": "TCGA",   "cut": 1},
        {"name": "HaeIII",  "site": "GGCC",   "cut": 2},
    ]

    sequence = load_fasta("tuberculosis.fasta")

    lane_fragments = []
    lane_names = []

    print("============ DIGEST RESULTS ============\n")

    for enzyme in enzymes:
        name = enzyme["name"]
        site = enzyme["site"]

        cuts = find_cuts(sequence, enzyme)
        fragments = compute_fragments(len(sequence), cuts)

        lane_fragments.append(fragments)
        lane_names.append(name)

        print(f"Enzyme: {name}")
        print(f"  Recognition site: {site}")
        print(f"  Number of cleavages: {len(cuts)}")
        print(f"  Cut positions: {cuts if cuts else 'None'}")
        print(f"  Fragment lengths: {fragments}\n")

    gel_output = simulate_gel(lane_fragments, lane_names)

    show_gel_window(lane_fragments, lane_names)


    print(gel_output)
