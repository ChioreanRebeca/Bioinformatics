import random
import time
import matplotlib.pyplot as plt
import os  # <-- added to extract filenames easily

# === Function to read a FASTA file ===
def read_fasta(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    # Remove header (lines starting with ">")
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    return seq


# === Function to compute GC content ===
def gc_content(sequence):
    g = sequence.count("G")
    c = sequence.count("C")
    return (g + c) / len(sequence) * 100


# === Function to create random samples ===
def sample_fragments(dna, num_samples=2000, min_len=100, max_len=150):
    seq_length = len(dna)
    seqs = []
    for _ in range(num_samples):
        length = random.randint(min_len, max_len)
        start = random.randint(0, seq_length - length)
        sequence = dna[start:start + length]
        seqs.append(sequence)
    return seqs


# === Function to find best overlap (>= 10 bases) ===
def find_best_overlap(current_seq, fragments, used_indices):
    best_overlap = 0
    best_fragment = None
    best_index = None

    for i, frag in enumerate(fragments):
        if i in used_indices:
            continue
        # Check decreasing overlap lengths
        for j in range(100, 10, -1):
            if current_seq.endswith(frag[:j]):
                if j > best_overlap:
                    best_overlap = j
                    best_fragment = frag
                    best_index = i
    return best_fragment, best_index, best_overlap


# === Function to assemble fragments ===
def assemble_fragments(seqs):
    reconstructed = seqs[0]
    used = {0}

    while True:
        next_frag, next_idx, overlap = find_best_overlap(reconstructed, seqs, used)
        if next_frag is None or overlap < 10:
            break
        reconstructed += next_frag[overlap:]
        used.add(next_idx)
    return reconstructed


# === MAIN SECTION ===

fasta_files = [
    'viruses\\bovine_stomatitis.fasta', 'viruses\\camelpox.fasta', "viruses\\canarypox.fasta",
    "viruses\\ectomelia.fasta", "viruses\\fowlpox.fasta", "viruses\\goatpox.fasta",
    "viruses\\lumpy_skin_disease.fasta", "viruses\\sheeppox.fasta", "viruses\\taterpox.fasta",
    "viruses\\yaba_monkey_tumor.fasta"
]

assembly_times = []
gc_percentages = []
virus_names = []

for fasta_path in fasta_files:
    dna = read_fasta(fasta_path)
    gc = gc_content(dna)
    gc_percentages.append(gc)

    seqs = sample_fragments(dna)

    start_time = time.time()
    assembled = assemble_fragments(seqs)
    end_time = time.time()

    elapsed_ms = (end_time - start_time) * 1000
    assembly_times.append(elapsed_ms)

    virus_name = os.path.basename(fasta_path).replace(".fasta", "")
    virus_names.append(virus_name)

    print(f"{virus_name}: GC% = {gc:.2f}, Time = {elapsed_ms:.2f} ms")

# === Plot chart ===
plt.figure(figsize=(10, 7))
plt.scatter(gc_percentages, assembly_times, color="blue", s=80)

# Add labels next to each point
for i, name in enumerate(virus_names):
    plt.text(gc_percentages[i] + 0.1, assembly_times[i] + 5, name, fontsize=9)

plt.xlabel("Overall C + G Percentage (%)")
plt.ylabel("Assembly Time (ms)")
plt.title("DNA Assembly Time vs GC Content (10 Viral Genomes)")
plt.grid(True)
plt.tight_layout()
plt.show()
