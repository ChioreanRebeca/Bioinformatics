
"""
Download from NCBI the FASTA files containing the COVID-19 genome and the influenza genome. Use AI to compare the codon frequencies
between the two. 
a) Make a chart that shows the top 10 most frequent codons for COVID-19
b) Make a chart that shows the top 10 most frequent codons for influenza
c) Compare the two results and show the most frequent codons between the two
d) Show in the output of the console top 3 amino acids for each genome.
e) What foods have less of the aminoacids found at previous points
"""
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np

genetic_code = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "AGU": "Ser", "AGC": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr", "CAU": "His", "CAC": "His",
    "CAA": "Gln", "CAG": "Gln", "AAU": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys", "GAU": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu", "UGU": "Cys", "UGC": "Cys",
    "UGG": "Trp", "CGU": "Arg", "CGC": "Arg", "CGA": "Arg",
    "CGG": "Arg", "AGA": "Arg", "AGG": "Arg", "GGU": "Gly",
    "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
    "UAA": "Stop", "UAG": "Stop", "UGA": "Stop"
}

def read_fasta(filename):
    """Read FASTA file and return the sequence as a single string (DNA → RNA)."""
    with open(filename, "r") as f:
        lines = f.readlines()
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    # Convert DNA to RNA
    return seq.replace("T", "U")

def codon_frequencies(sequence):
    """Count codon frequencies in an RNA sequence."""
    codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3) if len(sequence[i:i+3]) == 3]
    return Counter(codons)

def compare_codon_usage(file1, file2):
    """Compare codon frequencies between two FASTA files, plot top codons, and show amino acids."""
    seq1 = read_fasta(file1)
    seq2 = read_fasta(file2)

    freq1 = codon_frequencies(seq1)
    freq2 = codon_frequencies(seq2)

    total1 = sum(freq1.values())
    total2 = sum(freq2.values())

    print(f"\n{'Codon':<6} | {'COVID-19 (%)':>12} | {'Influenza (%)':>12}")
    print("-" * 38)

    all_codons = sorted(set(freq1.keys()) | set(freq2.keys()))
    for codon in all_codons:
        f1 = (freq1[codon] / total1 * 100) if codon in freq1 else 0
        f2 = (freq2[codon] / total2 * 100) if codon in freq2 else 0
        print(f"{codon:<6} | {f1:>10.3f}% | {f2:>10.3f}%")

    # --- Top 10 codons for each virus ---
    top10_covid = [codon for codon, _ in freq1.most_common(10)]
    top10_flu = [codon for codon, _ in freq2.most_common(10)]

    # Union of top codons
    top_codons_union = sorted(set(top10_covid + top10_flu))

    # Percentages for the union
    covid_perc = [(freq1[codon] / total1 * 100) if codon in freq1 else 0 for codon in top_codons_union]
    flu_perc = [(freq2[codon] / total2 * 100) if codon in freq2 else 0 for codon in top_codons_union]

    # --- Highlight most frequent codons between the two ---
    combined_freq = [(c, covid_perc[i] + flu_perc[i]) for i, c in enumerate(top_codons_union)]
    combined_freq.sort(key=lambda x: x[1], reverse=True)
    most_frequent_codons = [c for c, _ in combined_freq[:5]]
    print("\nMost frequent codons between COVID-19 and Influenza:", most_frequent_codons)

    # --- Top 3 amino acids for each genome ---
    def top_amino_acids(freq_counter, genome_name):
        aa_counter = Counter()
        for codon, count in freq_counter.items():
            aa = genetic_code.get(codon, "Unknown")
            aa_counter[aa] += count
        print(f"Top 3 amino acids in {genome_name}: {aa_counter.most_common(3)}")

    top_amino_acids(freq1, "COVID-19")
    top_amino_acids(freq2, "Influenza")

    # --- Plotting grouped bar chart ---
    x = np.arange(len(top_codons_union))
    width = 0.35

    plt.figure(figsize=(12, 6))
    plt.bar(x - width/2, covid_perc, width, label='COVID-19', color='salmon')
    plt.bar(x + width/2, flu_perc, width, label='Influenza', color='skyblue')

    plt.xlabel("Codon")
    plt.ylabel("Frequency (%)")
    plt.title("Top Codon Usage Comparison: COVID-19 vs Influenza")
    plt.xticks(x, top_codons_union)
    plt.legend()
    plt.show()


if __name__ == "__main__":
    compare_codon_usage("covid.fasta", "influenza.fasta")
    foods = "If you want foods low in the most frequent amino acids (like Leu, Val, Ile, Ser, Gly), focus on: \nFruits – apples, oranges, berries, melons\nVegetables – lettuce, cucumbers, tomatoes, carrots\nRefined grains – white rice, white bread\nStarches – potatoes, tapioca, corn\nFats and oils – olive oil, butter, coconut oil\nThese foods are carb- or fat-dominant, not protein-dominant, so they have low amino acid content overall"
    print(foods)