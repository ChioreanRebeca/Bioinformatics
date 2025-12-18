"""In the influenza folder you have 10 influeza genomes, adapt your application from the prevoius assignemnt in order
to scan each genome for possible motives. For each genome, make a chart that shows the signal
with most likely locations of realfunctional motives"""

import math
import os
import matplotlib.pyplot as plt


training_sequences = [
    "GAGGTAAAC", 
    "TCCGTAAGT", 
    "CAGGTTGGA", 
    "ACAGTCAGT",
    "TAGGTCATT", 
    "TAGGTACTG",
    "ATGGTAACT",
    "CAGGTATAC", 
    "TGTGTGAGT", 
    "AAGGTAAGT"  
]

MOTIF_LENGTH = len(training_sequences[0])
BACKGROUND_PROB = 0.25
FOLDER_PATH = "./influenza"



def get_count_matrix(seqs):
    counts = {'A': [0]*MOTIF_LENGTH, 'C': [0]*MOTIF_LENGTH, 'G': [0]*MOTIF_LENGTH, 'T': [0]*MOTIF_LENGTH}
    for seq in seqs:
        for i, nucleotide in enumerate(seq):
            if nucleotide in counts:
                counts[nucleotide][i] += 1
    return counts

def get_frequency_matrix(counts, total_seqs):
   
    freqs = {'A': [], 'C': [], 'G': [], 'T': []}
    for nuc in counts:
        freqs[nuc] = [count / total_seqs for count in counts[nuc]]
    return freqs

def get_log_likelihood_matrix(freqs, bg):
    pwm = {'A': [], 'C': [], 'G': [], 'T': []}
    for nuc in freqs:
        for p in freqs[nuc]:
            if p > 0:
                score = math.log(p / bg) 
            else:
                
                score = -50.0 
            pwm[nuc].append(score)
    return pwm

def scan_genome(sequence, pwm):

    scores = []
    
    for i in range(len(sequence) - MOTIF_LENGTH + 1):
        subseq = sequence[i : i + MOTIF_LENGTH]
        score = 0
        valid_window = True
        
        for pos, nuc in enumerate(subseq):
            if nuc not in pwm:
                score = -50
                valid_window = False
                break
                
            val = pwm[nuc][pos]
            if val < -40:
                score = -50 
               
            score += val
        
        scores.append(score)
    return scores



def read_genomes_from_folder(folder):

    genomes = {}
    
   
    if not os.path.exists(folder):
        print(f"Directory '{folder}' not found. Please create it and add files.")
        return genomes

    for filename in os.listdir(folder):
        filepath = os.path.join(folder, filename)
        if os.path.isfile(filepath):
            try:
                with open(filepath, 'r') as f:
                    lines = f.readlines()
                   
                    seq = "".join([line.strip() for line in lines if not line.startswith(">")]).upper()
                    if seq:
                        genomes[filename] = seq
            except Exception as e:
                print(f"Error reading {filename}: {e}")
    return genomes



def plot_signal(scores, filename):
    """
    Plots the score signal across the genome.
    """
    plt.figure(figsize=(12, 6))
    x_axis = range(len(scores))
    
    plt.plot(x_axis, scores, label='Motif Score', color='blue', linewidth=1)
    
   
    plt.axhline(y=0, color='r', linestyle='--', label='Threshold (0)')
    
    plt.title(f"Motif Signal Scan: {filename}")
    plt.xlabel("Genome Position (bp)")
    plt.ylabel("Log-Likelihood Score")
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    

    plt.show()


def main():

    print("Building Motif Model...")
    count_matrix = get_count_matrix(training_sequences)
    freq_matrix = get_frequency_matrix(count_matrix, len(training_sequences))
    pwm = get_log_likelihood_matrix(freq_matrix, BACKGROUND_PROB)
    print("Model Built.")


    print(f"Reading genomes from {FOLDER_PATH}...")
    genomes = read_genomes_from_folder(FOLDER_PATH)
    
    if not genomes:
        print("No genomes found. (Did you create the folder and add files?)")
      
        print("\n--- RUNNING ON DUMMY DATA FOR DEMO ---")
        genomes["Dummy_Influenza_1"] = "ATCG" * 50 + "CAGGTTGGA" + "ATCG" * 50 
        genomes["Dummy_Influenza_2"] = "GGGG" * 50 + "ACAGTCAGT" + "AAAA" * 50 


    for name, seq in genomes.items():
        print(f"\nProcessing {name} (Length: {len(seq)} bp)...")
        scores = scan_genome(seq, pwm)
        

        max_score = max(scores)
        max_pos = scores.index(max_score)
        
        print(f" -> Max Signal found at position {max_pos} with score {max_score:.4f}")
        
        if max_score > 0:
            print(" -> CONCLUSION: Strong candidate motif found.")
        else:
            print(" -> CONCLUSION: No strong motif found.")
            

        plot_signal(scores, name)

if __name__ == "__main__":
    main()