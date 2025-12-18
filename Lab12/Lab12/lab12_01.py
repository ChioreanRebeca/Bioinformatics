"""
A very early step in splice site recognition is exon definition, a process that is as yet poorly understood. Communication
 between the two ends of an exon is thought to be required for this step. Computational discovery of the exon-intron border 
 or the intron-exon border or the transcription factor binding sites 
(TFBS) is a challenging but important problem of bioinformatics. Implement a software application for DNA motif finding by following the steps below.

A total of 9 known motif sequences are given:



These sequences represent the exon-intron boundary.

1. make the count matrix

2. make the weight matrix

3. make the relative frequencies matrix

4. make the Log-likelihoods Matrix

5. Analize sequence S by using the Log-likelihoods Matrix:

S="CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

Calculate the score for each sliding window.

Do you have signals indicating that the S sequence contains an exon-intron border?

"""


import math

sequences = [
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


S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

motif_length = 9
num_sequences = len(sequences)
background_prob = 0.25 

def get_count_matrix(seqs):
    counts = {'A': [0]*motif_length, 'C': [0]*motif_length, 'G': [0]*motif_length, 'T': [0]*motif_length}
    
    for seq in seqs:
        for i, nucleotide in enumerate(seq):
            counts[nucleotide][i] += 1
    return counts
def get_frequency_matrix(counts, total):
    freqs = {'A': [], 'C': [], 'G': [], 'T': []}
    for nuc in counts:
        freqs[nuc] = [count / total for count in counts[nuc]]
    return freqs

def get_log_likelihood_matrix(freqs, bg):
    pwm = {'A': [], 'C': [], 'G': [], 'T': []}
    for nuc in freqs:
        for p in freqs[nuc]:
            if p > 0:
                score = math.log(p / bg) 
            else:
              
                score = -100.0 
            pwm[nuc].append(score)
    return pwm

count_matrix = get_count_matrix(sequences)
freq_matrix = get_frequency_matrix(count_matrix, num_sequences)
pwm = get_log_likelihood_matrix(freq_matrix, background_prob)


print("--- 1. Count Matrix ---")
print(f"{' ': <4} {'1': <6} {'2': <6} {'3': <6} {'4': <6} {'5': <6} {'6': <6} {'7': <6} {'8': <6} {'9': <6}")
for nuc in ['A', 'C', 'G', 'T']:
    row_str = f"{nuc}:  " + " ".join([f"{x: <6}" for x in count_matrix[nuc]])
    print(row_str)

print("\n--- 3. Relative Frequencies Matrix ---")
print(f"{' ': <4} {'1': <6} {'2': <6} {'3': <6} {'4': <6} {'5': <6} {'6': <6} {'7': <6} {'8': <6} {'9': <6}")
for nuc in ['A', 'C', 'G', 'T']:
    row_str = f"{nuc}:  " + " ".join([f"{x:.2f}  " for x in freq_matrix[nuc]])
    print(row_str)

print("\n--- 4. Log-Likelihoods Matrix (Natural Log) ---")
print(f"{' ': <4} {'1': <6} {'2': <6} {'3': <6} {'4': <6} {'5': <6} {'6': <6} {'7': <6} {'8': <6} {'9': <6}")
for nuc in ['A', 'C', 'G', 'T']:

    row_vals = []
    for x in pwm[nuc]:
        if x < -50: row_vals.append("-Inf  ")
        else: row_vals.append(f"{x:.2f}  ")
    row_str = f"{nuc}:  " + " ".join(row_vals)
    print(row_str)


print("\n--- 5. Scanning Sequence S ---")
print(f"Sequence S: {S}")
print(f"Length: {len(S)}")

results = []


for i in range(len(S) - motif_length + 1):
    subseq = S[i : i + motif_length]
    score = 0
    valid_motif = True
    
    for pos, nuc in enumerate(subseq):
        val = pwm[nuc][pos]
        if val < -50: 
            valid_motif = False
            score = -999
            break
        score += val
        
    if valid_motif:
        results.append((i, subseq, score))


results.sort(key=lambda x: x[2], reverse=True)

print(f"\nTop Scoring Windows (Threshold > 0):")
print(f"{'Pos': <5} {'Sequence': <12} {'Score': <10}")
print("-" * 30)

has_signal = False
for res in results:

    if res[2] > 0:
        print(f"{res[0]: <5} {res[1]: <12} {res[2]:.4f}")
        has_signal = True

print("-" * 30)
if has_signal:
    print("CONCLUSION: Yes, signals indicate S contains exon-intron border candidates.")
    print(f"The strongest signal is at index {results[0][0]}: {results[0][1]}")
else:
    print("CONCLUSION: No strong signals found.")