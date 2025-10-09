"""
Find the percentage for all the dinucleotide and trinucleotide combinations for the sequence:
S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA".

1. Build a brute force engine to generate all dinucleotide and trinucleotide combinations.
2. For each combination, find out the percentage inside the S sequence.
3. Show the percentage for each combination in the output of your implementation.


Nucleoids or "letters," from the genetic alphabet {A, C, G, T}.

Dinucleotide Combinations
Definition: A dinucleotide is a sequence of two nucleotides.  4 **2 = 16 unique combinations

Trinucleotide Combinations
Definition: A trinucleotide is a sequence of three nucleotides. 4 **3 = 64 unique combinations
"""


S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGAGG"
gen_alphabet =  {"A", "C", "G", "T"}
dinucleoids = {}
trinucleoids = {}
total_dn = 0
total_tn = 0
dinucleoids_perc = {}
trinucleoids_perc= {}

for n1 in gen_alphabet:
    for n2 in gen_alphabet:
        d = n1+n2
        dinucleoids[d] = 0


for n1 in gen_alphabet:
    for n2 in gen_alphabet:
        for n3 in gen_alphabet:
            t = n1+n2+n3
            trinucleoids[t] = 0


for d in  dinucleoids.keys():
    if d in S:
        dinucleoids[d] = S.count(d)
        total_dn += dinucleoids[d]


for d in dinucleoids.keys():
    dinucleoids_perc[d] = round(dinucleoids[d] / total_dn * 100.0, 2)
#print(f"Number of instances of dinucleoids in sequece S: {dinucleoids}")
print(f"Number of instances of dinucleoids percentage in sequece S: {dinucleoids_perc}")

for t in  trinucleoids.keys():
    if t in S:
        trinucleoids[t] = S.count(t)
        total_tn += trinucleoids[t]

for t in trinucleoids.keys():
    trinucleoids_perc[t] = round(trinucleoids[t] / total_tn * 100.0, 2)
#print(f"Number of instances of trinucleoids in sequece S: {trinucleoids}")
print(f"Number of instances of trinucleoids percentage in sequece S: {trinucleoids_perc}")