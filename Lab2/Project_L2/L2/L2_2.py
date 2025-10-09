"""
Find in seq S only the dinucleotides and trinucleotides that exist 
without the use of bruteforce engine. In order to achieve the results
one must verify this combination starting from the begining
"""

S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
dinucleoids = {}
trinucleoids = {}

for i in range(len(S)-1):
    dinucleoids[S[i:i+2]] = dinucleoids.get(S[i:i+2], 0) + 1
# dinucleoids = {S[i:i+2]: S.count(S[i:i+2]) for i in range(len(S) - 1)} # because I python can
print(dinucleoids)


for i in range(len(S)- 2):
    trinucleoids[S[i:i+3]]= trinucleoids.get(S[i:i+3], 0) + 1
# trinucleoids = {S[i:i+3]: S.count(S[i:i+2]) for i in range(len(S) - 2)}

print(trinucleoids)