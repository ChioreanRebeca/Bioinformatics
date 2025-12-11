"""
3. Formulate 3 scoring equations that are able to show the level of similarity between the two sequences in the second.
Implement each of the scoring equations in your current implementation.
"""

import numpy as np
import matplotlib.pyplot as plt
import math

def read_fasta(filename):
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        return "".join([line.strip() for line in lines if not line.startswith(">")]).upper()
    except FileNotFoundError:
        print(f"Error: {filename} not found.")
        return ""

class GenomeVisualizer:
    def __init__(self, seq1, seq2, window=100, step=50):
        self.seq1 = seq1
        self.seq2 = seq2
        self.window = window
        self.step = step

 
    def score_identity(self, s1, s2):
        length = min(len(s1), len(s2))
        if length == 0: return 0
        matches = 0
        for i in range(length):
            if s1[i] == s2[i]:
                matches += 1
        return (matches / length) * 100

    def score_weighted(self, s1, s2):
        length = min(len(s1), len(s2))
        if length == 0: return 0
        
        raw_score = 0

        purines = {'A', 'G'}
        pyrimidines = {'C', 'T'}
        
        for i in range(length):
            base1, base2 = s1[i], s2[i]
            
            if base1 == base2:
                raw_score += 1    
            elif (base1 in purines and base2 in purines) or \
                 (base1 in pyrimidines and base2 in pyrimidines):
                raw_score += 0     
            else:
                raw_score -= 1     
                

        normalized = ((raw_score + length) / (2 * length)) * 100
        return normalized


    def score_cosine(self, s1, s2):

        def get_vector(seq):
            return [seq.count('A'), seq.count('C'), seq.count('G'), seq.count('T')]
        
        v1 = get_vector(s1)
        v2 = get_vector(s2)
        

        dot_product = sum(a * b for a, b in zip(v1, v2))
        

        mag1 = math.sqrt(sum(x**2 for x in v1))
        mag2 = math.sqrt(sum(x**2 for x in v2))
        
        if mag1 == 0 or mag2 == 0: return 0
        
   
        return (dot_product / (mag1 * mag2)) * 100

    def generate_matrix(self, method='identity'):

        rows = (len(self.seq1) - self.window) // self.step + 1
        cols = (len(self.seq2) - self.window) // self.step + 1
        matrix = np.zeros((rows, cols))
        
        print(f"Calculating Matrix using method: {method.upper()}...")
        
        for r in range(rows):
            start1 = r * self.step
            sub1 = self.seq1[start1 : start1 + self.window]
            
            for c in range(cols):
                start2 = c * self.step
                sub2 = self.seq2[start2 : start2 + self.window]
          
                if method == 'identity':
                    val = self.score_identity(sub1, sub2)
                elif method == 'weighted':
                    val = self.score_weighted(sub1, sub2)
                elif method == 'cosine':
                    val = self.score_cosine(sub1, sub2)
                else:
                    val = 0
                
                matrix[r][c] = val
        return matrix

if __name__ == "__main__":

    s1 = read_fasta("covid1.fasta")
    s2 = read_fasta("A_virus_Korea.fasta")
    
    if s1 and s2:
   
        viz = GenomeVisualizer(s1, s2, window=100, step=50)
        
      
        methods = ['identity', 'weighted', 'cosine']
        
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        
        for i, method in enumerate(methods):
            matrix = viz.generate_matrix(method)
            
            
            im = axes[i].imshow(matrix, cmap='inferno', interpolation='nearest', aspect='auto')
            axes[i].set_title(f"Method: {method.capitalize()}")
            axes[i].set_xlabel("Influenza (Windowed)")
            if i == 0: axes[i].set_ylabel("Covid-19 (Windowed)")
           
            plt.colorbar(im, ax=axes[i], fraction=0.046, pad=0.04)

        plt.suptitle("Comparison of 3 Alignment Scoring Equations", fontsize=16)
        plt.tight_layout()
        plt.show()