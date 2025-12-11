"""
2)
    Download from NCBI, the influenza COVID19 genomes and align these 2 genomes 
    by using the local alignment method. 

    Note that you have to add an in-between layer solution in order to be able to align 
    the 2 genomes files.

    Hint: the alignment made step-by-step on big regions with the connection between the results

    Note that the sequence align algorithm does nto allow for the alignment of big sequences because
    the bigger the size, the bigger the scoring matrix.

    The main result should be the visualization of the similarities between the 2 genomes.

    Note: please do not use the shortcuts given by AI (use native code)
"""
import numpy as np
import matplotlib.pyplot as plt

def read_fasta(filename):
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        sequence = "".join([line.strip() for line in lines if not line.startswith(">")])
        return sequence.upper()
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return ""

class GenomeVisualizer:
    def __init__(self, seq1, seq2, window_size=100, step_size=50):
        self.seq1 = seq1
        self.seq2 = seq2
        self.window = window_size
        self.step = step_size

    def _compare_windows(self, s1_frag, s2_frag):
        length = min(len(s1_frag), len(s2_frag))
        if length == 0: return 0
        
        matches = 0
   
        for i in range(length):
            if s1_frag[i] == s2_frag[i]:
                matches += 1
                
        return (matches / length) * 100

    def generate_matrix(self):
        len1 = len(self.seq1)
        len2 = len(self.seq2)
        
        rows = (len1 - self.window) // self.step + 1
        cols = (len2 - self.window) // self.step + 1
        
        if rows <= 0 or cols <= 0:
            print("Error: Window size is larger than one of the sequences.")
            return None

        print(f"Generating Matrix: {rows} x {cols} windows...")
 
        matrix = np.zeros((rows, cols))
        
        for r in range(rows):
            start_row = r * self.step
            window_row = self.seq1[start_row : start_row + self.window]
            
     
            for c in range(cols):
                start_col = c * self.step
                window_col = self.seq2[start_col : start_col + self.window]

                score = self._compare_windows(window_row, window_col)
                matrix[r][c] = score
                
        return matrix

def plot_alignment(matrix, name1, name2):
    if matrix is None: return

    plt.figure(figsize=(10, 8))

    plt.imshow(matrix, cmap='viridis', interpolation='nearest', aspect='auto')
    
    plt.colorbar(label='Similarity Percentage (%)')
    plt.title(f"Alignment Visualization")
    plt.xlabel(f"{name2}")
    plt.ylabel(f"{name1}")
    plt.show()

if __name__ == "__main__":
    file1 = "covid1.fasta"
    #file2 = "A_virus_California.fasta"
    file2 = "A_virus_Korea.fasta"
    
    s1 = read_fasta(file1)
    s2 = read_fasta(file2)
    
    if s1 and s2:
        print(f"Loaded {file1}: {len(s1)} bp")
        print(f"Loaded {file2}: {len(s2)} bp")
        
        visualizer = GenomeVisualizer(s1, s2, window_size=50, step_size=20)
        
   
        similarity_matrix = visualizer.generate_matrix()
        
        plot_alignment(similarity_matrix, file1, file2)
    else:
        print("Could not load sequences. Please check if files exist in the folder.")