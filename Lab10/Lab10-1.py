"""
Implement a software application that makes DNA patterns from DNA sequences of gene promoters. Compute the C+G% and the Kappa Index of Coincidence values from each sliding window.

To ensure that the algorithm implementation works correctly, use the following steps:

1. Use the sequence: S=“CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG”.

2. Use a sliding window of length 30b.

3. Build a function to process the CpG content. The value that the function must return is: CG = 29.27

4. Build a function to process the Index of Coincidence. The value that the function must return is: IC = 27.53 

5. Plot the pattern on a chart.

6. Calculate the center of weight of the pattern.
7. Take the center of each pattern and plot it on a second chart.

8. Take the DNA sequence of a promoter and generate a pattern. Open PromKappa and see if your pattern is the same.

"""



import matplotlib.pyplot as plt

def get_counts(seq):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'total': 0}
    for char in seq:
        upper_char = char.upper()
        if upper_char in counts:
            counts[upper_char] += 1
            counts['total'] += 1
    return counts

def calculate_cg(counts):
    if counts['total'] == 0:
        return 0.0
    return ((counts['C'] + counts['G']) / counts['total']) * 100

def calculate_ic(counts):
    if counts['total'] <= 1:
        return 0.0
    
    numerator = (counts['A'] * (counts['A'] - 1)) + \
                (counts['C'] * (counts['C'] - 1)) + \
                (counts['G'] * (counts['G'] - 1)) + \
                (counts['T'] * (counts['T'] - 1))
    
    denominator = counts['total'] * (counts['total'] - 1)
    
    return (numerator / denominator) * 100

def analyze_sequence(sequence, window_size=30):
    clean_seq = "".join([c for c in sequence if c.upper() in "ACGT"])

    global_counts = get_counts(clean_seq)
    global_cg = calculate_cg(global_counts)
    global_ic = calculate_ic(global_counts)
    
    print("-" * 30)
    print(f"Global Sequence Stats (Validation)")
    print(f"Length: {len(clean_seq)}bp")
    print(f"Global CG: {global_cg:.2f}")
    print(f"Global IC:  {global_ic:.2f}")
    print("-" * 30)

    windows = []
    x_values = [] 
    y_values = []
    
    if len(clean_seq) < window_size:
        print("Error: Sequence is shorter than window size.")
        return

    print(f"Processing Sliding Windows (Size: {window_size}bp)...")
    
    for i in range(len(clean_seq) - window_size + 1):
        sub_seq = clean_seq[i : i + window_size]
        counts = get_counts(sub_seq)
        
        cg = calculate_cg(counts)
        ic = calculate_ic(counts)
        
        windows.append({
            'id': i + 1,
            'seq': sub_seq,
            'cg': cg,
            'ic': ic
        })
        x_values.append(cg)
        y_values.append(ic)
    avg_cg = sum(x_values) / len(x_values)
    avg_ic = sum(y_values) / len(y_values)
    
    print(f"Total Windows: {len(windows)}")
    print(f"Centroid (Center of Weight): CG={avg_cg:.2f}%, IC={avg_ic:.2f}")
    
    return x_values, y_values, avg_cg, avg_ic

def plot_results(x_vals, y_vals, cx, cy):

    plt.figure(figsize=(12, 5))

    plt.subplot(1, 2, 1)
    plt.scatter(x_vals, y_vals, c='blue', alpha=0.6, label='Windows')
    plt.title('Step 5: DNA Pattern (CG% vs IC)')
    plt.xlabel('C+G %')
    plt.ylabel('Kappa IC (Index of Coincidence)')
    plt.grid(True, linestyle='--', alpha=0.7)

    plt.subplot(1, 2, 2)
    plt.scatter(x_vals, y_vals, c='gray', alpha=0.1, s=10)
    plt.scatter([cx], [cy], c='red', s=200, marker='X')
    plt.title('Step 7: Center of Weight')
    plt.xlabel('Average C+G %')
    plt.ylabel('Average IC')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    S="CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
    
    WINDOW_SIZE = 30
    
    results = analyze_sequence(S, WINDOW_SIZE)
    
    if results:
        x, y, center_x, center_y = results
        plot_results(x, y, center_x, center_y)