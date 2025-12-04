"""
On the GitHub account you'll find a version of prompt Kappa made in JavaScript (@Gagniuc)

2. Download 10 influenza genomes and 10 COVID-19 genomes. Use these FASTA files to plot their objective digital staints. 
On a 2nd type of chart, plot the center weight for each ODS. Label each of the 20 centers of weight (plotted as a small circle)
with the name and variant of the virus.


3. Enter @gagniuc on github and vote the projects that you think are importnat/useful for your younger colleagues.
"""
import matplotlib.pyplot as plt
import os
import glob


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



def read_fasta_file(filepath):
    filename = os.path.basename(filepath)
    
    label = os.path.splitext(filename)[0]
    
    sequence = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            sequence.append(line)
    
    return label, "".join(sequence)

def analyze_sequence_data(sequence, window_size=30):
    
    clean_seq = "".join([c for c in sequence if c.upper() in "ACGT"])
    
    x_values = []
    y_values = []
    
    if len(clean_seq) < window_size:
        return [], [], 0, 0

    
    for i in range(len(clean_seq) - window_size + 1):
        sub_seq = clean_seq[i : i + window_size]
        counts = get_counts(sub_seq)
        
        cg = calculate_cg(counts)
        ic = calculate_ic(counts)
        
        x_values.append(cg)
        y_values.append(ic)

    if len(x_values) > 0:
        avg_cg = sum(x_values) / len(x_values)
        avg_ic = sum(y_values) / len(y_values)
    else:
        avg_cg, avg_ic = 0, 0
        
    return x_values, y_values, avg_cg, avg_ic

def process_folders(base_folders, window_size):
    dataset = []
    
    for folder in base_folders:
        if not os.path.exists(folder):
            print(f"Warning: Folder '{folder}' not found. Skipping.")
            continue
            
        files = glob.glob(os.path.join(folder, "*.fasta")) + \
                glob.glob(os.path.join(folder, "*.fna")) + \
                glob.glob(os.path.join(folder, "*.txt"))
        
        print(f"Processing folder '{folder}': Found {len(files)} files.")
        
        for filepath in files:
            label, seq = read_fasta_file(filepath)
            
            x, y, cx, cy = analyze_sequence_data(seq, window_size)
            
            if x:
                dataset.append({
                    'label': label,
                    'group': folder,
                    'x': x,
                    'y': y,
                    'cx': cx,
                    'cy': cy
                })
                print(f"  -> Analyzed {label} ({len(x)} windows)")
                
    return dataset

def plot_batch_results(dataset):
    if not dataset:
        print("No data available to plot.")
        return

    plt.figure(figsize=(14, 12))
    
    colors = {'influenza': '#1f77b4', 'covid': '#d62728'} 
    
    plt.subplot(2, 1, 1)
    plt.title("Objective Digital Stains (DNA Patterns)")
    
    for data in dataset:
        c = colors.get(data['group'], 'gray')
        plt.scatter(data['x'], data['y'], s=1, alpha=0.05, color=c)

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=colors['influenza'], label='Influenza (Blue)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=colors['covid'], label='Covid-19 (Red)')
    ]
    plt.legend(handles=legend_elements, loc='upper right')
    plt.xlabel("C+G %")
    plt.ylabel("Kappa IC")
    plt.grid(True, linestyle='--', alpha=0.3)

    plt.subplot(2, 1, 2)
    
    
    for data in dataset:
        c = colors.get(data['group'], 'gray')
        plt.scatter(data['cx'], data['cy'], s=150, color=c, edgecolors='black', alpha=0.9, zorder=2)
        
        plt.annotate(data['label'], 
                     (data['cx'], data['cy']), 
                     xytext=(5, 5), 
                     textcoords='offset points', 
                     fontsize=8, 
                     rotation=15,
                     alpha=0.8)

    plt.xlabel("Average C+G %")
    plt.ylabel("Average IC")
    plt.grid(True, linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plt.show()


def plot_batch_results_ver2(dataset):
    if not dataset:
        print("No data available to plot.")
        return

    flu_data = [d for d in dataset if d['group'] == 'influenza']
    cov_data = [d for d in dataset if d['group'] == 'covid']
    
    print(f"\nGenerating Plots for {len(flu_data)} Influenza and {len(cov_data)} Covid-19 samples...")

   
    fig = plt.figure(figsize=(16, 12))
    
 
    color_flu = '#1f77b4'  # Blue
    color_cov = '#d62728'  # Red
    
 
    ax1 = plt.subplot(2, 2, 1)
    ax1.set_title(f"Influenza Stains ({len(flu_data)} Genomes)")
    for data in flu_data:
       
        ax1.scatter(data['x'], data['y'], s=1, alpha=0.05, color=color_flu)
    ax1.set_xlabel("C+G %")
    ax1.set_ylabel("Kappa IC")
    ax1.set_xlim(0, 100) 
    ax1.set_ylim(0, 100) 
    ax1.grid(True, linestyle='--', alpha=0.3)

    ax2 = plt.subplot(2, 2, 2)
    ax2.set_title(f"Covid-19 Stains ({len(cov_data)} Genomes)")
    for data in cov_data:
        ax2.scatter(data['x'], data['y'], s=1, alpha=0.05, color=color_cov)
    ax2.set_xlabel("C+G %")
    ax2.set_ylabel("Kappa IC")
    ax2.set_xlim(0, 100) 
    ax2.set_ylim(0, 100) 
    ax2.grid(True, linestyle='--', alpha=0.3)

   
    ax3 = plt.subplot(2, 1, 2)
    ax3.set_title("Centers of Weight (Comparison)")
    
    
    for data in flu_data:
        ax3.scatter(data['cx'], data['cy'], s=150, color=color_flu, edgecolors='black', alpha=0.9, label='Influenza')
        ax3.annotate(data['label'], (data['cx'], data['cy']), xytext=(5, 5), textcoords='offset points', fontsize=8, rotation=15, alpha=0.7)

    
    for data in cov_data:
        ax3.scatter(data['cx'], data['cy'], s=150, color=color_cov, edgecolors='black', alpha=0.9, label='Covid-19')
        ax3.annotate(data['label'], (data['cx'], data['cy']), xytext=(5, 5), textcoords='offset points', fontsize=8, rotation=15, alpha=0.7)

  
    handles, labels = ax3.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax3.legend(by_label.values(), by_label.keys(), loc='upper right')

    ax3.set_xlabel("Average C+G %")
    ax3.set_ylabel("Average IC")
    ax3.grid(True, linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    
    WINDOW_SIZE = 30
    folders_to_scan = ['influenza', 'covid']
    
    if any(os.path.exists(f) for f in folders_to_scan):
        print("Starting Batch Analysis...")
        results = process_folders(folders_to_scan, WINDOW_SIZE)
        plot_batch_results_ver2(results)
    else:
        print("Folders 'influenza' and 'covid' not found.")
        print("Running Demo on single sequence instead...\n")
        
        S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
        x, y, cx, cy = analyze_sequence_data(S, WINDOW_SIZE)
        
        demo_data = [{
            'label': 'Demo Sequence',
            'group': 'influenza',
            'x': x, 'y': y, 'cx': cx, 'cy': cy
        }]
        
        print("-" * 30)
        print(f"Validation Target CG: 29.27 | Calculated: {cx:.2f}")
        print(f"Validation Target IC: 27.53 | Calculated: {cy:.2f}")
        print("-" * 30)
        
        plot_batch_results(demo_data)