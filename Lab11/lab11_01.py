"""
1)
    Implement a software application that aligns two DNA sequences by using the Needleman-Wunsch algorithm:
    S1="ACCGTGAAGCCAATAC"
    S2="AGCGTGCAGCCAATAC"


"""

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def needleman_wunsch_extended(seq1, seq2, match=1, mismatch=-1, gap=0):
    n, m = len(seq1), len(seq2)
    score_matrix = np.zeros((n + 1, m + 1), dtype=int)
    

    for i in range(n + 1):
        score_matrix[i][0] = i * gap
    for j in range(m + 1):
        score_matrix[0][j] = j * gap
        

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diagonal = score_matrix[i - 1][j - 1] + match
            else:
                diagonal = score_matrix[i - 1][j - 1] + mismatch
            up = score_matrix[i - 1][j] + gap
            left = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(diagonal, up, left)
            
 
    align1, align2 = "", ""
    i, j = n, m
    path = [] 
    path.append((i, j))
    
    while i > 0 or j > 0:
        current_score = score_matrix[i][j]
        
        if i > 0 and j > 0:
            is_match = (seq1[i - 1] == seq2[j - 1])
            diag_score = match if is_match else mismatch
            if current_score == score_matrix[i - 1][j - 1] + diag_score:
                align1 += seq1[i - 1]
                align2 += seq2[j - 1]
                i -= 1
                j -= 1
                path.append((i, j))
                continue
                
        if i > 0 and current_score == score_matrix[i - 1][j] + gap:
            align1 += seq1[i - 1]
            align2 += "-"
            i -= 1
            path.append((i, j))
            continue
            
        if j > 0 and current_score == score_matrix[i][j - 1] + gap:
            align1 += "-"
            align2 += seq2[j - 1]
            j -= 1
            path.append((i, j))
            
    return align1[::-1], align2[::-1], score_matrix, path

class AlignmentApp:
    def __init__(self, root):
        self.root = root
        self.root.title("app")
        self.root.geometry("1000x600")

        self.top_frame = tk.Frame(root)
        self.top_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.left_panel = tk.Frame(self.top_frame, width=250)
        self.left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=5)
        
        self.mid_panel = tk.Frame(self.top_frame) 
        self.mid_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5)

        self.right_panel = tk.Frame(self.top_frame) 
        self.right_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5)

        self.bottom_frame = tk.LabelFrame(root, text="Show Alignment:")
        self.bottom_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=10, ipady=10)

        self._init_controls()
        self._init_plots()
        self._init_output()

    def _init_controls(self):

        lf_seq = tk.LabelFrame(self.left_panel, text="Sequences")
        lf_seq.pack(fill=tk.X, pady=5)
        
        tk.Label(lf_seq, text="Sq 1 =").grid(row=0, column=0)
        self.entry_s1 = tk.Entry(lf_seq)
        self.entry_s1.insert(0, "ACCGTGAAGCCAATAC")
        self.entry_s1.grid(row=0, column=1, padx=5, pady=2)

        tk.Label(lf_seq, text="Sq 2 =").grid(row=1, column=0)
        self.entry_s2 = tk.Entry(lf_seq)
        self.entry_s2.insert(0, "AGCGTGCAGCCAATAC")
        self.entry_s2.grid(row=1, column=1, padx=5, pady=2)

       
        lf_param = tk.LabelFrame(self.left_panel, text="Parameters")
        lf_param.pack(fill=tk.X, pady=5)

        tk.Label(lf_param, text="Gap =").grid(row=0, column=0, sticky="e")
        self.entry_gap = tk.Entry(lf_param, width=5)
        self.entry_gap.insert(0, "0")
        self.entry_gap.grid(row=0, column=1)

        tk.Label(lf_param, text="Mach =").grid(row=1, column=0, sticky="e")
        self.entry_match = tk.Entry(lf_param, width=5)
        self.entry_match.insert(0, "1")
        self.entry_match.grid(row=1, column=1)

        tk.Label(lf_param, text="MMach =").grid(row=2, column=0, sticky="e")
        self.entry_mmatch = tk.Entry(lf_param, width=5)
        self.entry_mmatch.insert(0, "-1")
        self.entry_mmatch.grid(row=2, column=1)

        
        self.btn_align = tk.Button(self.left_panel, text="Align", command=self.run_alignment, height=2)
        self.btn_align.pack(fill=tk.X, pady=10)


    def _init_plots(self):
       
        self.lbl_map = tk.Label(self.mid_panel, text="Graphic representation (Heatmap)")
        self.lbl_map.pack()
        self.fig, self.ax = plt.subplots(figsize=(3,3))
        self.fig.patch.set_facecolor('#F0F0F0')
        self.canvas_plot = FigureCanvasTkAgg(self.fig, master=self.mid_panel)
        self.canvas_plot.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.ax.axis('off')

        
        self.lbl_grid = tk.Label(self.right_panel, text="Traceback path (Diagonal)")
        self.lbl_grid.pack()
        self.canvas_grid = tk.Canvas(self.right_panel, bg="white", width=250, height=250)
        self.canvas_grid.pack(fill=tk.BOTH, expand=True)

    def _init_output(self):
        self.txt_out = tk.Text(self.bottom_frame, height=8, font=("Courier New", 10))
        self.txt_out.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)
        
        
        scr = tk.Scrollbar(self.bottom_frame, command=self.txt_out.yview)
        scr.pack(side=tk.RIGHT, fill=tk.Y)
        self.txt_out.config(yscrollcommand=scr.set)

    def run_alignment(self):
        
        try:
            s1 = self.entry_s1.get().upper()
            s2 = self.entry_s2.get().upper()
            gap = int(self.entry_gap.get())
            match = int(self.entry_match.get())
            mmatch = int(self.entry_mmatch.get())
        except ValueError:
            messagebox.showerror("Error", "Please enter valid integers for parameters.")
            return

        a1, a2, matrix, path = needleman_wunsch_extended(s1, s2, match, mmatch, gap)

        self.txt_out.delete(1.0, tk.END)
        
        pipes = ""
        match_count = 0
        for i in range(len(a1)):
            if a1[i] == a2[i]:
                pipes += "|"
                match_count += 1
            else:
                pipes += " "
        
        length = len(a1)
        sim_pct = (match_count / length * 100) if length > 0 else 0

        res_str = f"{a1}\n{pipes}\n{a2}\n\n"
        res_str += f"Matches = {match_count}\n"
        res_str += f"Length  = {length}\n"
        res_str += f"Similarity = {int(sim_pct)} %\n"
        res_str += f"Tracing back: M[{len(s1)},{len(s2)}]"
        
        self.txt_out.insert(tk.END, res_str)

        self.ax.clear()
        self.ax.imshow(matrix, cmap='inferno', aspect='auto')
        self.ax.axis('off')
        self.canvas_plot.draw()

        self.canvas_grid.delete("all")
        rows, cols = matrix.shape
        
        w = self.canvas_grid.winfo_width()
        h = self.canvas_grid.winfo_height()
        
        cell_w = w / cols
        cell_h = h / rows
        
        for r in range(rows):
            for c in range(cols):
                x1 = c * cell_w
                y1 = r * cell_h
                x2 = x1 + cell_w
                y2 = y1 + cell_h
                
                if (r, c) in path:
                    color = "#C0392B" 
                else:
                    color = "#FFF9C4"
                
                self.canvas_grid.create_rectangle(x1, y1, x2, y2, fill=color, outline="black")

if __name__ == "__main__":
    root = tk.Tk()
    app = AlignmentApp(root)
    root.mainloop()