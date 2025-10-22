"""
Implement an application that converts the coding region of a gene into an amino acid sequence. 
Use the genetic code from from moodle.
"""

GENETIC_CODE_not_very_usefull = { 
    "Phe": ["UUU", "UUC"],
    "Ser": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
    "Leu": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
    "Tyr": ["UAU", "UAC"],
    "Cys": ["UGU", "UGC"],
    "Trp": ["UGG"],
    "Pro": ["CCU", "CCC", "CCA", "CCG"],
    "His": ["CAU", "CAC"],
    "Gln": ["CAA", "CAG"],
    "Arg": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Ile": ["AUU", "AUC", "AUA"],
    "Met": ["AUG"],
    "Thr": ["ACU", "ACC", "ACA", "ACG"],
    "Asn": ["AAU", "AAC"],
    "Lys": ["AAA", "AAG"],
    "Val": ["GUU", "GUC", "GUA", "GUG"],
    "Ala": ["GCU", "GCC", "GCA", "GCG"],
    "Asp": ["GAU", "GAC"],
    "Glu": ["GAA", "GAG"],
    "Gly": ["GGU", "GGC", "GGA", "GGG"],
    "STOP": ["UAA", "UAG", "UGA"]
}

genetic_code = {
    "UUU": "Phe", 
    "UUC": "Phe",
    "UUA": "Leu", 
    "UUG": "Leu",
    "CUU": "Leu", 
    "CUC": "Leu", 
    "CUA": "Leu", 
    "CUG": "Leu",
    "AUU": "Ile", 
    "AUC": "Ile", 
    "AUA": "Ile",
    "AUG": "Met",
    "GUU": "Val", 
    "GUC": "Val", 
    "GUA": "Val", 
    "GUG": "Val",
    "UCU": "Ser", 
    "UCC": "Ser", 
    "UCA": "Ser", 
    "UCG": "Ser", 
    "AGU": "Ser", 
    "AGC": "Ser",
    "CCU": "Pro", 
    "CCC": "Pro", 
    "CCA": "Pro", 
    "CCG": "Pro",
    "ACU": "Thr", 
    "ACC": "Thr", 
    "ACA": "Thr", 
    "ACG": "Thr",
    "GCU": "Ala", 
    "GCC": "Ala", 
    "GCA": "Ala", 
    "GCG": "Ala",
    "UAU": "Tyr", 
    "UAC": "Tyr",
    "CAU": "His",
    "CAC": "His",
    "CAA": "Gln", 
    "CAG": "Gln",
    "AAU": "Asn", 
    "AAC": "Asn",
    "AAA": "Lys", 
    "AAG": "Lys",
    "GAU": "Asp", 
    "GAC": "Asp",
    "GAA": "Glu", 
    "GAG": "Glu",
    "UGU": "Cys", 
    "UGC": "Cys",
    "UGG": "Trp",
    "CGU": "Arg", 
    "CGC": "Arg", 
    "CGA": "Arg", 
    "CGG": "Arg", 
    "AGA": "Arg", 
    "AGG": "Arg",
    "GGU": "Gly", 
    "GGC": "Gly", 
    "GGA": "Gly", 
    "GGG": "Gly",
    "UAA": "Stop", "UAG": "Stop", "UGA": "Stop"
}



S = "AGAAUGGAAUUUUGA"

def start_frame(S):
    for i in range(len(S) - 2):
        codon = S[i:i+3]
        if codon == "AUG": 
            return S[i:]
    return ""


def translate_rna_to_protein(mrna_seq):
    protein = []
    
    for i in range(0, len(mrna_seq), 3):
        codon = mrna_seq[i:i+3]
        if len(codon) < 3:
            break
        amino_acid = genetic_code.get(codon)
        if amino_acid == "Stop":
            break
        protein.append(amino_acid)
    
    return "-".join(protein)


def main():
    amino_acids = translate_rna_to_protein(start_frame(S))
    print("mRNA sequence:", S)
    print("Amino acid sequence:", amino_acids)


if __name__ == "__main__":
    main()