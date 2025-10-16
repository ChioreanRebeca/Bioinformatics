"""
The melting temperature (Tm) is the temperature at which one-half of a particular DNA duplex will dissociate 
and become a single strand of DNA. Primer length and sequence are of critical importance in designing the parameters 
of a successful amplification. The melting temperature of a nucleic acid duplex increases both with its length, and 
with increasing GC content. A simple formula for calculation of the (Tm) is:  Tm = 4(G + C) + 2(A + T) °C

The actual Tm is influenced by the concentration of Mg2+ , K+ , and cosolvents. An alternative formula is:

Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) - 600/length

Where Na+ is the ion concentartion of the solution and can take a value of .001


Implement an application that calculates the melting temperature of a DNA sequence using one of these formulas or both.

Input = a string of DNA

Output = temperature in celsius


"""

import math

#S1 = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"

S = "TACGTGCGCGCG"

Na_plus = 0.001 #this is just a random value

def melting_temperature_simple(G: int, C: int, A: int, T: int):
    Tm = 4 * (G + C) + 2 * (A + T)

    return Tm, f"{Tm} °C"

def melting_temperature_alternative(CG_perc, length):
    Tm = (81.5 + 16.6 * (math.log(Na_plus)) + 0.41*(CG_perc) - 600/length) * -1

    return Tm, f"{Tm} °C"


def main():
    G = 0
    C = 0
    A = 0
    T =0
    for l in S:
        if l == 'A':
            A+=1
        elif l == 'C':
            C+=1
        elif l == 'G':
            G+=1
        else:
            T+=1

    CG_perc = (C+G)/(len(S)) * 100

    Tm1, string1 = melting_temperature_simple(G, C, A, T)
    print(f"By using the formula: Tm = 4 * (G + C) + 2 * (A + T) we obtain Tm= {string1}")

    Tm2, string2 = melting_temperature_alternative(CG_perc, len(S))
    print(f"By using the formula: Tm = 81.5 + 16.6(log10([Na+])) + .41*(%GC) - 600/length we obtain Tm= {string2}")

if __name__ == "__main__":
    main()