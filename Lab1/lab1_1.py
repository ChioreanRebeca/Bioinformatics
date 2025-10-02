"""
Design an application which is able to find the alphabet of a given sequence. The alphabet
means the unique symbols from which the sequence is made.

Eg: sequence  S = "ATTGCCCCGAAAT"
    find the alph. of the sequence
"""

SEQENCE_EXAMPLE = "ATTGCCCCGAAAT"
alphabet_dict = {}

def main():
    
    sequence = input("Write the sequence here: ")
    for s in sequence:
        alphabet_dict[s] = alphabet_dict.get(s, 0) + 1

    for letter in sorted(alphabet_dict.keys()):
        print(letter)



if __name__ == "__main__":
    main()