"""
Design an application which is able to calculate the relative frequences of each symbol from the alphabet of the sequence. 
The alphabet means the unique symbols from which the sequence is made.

Use the same seqence as before.


relative freq the percentage of the total number of symbols that a particular symbol represents.
"""



SEQENCE_EXAMPLE = "ATTGCCCCGAAAT"
alphabet_dict = {}

def main():
    
    sequence = input("Write the sequence here: ")
    for s in sequence:
        alphabet_dict[s] = alphabet_dict.get(s, 0) + 1

    total = 0

    for letter_cnt in alphabet_dict.values():
        total += letter_cnt
    
    for letter in sorted(alphabet_dict.keys()):
        print(f"{letter}: {alphabet_dict[letter]/total *100}%")



if __name__ == "__main__":
    main()