'''
A python script to split a multisequence FASTA file
into several single sequence files
FORMAT: python split_fasta.py [file name]
EXAMPLE: python split_fasta.py ./fasta.fas
'''
import sys
from Bio import SeqIO


def parse_fasta(fname):
    '''
    The function to parse a FASTA file
    '''
    with open(fname, "rU") as fhandle:
        seqs = SeqIO.parse(fhandle, "fasta")
        for seq in seqs:
            print("processing", seq.id)
            with open(seq.id+".fas", "a") as fnew:
                SeqIO.write(seq, fnew, "fasta")
    print("done")


def main():
    '''
    The main function
    '''
    if len(sys.argv) == 2:
        fname = sys.argv[1]
        parse_fasta(fname)
    else:
        print("FORMAT: python split_fasta.py [file name]")
        print("EXAMPLE: python split_fasta.py ./fasta.fas")


if __name__ == "__main__":
    main()
