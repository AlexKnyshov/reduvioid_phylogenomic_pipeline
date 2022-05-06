'''
A script to covert between sequence formats
Relies on Biopython
'''

import sys
import glob
import os
import logging
from Bio import AlignIO
from Bio import SeqIO
from Bio.Alphabet import IUPAC, Gapped


def check_alphabet(filename, fformat):
    '''
    The function to determine the sequence alphabet.
    The code for this function is modified from
    https://stackoverflow.com/questions/41588725/estimate-
    alphabet-in-biopython-from-fasta-file/41596563#41596563
    '''
    alphabets = [IUPAC.ambiguous_dna, IUPAC.protein]
    first_record = list(SeqIO.parse(filename, fformat))[0]
    # check DNA first:
    leftover = set(str(first_record.seq).upper()) - \
        set(alphabets[0].letters) - set(["-", "?"])
    if len(leftover) == 0:
        detected = "DNA"
        logging.info("Detected alphabet: DNA")
    else:
        leftover = set(str(first_record.seq).upper()) - \
            set(alphabets[1].letters) - set(["-", "?", "X"])
        if len(leftover) == 0:
            detected = "Prot"
            logging.info("Detected alphabet: Protein")
        else:
            logging.error("Unknown alphabet, problematic symbols: %s",
                          leftover)
            sys.exit()
    return detected


def parse_files(files, option, extlen,
                inputformat, outputfolder,
                outputformat, outputext):
    '''
    Parses files, converts, writes output
    '''
    count = 0
    for absfile in files:
        infile = absfile.split("/")[-1]
        basename = infile[:(len(infile)-extlen)]
        logging.info("basename detected: %s", basename)
        logging.info("output file: " + outputfolder + "/" +
                     basename + outputext)
        with open(infile, "rU") as input_handle:
            with open(outputfolder+
                      "/"+
                      basename+
                      outputext, "w") as output_handle:
                if option == "-a":
                    if outputformat == "nexus" or inputformat == "nexus":
                        alph = check_alphabet(input_handle, inputformat)
                        input_handle.seek(0)
                        if alph == "DNA":
                            alignments = AlignIO.parse(input_handle,
                                                       inputformat,
                                                       alphabet=Gapped(IUPAC.ambiguous_dna))
                        elif alph == "Prot":
                            alignments = AlignIO.parse(input_handle,
                                                       inputformat,
                                                       alphabet=Gapped(IUPAC.protein, '-'))
                    else:
                        alignments = AlignIO.parse(input_handle, inputformat)
                    AlignIO.write(alignments, output_handle, outputformat)
                elif option == "-s":
                    if outputformat == "nexus" or inputformat == "nexus":
                        alph = check_alphabet(input_handle, inputformat)
                        input_handle.seek(0)
                        if alph == "DNA":
                            sequences = SeqIO.parse(input_handle,
                                                    inputformat,
                                                    alphabet=Gapped(IUPAC.ambiguous_dna))
                        elif alph == "Prot":
                            sequences = SeqIO.parse(input_handle,
                                                    inputformat,
                                                    alphabet=Gapped(IUPAC.protein, '-'))
                    else:
                        sequences = SeqIO.parse(input_handle, inputformat)
                    SeqIO.write(sequences, output_handle, outputformat)
                elif option == "-print" and outputformat == "fasta":
                    sequences = SeqIO.parse(input_handle, inputformat)
                    for seq in sequences:
                        print(">"+seq.id, "\n", seq.seq, file=output_handle)
                elif option == "-print" and outputformat == "phylip-relaxed":
                    alignments = AlignIO.read(input_handle, inputformat)
                    print(str(len(alignments)) + " " +
                          str(alignments.get_alignment_length()), file=output_handle)
                    for seq in alignments:
                        print(str(seq.id)+" "+str(seq.seq), file=output_handle)
        count += 1
    logging.info("Converted %i records", count)
    logging.info("Done")


def main():
    '''
    The main function
    '''
    logging.basicConfig(format='%(asctime)s | %(levelname)s: %(message)s',
                        level=logging.NOTSET)
    if len(sys.argv) == 6 or len(sys.argv) == 8:
        if len(sys.argv) == 8:
            option = sys.argv[1]
            inputfolder = sys.argv[2]
            inputformat = sys.argv[3]
            inputext = sys.argv[4]
            outputfolder = sys.argv[5]
            outputformat = sys.argv[6]
            outputext = sys.argv[7]
            files = glob.glob(inputfolder+"/*"+inputext)
            extlen = len(inputext)
            if not os.path.exists("./"+outputfolder):
                os.makedirs("./"+outputfolder)
        else:
            option = sys.argv[1]
            infilename = sys.argv[2]
            inputformat = sys.argv[3]
            outputformat = sys.argv[4]
            outputext = sys.argv[5]
            files = [infilename]
            extlen = 0
            outputfolder = "."
        parse_files(files, option, extlen,
                    inputformat, outputfolder,
                    outputformat, outputext)
    else:
        print("fconv.py script for converting aligned or unaligned " +
              "sequence files")
        print("-----------folder input------------")
        print("FORMAT: python fconv.py [option: -a, -s, -print] " +
              "[inputfolder] [inputformat] [inputext] [outputfolder] " +
              "[outputformat] [outputext]")
        print("EXAMPLE: python fconv.py -a ./fasta fasta .fas ./phylip " +
              "phylip-relaxed .phy")
        print("------------file input-------------")
        print("FORMAT: python fconv.py [option: -a, -s, -print] [inputfile] " +
              "[inputformat] [outputformat] [outputext]")
        print("EXAMPLE: python fconv.py -a ./test.fas fasta " +
              "phylip-relaxed .phy")
        print("--------general manual--------")
        print("fasta - fasta format")
        print("phylip - basic phylip with truncated names")
        print("phylip-relaxed - extended interleaved phylip (only in -a mode)")
        print("option -print in conjunction with output format set to " +
              "fasta writes fasta file directly")
        print("option -print in conjunction with output format set to " +
              "phylip-relaxed writes non-interleaved phylip-relaxed " +
              "file directly")
        print("For options -a (AlignIO) and -s (SeqIO) see format options " +
              "in the corresponding biopython module manual. " +
              "Available formats are:")
        print("SeqIO formats:", ", ".join(SeqIO._FormatToIterator.keys()))
        print("AlignIO formats:", ", ".join(AlignIO._FormatToIterator.keys()))


if __name__ == "__main__":
    main()
