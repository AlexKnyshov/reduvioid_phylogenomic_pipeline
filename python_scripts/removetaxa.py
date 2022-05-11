'''
A script to filter FASTA files based on
missing data, filter in/out specific taxa,
filter one set of files based on another,
as well as rename taxa in files based
on a provided template
'''
import csv
import sys
import glob
import os
import shutil
import logging
from Bio import SeqIO


def return_len(inputseq):
    '''
    A function to obtain gapless sequence length for filtering
    '''
    strseq = str(inputseq)
    dnaiupac = set("GATCRYWSMKHBVDN-?")
    leftover = set(strseq.upper()) - dnaiupac
    if len(leftover) == 0:
        return float(len(strseq.replace("-", "")
                         .upper().replace("?", "")
                         .replace("N", "")))/len(strseq)
    else:
        return float(len(strseq.replace("-", "")
                         .upper().replace("?", "")
                         .replace("X", "")))/len(strseq)


def mkdir():
    '''
    A function to prepare the output directory
    '''
    logging.info("creating an output folder...")
    if not os.path.exists("./rmtaxaout/"):
        # creating folder if necessary
        os.makedirs("./rmtaxaout")
    else:
        # removing old files
        shutil.rmtree("./rmtaxaout/")
        # creating folder
        os.makedirs("./rmtaxaout")


def read_taxa_list(list_file):
    '''
    A function to parse taxa list for filtering
    '''
    logging.info("reading taxalist...")
    taxalist = []
    with open(list_file, "r") as list_file_handle:
        for line in list_file_handle:
            taxalist.append(line.strip())
    logging.info("read %i records", len(taxalist))
    return taxalist


def read_taxa_csv(list_file):
    '''
    A function to parse taxa csv for renaming
    '''
    logging.info("reading taxalist...")
    taxalist = {}
    with open(list_file, "r") as list_file_handle:
        reader = csv.reader(list_file_handle)
        for row in reader:
            taxalist[row[0]] = row[1]
    logging.info("read %i records", len(taxalist))
    return taxalist


def get_file_list(inputfolder):
    '''
    A function to get input file list
    '''
    files = glob.glob(inputfolder + "/*.fas")
    if len(files) == 0:
        logging.error("no fasta (*.fas) files in the directory")
        sys.exit()
    else:
        return files


def get_trim_list(trimfolder):
    '''
    A function to get filtering file list
    '''
    trimfiles = []
    for trimfile in glob.glob(trimfolder + "/*.fas"):
        trimfiles.append(os.path.basename(trimfile))
    if len(trimfiles) == 0:
        logging.error("no fasta (*.fas) files " +
                      "in the filtering directory")
        sys.exit()
    else:
        return trimfiles


def processing_input(inputfolder, trimopt, taxalist=None,
                     trimfolder=None, threshold=0.0):
    '''
    A function to process input and perform filtering
    '''
    mkdir()
    logging.info("parsing the files...")
    files = get_file_list(inputfolder)
    for infile in files:
        absfilename = infile.split("/")[-1]
        prog = "working on file "+absfilename
        sys.stdout.write(prog+"\r")
        sys.stdout.flush()
        if trimopt == "-ll":
            if len(list(SeqIO.parse(infile, "fasta"))) >= threshold:
                shutil.copy2(inputfolder + "/" +
                             absfilename, "./rmtaxaout")
        elif trimopt == "-la":
            if len(SeqIO.parse(infile, "fasta")
                   .next().seq) >= threshold:
                shutil.copy2(inputfolder + "/" +
                             absfilename, "./rmtaxaout")
        else:
            with open("./rmtaxaout/" + absfilename, "w") as outputfile:
                count = 0
                if trimopt == "-m":
                    keeplist = []
                    trimfiles = get_trim_list(trimfolder)
                    if absfilename in trimfiles:
                        with open(trimfolder + "/" + absfilename) as trimhandle:
                            for trimseq in SeqIO.parse(trimhandle, "fasta"):
                                keeplist.append(trimseq.id)
                for seq in SeqIO.parse(infile, "fasta"):
                    if trimopt == "-e" and seq.id not in taxalist:
                        print(">" + seq.id + "\n" + seq.seq,
                              file=outputfile)
                        count += 1
                    elif trimopt == "-a" and seq.id in taxalist:
                        print(">" + seq.id + "\n" + seq.seq,
                              file=outputfile)
                        count += 1
                    elif trimopt == "-ar" and seq.id in taxalist:
                        print(">" + taxalist[seq.id] + "\n" + seq.seq,
                              file=outputfile)
                        count += 1
                    elif trimopt == "-r":
                        if seq.id in taxalist:
                            print(">" + taxalist[seq.id] + "\n" + seq.seq,
                                  file=outputfile)
                        else:
                            print(">" + seq.id + "\n" + seq.seq,
                                  file=outputfile)
                        count += 1
                    elif trimopt == "-l":
                        if return_len(seq.seq) > threshold:
                            print(">" + seq.id + "\n" + seq.seq,
                                  file=outputfile)
                            count += 1
                    elif trimopt == "-m":
                        if seq.id in keeplist:
                            print(">" + seq.id + "\n" + seq.seq,
                                  file=outputfile)
                            count += 1
            if count == 0:
                os.remove("./rmtaxaout/"+absfilename)
    print("\ndone")


def main():
    '''
    The main function
    '''
    logging.basicConfig(format='%(asctime)s | %(levelname)s: %(message)s',
                        level=logging.NOTSET)
    if len(sys.argv) == 4:
        inputfolder = sys.argv[1]
        trimopt = sys.argv[2]
        if trimopt in ("-a", "-e", "-ar", "-r"):
            if trimopt in ("-a", "-e"):
                taxalist = read_taxa_list(sys.argv[3])
            else:
                taxalist = read_taxa_csv(sys.argv[3])
            processing_input(inputfolder, trimopt, taxalist=taxalist)
        elif trimopt in ("-l", "-ll", "-la"):
            threshold = float(sys.argv[3])
            processing_input(inputfolder, trimopt, threshold=threshold)
        elif trimopt == "-m":
            trimfolder = sys.argv[3]
            processing_input(inputfolder, trimopt, trimfolder=trimfolder)
        else:
            logging.error("incorrect command line parameters")
    else:
        print("FORMAT: python removeTaxa.py [folder with fasta] " +
              "[option: -a (leave specified taxa), -ar (leave and " +
              "rename specified taxa), -r (rename only specified taxa " +
              "leaving all the rest unchanged), -e (exclude specified " +
              "taxa), -l (exclude short seq taxa), -ll (exclude loci " +
              "with few taxa), -m (exclude taxa which not present in " +
              "a companion file)] [taxalist (or csv) or lenght percent " +
              "threshold or folder to companion files]")
        print("EXAMPLE: python removeTaxa.py ./fasta -l 0.75")
        print("EXAMPLE: python removeTaxa.py ./fasta -a list.lst")
        print("EXAMPLE: python removeTaxa.py ./fasta -r list.csv")
        print("EXAMPLE: python removeTaxa.py ./fasta -m ./trimmedfasta")


if __name__ == "__main__":
    main()
