#!python3

from Bio import SeqIO
from Bio.Seq import Seq
import os
import textwrap
import yaml
import sys
from Bio.Alphabet import IUPAC
import random
import re
import numpy as np
from itertools import product

def config():
    with open("ASM_CONFIG.yml", 'r') as ymlfile:
        cfg = yaml.load(ymlfile)
    
    #for section in cfg:
    #    print(section)

    return cfg

cfg_items=config()
data_loc=cfg_items['initial']['fastq_dir']
start_text='\n'+cfg_items['initial']['start_text']
overlap=cfg_items['settings']['overlap_length']

def trim_by_qual(sequences):
    quality_reads = [record for record in sequences \
              if np.mean(record.letter_annotations["phred_quality"]) >= 20]
    for record in quality_reads:
        record.letter_annotations={}
    count = SeqIO.write(quality_reads, data_loc+"good_quality.fasta", "fasta")
    print("\n\t\tSaved %i reads after removing reads below a mean quality of 20 and stripped extraneous quality info." % count)
    return (record for record in quality_reads)


def trim_by_N(sequences):
    N_reads=[]
    for record in sequences:
        if record.seq[0]=='N' or record.seq[-1::]=='N':
            front_strip=re.sub('^N+','',str(record.seq))
            end_strip=re.sub('N+$','',str(record.seq))
            record.seq=Seq(end_strip)
        if (record.seq.count('N') <= len(record.seq)/5):
            N_reads.append(record)

    count = SeqIO.write(N_reads, data_loc+"N_reads.fasta", "fasta")
    print("\n\tSaved %i reads after filtering out N-rich reads and stripping poly-N segments from start and end of reads." % count)
    return (record for record in N_reads)

def trim_by_complexity(sequences):
    Comp=[]
    for record in sequences:
        if (record.seq.count('A') >= len(record.seq)/2) or (record.seq.count('C') >= len(record.seq)/2) or \
            (record.seq.count('G') >= len(record.seq)/2) or (record.seq.count('T') >= len(record.seq)/2):
            continue
        seq=str(record.seq)
        comp_red = re.sub('A{5,}','AAA',seq)
        comp_red = re.sub('G{5,}','GGG',comp_red)
        comp_red = re.sub('C{5,}','CCC',comp_red)
        comp_red = re.sub('T{5,}','TTT',comp_red)
        record.seq=comp_red
        Comp.append(record)
    print('\nSaved %i reads after removing homopolymers and masking simple repeats.' % len(Comp))
    return (record for record in Comp)



def main(flags):
    flags_list=flags.strip('').split('-')
    
    print(textwrap.fill(start_text.strip(), 100))
    print('\nYou\'ve chosen -->', data_loc ,'<-- as your assembly library.')
    if len(flags_list)>1:
        sep_flags=[item for item in flags_list if item != '']
        for item in sep_flags:
            pair=item.strip().split(' ')
            if len(pair) > 2:
                print('\n***Format specification for flag was wrong. Please adhere to one flag one value. Dropping extra values.***')
            if pair[0]=='help':
                sys.exit('Help text here later. List of flags. Modifiable settings possible.')
            if pair[0]=='o':
                try:
                    overlap=int(pair[1])
                    print('\n^You flagged an overlap specification: %i' % overlap)
                    if overlap > 50:
                        print('\n\t***NOTE: This overlap is very unlikely to occur. Alignment will be very rare.***')
                    
                except ValueError:
                    pair[1]=input('\nYou input an invalid value for minimum overlap. Would you like to re-input or use default?\n \
                            Type default or a number < size of average read from your data. ~-> ')
                    if pair[1]=='default':
                        overlap=5
                        pass
                    else:
                        try:
                            overlap=int(pair[1])
                            print('\n^You flagged an overlap specification: %i' % overlap)
                            if overlap > 50:
                                print('\n\t***NOTE: This overlap is very unlikely to occur. Alignment will be very rare.***')
                        except ValueError:
                            sys.exit('\nTwo incorrect inputs. Not bothering to make more error handling. Please run from the \
beginning with correct flags.')

    else:
        print('You didn\'t flag anything so default settings will be used.')
    files_to_parse=input('\nInput file names with extensions. Seperate files with commas. ~> ') # error handling later

    # Ensure files are valid data files, if not send to invalid file list and save for verbose output
    file_list=files_to_parse.split(',')
    file_list=[item.strip(' ') for item in file_list]
    valid_list=[item for item in file_list if os.path.isfile(data_loc+item)]
    invalid_list=[item for item in file_list if item not in valid_list]
    print('\nValid filenames: '+' '.join(valid_list)+
                        '\n\tInvalid filenames saved for later.')
    if len(valid_list)==0:
        sys.exit('No valid files. Exiting.')

    #Parse files into lists of sequence records
    sequence=[]
    for filenm in valid_list:
        print('Parsing '+data_loc+filenm)
        sequence.append(list(SeqIO.parse(data_loc+filenm,'fastq',IUPAC.ambiguous_dna)))
    sequence=[item for sublist in sequence for item in sublist] 
#    starting_reads = SeqIO.write(sequence, data_loc+"starting_reads.fastq", "fastq")
#    print("\nParsing done! There are %i reads in the genome. It was saved to starting_reads.fastq in the ../data/ folder." % starting_reads)
    print('\n\n\t\tParsing done! There are %i reads in the genome.' % len(sequence))
    qual_r=trim_by_qual(sequence)
    trimmed_N=trim_by_N(qual_r)
    qual_r=trim_by_complexity(trimmed_N)
# Start of assembly
    sequence=list(qual_r)
    sys.exit('\n\tProgram finished. Check the data folder for output files.')







if __name__ == '__main__':
#   cfg_items=config()
    flags=''
    if len(sys.argv)>=2: flags=' '.join(sys.argv[1:]) #error handling later
    main(flags)
