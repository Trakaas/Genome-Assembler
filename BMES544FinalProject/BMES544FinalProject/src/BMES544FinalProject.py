#!python3

from Bio import SeqIO
from Bio.Seq import Seq
import os
import textwrap
import yaml
import sys
from Bio.Alphabet import IUPAC
import random

def config():
    with open("ASM_CONFIG.yml", 'r') as ymlfile:
        cfg = yaml.load(ymlfile)
    
    #for section in cfg:
    #    print(section)

    return cfg

cfg_items=config()
data_loc=cfg_items['initial']['fastq_dir']
start_text='\n'+cfg_items['initial']['start_text']


def trim_by_qual(sequences):
    quality_reads = (record for record in sequences \
              if min(record.letter_annotations["phred_quality"]) >= 20)
    count = SeqIO.write(quality_reads, data_loc+"good_quality.fastq", "fastq")
    print("Saved %i reads" % count)
    return quality_reads


def trim_by_N(sequences):
    N_reads = (record for record in sequences \
              if (record.seq.count('N') <= len(record.seq)/5) \
              or (record.seq[:10]=='N'*10) or (record.seq[-10:]=='N'*10))
    for record in sequence:
        if record.seq[0]=='N' or record.seq[-1::]=='N':
            ind=[record.seq.find()]


    count = SeqIO.write(N_reads, data_loc+"N_reads.fastq", "fastq")
    print("Saved %i reads" % count)
    return N_reads


def trim_by_complexity(sequences):
    return 0


def main(cfg_items, flags):
    flags_list=flags.strip('').split('-')
    
#    start_text='\n'+cfg_items['initial']['start_text']
#    data_loc=cfg_items['initial']['fastq_dir']
    print(textwrap.fill(start_text.strip(), 100))
    print('\nYou\'ve chosen -->', data_loc ,'<-- as your assembly library.')
    print('You flagged -->', [item for item in flags_list if item != ''], '<-- as your CLI options.')
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
    starting_reads = SeqIO.write(sequence, data_loc+"starting_reads.fastq", "fastq")
    print("\nParsing done! There are %i reads in the genome. It was saved to starting_reads.fastq in the ../data/ folder." % starting_reads)
#   qual_r=trim_by_qual(sequence)
#    qual_r,new_count=trim_by_N(sequence)
    new_count=trim_by_N((record for record in sequence))
#   qual_r=trim_by_qual(sequence)

    return 0







if __name__ == '__main__':
    cfg_items=config()
    flags=''
    if len(sys.argv)>=2: flags=' '.join(sys.argv[1:]) #error handling later
    main(cfg_items,flags)
