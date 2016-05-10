#!python3

from Bio import Seq, SeqIO
import os
import textwrap
import yaml
import sys
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
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
              if min(rec.letter_annotations["phred_quality"]) >= 30)
    
    print("Saved %i reads" % count)
    return quality_reads
def trim_by_N(sequences):
    return 0
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

    #Parse files into lists of sequence records
    sequence=[]
    for filenm in valid_list:
        print('Parsing '+data_loc+filenm)
        sequence.append(list(SeqIO.parse(data_loc+filenm,'fastq')))
    sequence=[item for sublist in sequence for item in sublist] 
    


    return 0







if __name__ == '__main__':
    cfg_items=config()
    flags=''
    if len(sys.argv)>=2: flags=' '.join(sys.argv[1:]) #error handling later
    main(cfg_items,flags)
