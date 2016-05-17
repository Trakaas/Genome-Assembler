#! /usr/bin/env python3.5

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
min_overlap=cfg_items['settings']['overlap_length']

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



def main(flags,cfg_items):
    flags_list=flags.strip('').split('-')

    min_overlap=cfg_items['settings']['overlap_length']

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
                    min_overlap=int(pair[1])
                    print('\n^You flagged an overlap specification: %i' % min_overlap)
                    if min_overlap > 50:
                        print('\n\t***NOTE: This overlap is very unlikely to occur. Alignment will be very rare.***')
                    
                except ValueError:
                    pair[1]=input('\nYou input an invalid value for minimum overlap. Would you like to re-input or use default?\n \
                            Type default or a number < size of average read from your data. ~-> ')
                    if pair[1]=='default':
                        min_overlap=5
                        pass
                    else:
                        try:
                            min_overlap=int(pair[1])
                            print('\n^You flagged an overlap specification: %i' % min_overlap)
                            if min_overlap > 50:
                                print('\n\t***NOTE: This overlap is very unlikely to occur. Alignment will be very rare.***')
                        except ValueError:
                            sys.exit('\nTwo incorrect inputs. Not bothering to make more error handling. Please run from the \
beginning with correct flags.')

    else:
        print('You didn\'t flag anything so default settings will be used.')
    files_to_parse=input('\nInput file names with extensions. Seperate files with commas. ~> ') # error handling later
#    files_to_parse='custom.fastq'
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
    comp_s=trim_by_complexity(trimmed_N)
# Start of assembly
    sequence=[item for item in comp_s]
    possible=True
    merges=0
    contained=[]
    iteration_cutoff=False
    read_index=0
# outline
    # take a read search sequence list for overlap
    while possible:
        read=sequence[read_index]
        start_overlap_region=read[:min_overlap] # overlap region is at the beginning of the read
        end_overlap_region=read[-min_overlap:] # new overlap region is at the end of the read
#        print('Current index: ',read_index)
        if merges > 0: # if merges have occurred check if currently selected read is contained in [-merges:]
            for other_reads in sequence[-merges:]:
                if str(read.seq) in other_reads.seq:
                    contained.append(read) # if contained add it to contained list and delete from master list
                    sequence.pop(read_index)
                    read_index = 0
            # decide what to do with this if I have time, probably used after contig alignment
                    continue

        for other_loc,other_reads in enumerate(sequence[read_index:]): # check begin-end (does the read overlap on the end of the other reads)
            if str(start_overlap_region.seq) in other_reads[-min_overlap:].seq:
                new_read=other_reads[:-min_overlap]+read
                sequence.append(new_read) # append to list delete original sequences
                sequence.pop(read_index)
                sequence.pop(read_index+other_loc)
                read_index = 0 # reset the index, since there might be new overlaps or contained reads
                merges += 1 # merges increase by 1
                break # exit back to while loop and begin again

            if str(end_overlap_region.seq) in other_reads[:min_overlap].seq:
                new_read=read+other_reads[:-min_overlap]
                sequence.append(new_read)# append to list delete original sequences
                sequence.pop(read_index)
                sequence.pop(read_index+other_loc)
                read_index = 0 # reset the index, since there might be new overlaps or contained reads
                merges += 1 # merges increase by 1
                break # exit back to while loop and begin again
        # didn't find begin-end

#            overlap_region=read[-min_overlap:] # new overlap region is at the end of the read
#        for other_reads in sequence[read_index:]: # check end-begin (does the read overlap on the beginning of the other reads)
#            if str(end_overlap_region.seq) in other_reads[:min_overlap].seq:
#                print('end-begin')
#                new_read=read+other_reads[:-min_overlap]
#                sequence.append(new_read)# append to list delete original sequences
#                sequence.pop(read_index)
#                read_index = 0 # reset the index, since there might be new overlaps or contained reads
#                merges += 1 # merges increase by 1
#                break # exit back to while loop and begin again

        read_index+=1 # special condition, no overlap using the the current index. increment up

        if read_index==len(sequence)-2: # eventually, no more reads to compare so end the whole thing and start aligning contigs
            possible=False
            
        
    # restart loop

    print(len(sequence))
    for the_seq in sequence: # cleaning up for output, sequence records are strings and not SeqRecords
        the_seq.seq=Seq(the_seq.seq,IUPAC.ambiguous_dna)
    
    SeqIO.write(sequence, data_loc+"contigs.fasta",'fasta')


    sys.exit('\n\tProgram finished. Check the data folder for output files.')







if __name__ == '__main__':
#   cfg_items=config()
    flags=''
    if len(sys.argv)>=2: flags=' '.join(sys.argv[1:]) #error handling later

    cfg_items=config()
    main(flags,cfg_items)
