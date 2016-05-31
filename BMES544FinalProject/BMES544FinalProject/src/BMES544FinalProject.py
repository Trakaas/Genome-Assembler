#! /usr/bin/env python3.5

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
import os
import textwrap
import yaml
import sys
from Bio.Alphabet import IUPAC
import numpy as np
import re
import shutil

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

def filtering_s(compact_sequences):
    quality_reads = []
    cur_per=0.0; last_per=0.0;
    for count,record in enumerate(compact_sequences):
        term_size=shutil.get_terminal_size()
        cur_per = count/len(compact_sequences)
        if last_per != cur_per:
            p_text='  {:.2%} done.  '.format(cur_per)
            print(p_text.center(term_size[0],'~'),end='\r')
            last_per = cur_per
        
        if np.mean(record.letter_annotations["phred_quality"]) >= 10:
            record.letter_annotations={}
            if (record.seq.count('N') <= len(record.seq)/5):
                if str(record.seq[0])=='N' or str(record.seq[-1:])=='N':
                    front_strip=re.sub('^N+','',str(record.seq))
                    end_strip=re.sub('N+$','',str(record.seq))
                    record.seq=Seq(end_strip,IUPAC.ambiguous_dna)

                if (record.seq.count('A') >= len(record.seq)/2) or (record.seq.count('C') >= len(record.seq)/2) or \
                    (record.seq.count('G') >= len(record.seq)/2) or (record.seq.count('T') >= len(record.seq)/2):
                    continue

                if (len(record)<(min_overlap*1.5)):
                    continue

                else:
                    seq=str(record.seq)
                    comp_red = re.sub('A{5,}','AAA',seq)
                    comp_red = re.sub('G{5,}','GGG',comp_red)
                    comp_red = re.sub('C{5,}','CCC',comp_red)
                    comp_red = re.sub('T{5,}','TTT',comp_red)
                    record.seq=Seq(comp_red,IUPAC.ambiguous_dna)
                    quality_reads.append(record)
    print('\n\nAfter compact filtering there are %i reads.' % len(quality_reads))
    return (item for item in quality_reads)

def coverage_stats(sequence_list,genome_size_estimate):
    if genome_size_estimate=='pass':
        lengths=[len(item) for item in sequence_list]
        mean_length=np.mean(lengths)
        read_count=len(sequence_list)
    else:
        lengths=[len(item) for item in sequence_list]
        mean_length=np.mean(lengths)
        read_count=len(sequence_list)
        coverage=(read_count*mean_length)/genome_size_estimate 
        print('\nCoverage: %ix coverage of genome' % int(coverage))

    print('\nMean length: %i' % int(mean_length))
    return int(mean_length)

def contig_stats(record):
    return 0
    
def gen_s():
    estimate=input('Do you know the estimated genome size for your reads(Y/N)? ')
    while estimate!='pass':
        if estimate in ['y','Yes','Y','yes']:
            estimate=input('Input size (must be int, will assume 3mil bps if number not provided): ')
            try:
                estimate=int(estimate)
                genome_size=estimate
                break
            except (ValueError,TypeError):
                print('Not an int. Assumming 3 million bp.')
                genome_size=3000000
                break
        elif estimate in ['n','N','no','No']:
            estimate=input('Do you know what kind of genome this is(Y/N)? ')
            if estimate in ['y','Yes','Y','yes']:
                estimate=input('Is it bacterial, viral, or neither? ')
                if estimate in ['bacterial','b','bac']:
                    genome_size=3000000
                    break
                elif estimate in ['viral','Viral','VIRAL','v','V','vir','Vir']:
                    genome_size=1000000
                    break
                else:
                    print('haven\'t written this section in yet. Changing to pass.')
                    print('No information known, coverage stats will be skipped. Only sample features will be displayed.')
                    estimate='pass'
                    genome_size=estimate
                    break
            elif estimate in ['n','N','no','No']:
                estimate = 'pass'
                genome_size=estimate
                print('No information known, coverage stats will be skipped. Only sample features will be displayed.')
                break
        else:
            estimate=input('Not a valid input. Enter Y/N if you have information or pass if you don\'t know.')
            if estimate=='pass':
                break
    return genome_size

def choose_gen(data_loc):
    choice=input('Do you have a reference genome to align your contigs against(Y/N)? ')
    if choice in ['y','Yes','Y','yes']:
        choice=input('What\'s the file name (ensure the genome is a single file)? ' )

        if os.path.isfile(data_loc+choice):
            gen=SeqIO.read(data_loc+choice,'fasta')
        else:
            sys.exit('No valid files. Exiting.')
    return gen

def main(flags,cfg_items):
    flags_list=flags.strip('').split('-')
    term_size=shutil.get_terminal_size()
    min_overlap=cfg_items['settings']['overlap_length']
    mode=cfg_items['settings']['al_type']

    print(textwrap.fill(start_text.strip(), 100))
    print('\nYou\'ve chosen -->', data_loc ,'<-- as your assembly library.')
    if len(flags_list)>1:
        sep_flags=[item for item in flags_list if item != '']
        for item in sep_flags:
            pair=item.strip().split(' ')
            if pair[0]=='h' or pair[0]=='hash':
                mode='h'
                print(('  You flagged hash table mode.   ').ljust(term_size[0],'<'))
                if len(pair) > 1:
                    print('\n***Format specification for mode was wrong. Please adhere to one flag for mode. Dropping extra values.***')

            if len(pair) > 2:
                print('\n***Format specification for flag was wrong. Please adhere to one flag one value. Dropping extra values.***')

            if pair[0]=='help':
                sys.exit('Help text here later. List of flags. Modifiable settings possible.')

            if pair[0]=='o':
                try:
                    min_overlap=int(pair[1])
                    print(('  You flagged an overlap specification: %i   ' % min_overlap).ljust(term_size[0],'<'))
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
                            print(('  You flagged an overlap specification: %i   ' % min_overlap).ljust(term_size[0],'<'))
                            if min_overlap > 50:
                                print('\n\t***NOTE: This overlap is very unlikely to occur. Alignment will be very rare.***')
                        except ValueError:
                            sys.exit('\nTwo incorrect inputs. Not bothering to make more error handling. Please run from the \
beginning with correct flags.')

    else:
        print('You didn\'t flag anything so default settings will be used.')

    ass_file='o_%i_contigs%s.fasta' % (min_overlap, mode)
    print(data_loc+ass_file)
    mode_2='re'
    if os.path.isfile(data_loc+ass_file):
        mode_2=input('Found an assembly for these specifications. Would you like to re-assemble or skip to statistics(re/st)? ')
    
    if mode_2 in ['re','RE','Re','rE','re-assemble','reassemble']:
    # start file input
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
        
        genome_size=gen_s()
        coverage_stats(sequence,genome_size)
        compact_seq_trimming=filtering_s(sequence)
    #    qual_r=trim_by_qual((record for record in sequence))
    #    trimmed_N=trim_by_N(qual_r)
    #    comp_s=trim_by_complexity(trimmed_N)
    # Start of assembly
        sequence=[str(item.seq) for item in compact_seq_trimming]
        avg_length=coverage_stats(sequence,genome_size)
    # check file presence for contigs
    
        total_len=len(sequence)
        possible=True
        merges=0
        skips=0
        contained=[]
        iteration_cutoff=False
        read_index=0
        last_per = 0.0
        term_size=shutil.get_terminal_size()
        text='  {:.2%} done.  '.format(last_per)
    #    print(text.center(term_size[0],'~'),end='\r')
    
        if mode == 'h':
        # making my hash table of overlaps
            ov_hashes={}
            for read in sequence:
                overlap_region=read[:min_overlap] # overlap region is at the beginning of the read
        
                if overlap_region not in ov_hashes:
                    ov_hashes[overlap_region]=[read]
                else:
                    ov_hashes[overlap_region].append(read)
            total_len=len(ov_hashes.keys())
            count=0
            for current_hash, current_sequences in ov_hashes.items(): # iterate over table and retrieve current bucket
                possible=True
                count+=1
                for current_index, current_read in enumerate(current_sequences):
                    possible = True
                   
                    while possible:
                        target_hash=current_read[-min_overlap:]
                        try:
                            target_sequences=[pos_seq for pos_seq in ov_hashes[target_hash] if pos_seq != current_read] # ensure that target seq is not current seq
                            if target_sequences == []: # if target bucket empty, next sequence in first bucket
                                term_size=shutil.get_terminal_size()
                                possible = False
                                continue
                            
                            other_read=target_sequences[len(target_sequences)-1] # get the last sequence in the hash to merge
                            new_read=current_read[:-min_overlap] + other_read # merge the current sequence and the target sequence
    
                            if current_index == len(current_sequences) and ov_hashes[current_hash][current_index-1] in new_read:
                                ov_hashes[current_hash][current_index-1] = new_read # merge sequences in place inside current bucket
                            else:
                                ov_hashes[current_hash][current_index] = new_read
    
                            try:
                                ov_hashes[target_hash].pop(ov_hashes[target_hash].index(other_read)) # delete merged sequence from target bucket
                            except (ValueError,KeyError):
                                possible=False
                                continue
    
                        except KeyError: # if the end can't be mapped to a hash move to the next sequence
                            possible = False
                        
                cur_per = count/total_len
                term_size=shutil.get_terminal_size()
                if last_per != cur_per:
                    text='  {:.2%} done.  '.format(cur_per)
                    print(text.center(term_size[0],'~'),end='\r')
                    last_per = cur_per
            sequence=[item for key,item in ov_hashes.items()]
            sequence=[item for sublist in sequence for item in sublist if sublist != [] and item != '' and len(item)>(1.5*avg_length)]
    
            sequence.sort(key=lambda s: len(s))
            records = [SeqRecord(Seq(seq, IUPAC.ambiguous_dna), str(index)) for index,seq in enumerate(sequence) ]
            del sequence
            ass_file='o_%i_contigs%s' % (min_overlap, mode)
            SeqIO.write(records, data_loc+ass_file+".fasta",'fasta')
    
    
        if mode == 'g':
            while possible:
                read=sequence[read_index]
                start_overlap_region=read[:min_overlap] # overlap region is at the beginning of the read
                end_overlap_region=read[-min_overlap:] # new overlap region is at the end of the read
    
                cur_per = (total_len-len(sequence))/total_len
                term_size=shutil.get_terminal_size()
                if last_per != cur_per:
                    text='  {:.2%} done.  '.format(cur_per)
                    print(text.center(term_size[0],'~'),end='\r')
                    last_per = cur_per
    
    #            if merges > 0: # if merges have occurred check if currently selected read is contained in [-merges:]
    #                for other_reads in sequence[-merges:]:
    #                    if read in other_reads:
    #                        contained.append(read) # if contained add it to contained list and delete from master list
    #                        sequence.pop(read_index)
    #                        read_index = 0
    #                # decide what to do with this if I have time, probably used after contig alignment
    #                        continue
        
                for other_loc,other_reads in enumerate(sequence[read_index:]): # check begin-end (does the read overlap on the end of the other reads)
                    if start_overlap_region in other_reads[-min_overlap:]:
                        new_read=other_reads[:-min_overlap]+read
                        sequence.append(new_read) # append to list delete original sequences
                        sequence.pop(read_index)
                        sequence.pop(read_index+other_loc)
                        read_index = 0 # reset the index, since there might be new overlaps or contained reads
                        merges += 1 # merges increase by 1
                        break # exit back to while loop and begin again
        
                    if end_overlap_region in other_reads[:min_overlap]:
                        new_read=read+other_reads[:-min_overlap]
                        sequence.append(new_read)# append to list delete original sequences
                        sequence.pop(read_index)
                        sequence.pop(read_index+other_loc)
                        read_index = 0 # reset the index, since there might be new overlaps or contained reads
                        merges += 1 # merges increase by 1
                        break # exit back to while loop and begin again
                # didn't find begin-end
        
                read_index+=1 # special condition, no overlap using the the current index. increment up
        
                if read_index==len(sequence)-1: # eventually, no more reads to compare so end the whole thing and start aligning contigs
                    possible=False
                    text='  Assembly done.  '.format(last_per)
                    print(text.center(term_size[0],'!'))
    
    #        for index,the_seq in enumerate(sequence): # cleaning up for output, sequence records are strings and not SeqRecords
    #            seq=Seq(the_seq,IUPAC.ambiguous_dna)
    #            print(type(seq))
    #            the_seq=SeqRecord(Seq(the_seq,IUPAC.ambiguous_dna),str(index))
            sequence.sort(key=lambda s: len(s))
            records = [SeqRecord(Seq(seq, IUPAC.ambiguous_dna), str(index)) for index,seq in enumerate(sequence) if len(seq)>(1.5*avg_length)]
            del sequence
            ass_file='o_%i_contigs%s' % (min_overlap, mode)
            SeqIO.write(records, data_loc+ass_file+".fasta",'fasta')
            # restart loop
            
    # start statistics here (N50, alignment, coverage) make sure to do comparisons to before and after filtering

        print('\n\n Calculating N50, L50, and starting global alignment of contigs.')
        
        print('There are %i contigs.' % len(records))
        coverage_stats([str(item.seq) for item in records],genome_size)
        records.sort(key=lambda s: len(str(s.seq)))
        N50=0; cur_len=0; L50=0;

        if genome_size=='pass':
            print('Without a genome size estimate, this is basically useless. Assumming 3mil BP genome.')
            genome_size=3000000

        for ind,seq in enumerate(records):
            cur_len+=len(seq)
            if cur_len>=int(genome_size/2):
                print('N50 (sequence length at 50% genome size): %i' % len(seq))
                print('L50 (number of sequences to reach 50% genome size): %i' % (ind+1))
                break
        
    # ask if want to assemble, or just get stats
    else:
        genome_size=gen_s()
        records=list(SeqIO.parse(data_loc+ass_file,'fasta',IUPAC.ambiguous_dna)) 
        print('\n\n Calculating N50, L50, and starting global alignment of contigs.')
        
        print('There are %i contigs.' % len(records))
        coverage_stats([str(item.seq) for item in records],genome_size)
        records.sort(key=lambda s: len(str(s.seq)))
        N50=0; cur_len=0; L50=0;

        if genome_size=='pass':
            print('Without a genome size estimate, this is basically useless. Assumming 3mil BP genome.')
            genome_size=3000000

        for ind,seq in enumerate(records):
            cur_len+=len(seq)
            if cur_len>=int(genome_size/2):
                print('N50 (sequence length at 1/2 genome size): %i' % len(seq))
                print('L50 (number of sequences to reach 1/2 genome size): %i' % (ind+1))
                break
        
    
    # locally align to reference genome, free end gap penalties, ambiguous dna
    score=[]
    genome_seq=choose_gen(data_loc)
    print(str(genome_seq.seq[0:10000]))
    print('Aligning. This is going to take a while. Get some coffee. Results will be written to a file.')
    alignment = pairwise2.align.localxx(str(genome_seq.seq), str(records[0].seq))
    print(format_alignment(alignment))
#    for seq in records[:5000:5]:
#        alignment = pairwise2.align.localxx(genome_seq.seq, seq.seq)
#        print(format_alignment(alignment))
    
    
    
    sys.exit('\n\tProgram finished. Check the data folder for output files.')






if __name__ == '__main__':
#   cfg_items=config()
    flags=''
    if len(sys.argv)>=2: flags=' '.join(sys.argv[1:]) #error handling later

    cfg_items=config()
    main(flags,cfg_items)
