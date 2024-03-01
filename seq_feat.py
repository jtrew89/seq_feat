#!/usr/bin/env python

##Sliding window that looks for areas with different genomic coponents

##load modules required for script
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import pandas as pd
import math
import os
import re

##set up arguments
parser = argparse.ArgumentParser(description=(
	'Sequence '
						)
				)
parser.add_argument('-d', '--working_directory', dest='directory', help='directory input file(s) are kept in', required=True)
parser.add_argument('-i', '--input', dest='in_filename', help='name of input file, including extension', required=True)
parser.add_argument('-t', '--input_format', dest='in_form', help= "alignment format of input file ('fasta', 'phylip', 'nexus')", required =True)
parser.add_argument('-od', '--output_directory', dest='out_directory', help='specify directory for output file if different')
parser.add_argument('-ws', '--window_size', dest='wind_size', help='size of window you wish to use for  stats',type=int, required=True)
parser.add_argument('-ss', '--step_size', dest='step', help='size of step you wish to take from the first position of the previos window to comence the next window',type=int, required=True)
parser.add_argument('-1', '--op', dest='op', help='specify operation you wish to be carried out (GC count = gc_cnt; start + end codon count = cd_cnt; if aligned, indel count = ind_cnt; if aligned, substitution count = sub_cnt)')
parser.add_argument('-o', '--output', dest='out_filename', help='name of output file, including extension', required = False) #conditional argument based on prior arguments
parser.add_argument('-a', '--alingned_input', dest='aligned', help='use if input file is aligned', action='store_true')

args = parser.parse_args()

##set working directory
os.chdir(args.directory)

##function for GC count
def GC_count(seq):

	a = seq.count('G')
	b = seq.count('C')
	c = a+b
	return(c/(len(seq)/100))
	
##Function to count start codons  and end codons (imperfect proxy for gene count)
def codon_count(seq):

	start_codon_a = 'ATG'
	start_codon_b = 'GTG'
	stop_codon_a = 'TAA'
	stop_codon_b = 'TAG'
	stop_codon_c = 'TGA'
	#custom_codon = f'{args.in_codon}'

	start_a = seq.count(start_codon_a)
	start_b = seq.count(start_codon_b)
	stop_a = seq.count(stop_codon_a)
	stop_b = seq.count(stop_codon_b)
	stop_c = seq.count(stop_codon_c)

	return start_a, start_b, stop_a, stop_b, stop_c
	
##function for indel and substitution count and location (for indel)
def indel_sub_count(seq):

	ind = 0
	ind_pos = []
	sub = 0
	## think about where to open alignment align = AlignIO.read(,'fasta')

	for i in range(seq.get_alignment_length()):
		pos = seq[:, i]
		cur_base = ''.join(set(pos)) #to check for substitutions take only unique characters in string for this position
		cur_base = re.sub('-','',cur_base) #remove possible indel characater for next step

		if '-' in pos:
			ind = ind + 1
			ind_pos.append(str(i+1))

		elif len(cur_base) > 1: #check to see if there is more than on base, if there is there has been a substitution in this spot, as the bases are not the same
			sub = sub + 1

		else:
			pass

	return ind, ind_pos, sub

##function for sliding window
def sliding_window(seq, window_size, step):

	pos_start = 1 #position start for output
	i_start = 0 #counter to move window along sequence

	if step > window_size:
		print('Step size is great than window size, you will miss information on your serquence')

	elif len(seq) >= window_size:

		#results to return from function
		gc_results = []
		ind_results = []
		sub_results = []
		codon_results = []

		for i in range(math.floor((len(seq) - window_size + step)/step)): #to carry out as many loops through windows as can fit into the length of the sequence fully (there will be a remainder). Math floor rounds down to ensure a whole number is used

			pos_end = window_size + pos_start -1 #position end for output
			gc_results.append(str(pos_start) + '-' + str(pos_end) + '_' +str(GC_count(seq[i_start:i_start+window_size])))
			start_a, start_b, stop_a, stop_b, stop_c = codon_count(seq[i_start:i_start+window_size])
			codon_out = '_'.join([str(start_a), str(start_b), str(stop_a), str(stop_b), str(stop_c)])
			codon_results.append(str(pos_start) + '-' + str(pos_end) + '_' +str(codon_out))
			pos_start = pos_start + step #add to the start position per loop the value of the step
			i_start = i_start + step #add to window start counter so it moves for based on step size

	elif len(seq) <= window_size:
		print('Window size is greater than sequence length, please enter smaller value')

	return gc_results, codon_results

##Function to return results of operations carried out on a dataframe
def return_results(gc_count,codon_results,ind_pos=None):

	results_df = pd.DataFrame() #create dataframe to add results to
	
	results_df['Loc'] = [i.split('_')[0] for i in gc_count] #add genome location data
	results_df['GC_count'] = [i.split('_')[1] for i in gc_count] #add gc_count data

	results_df['Codon_count'] = [i.split('_')[1:6] for i in codon_results] #add codon count for all 5 codonds check for in script
	
	if args.aligned: 
		results_df['INDEL_POS'] = '' #create indel column
		inde_pos = [0] * len(results_df.index) #list the length of the datadrame filled will 0s to add actual values later, if they deviate from 0
		
		for inde in results_df.index: #check if indel is in range and add position to dataframe
			loc = results_df.loc[inde]['Loc']
			for ind in ind_pos:
				if int(loc.split('-')[0]) <= int(ind) <= int(loc.split('-')[1]): #condition for if indel position falls within the current genome location
					inde_pos[inde] += 1 #if it and others do, add to the INDEL count for the window range at the same index position in the list
				else:
					pass


	else:
		pass
	
	results_df['INDEL_POS'] = inde_pos
	results_df.set_index('Loc', inplace=True) #set index to loc for posible aligned input
	return results_df
	
##variables used to hold results
#perhaprs serieses with contig ID, postition in contig and results value

##Loop through sequences to carry out requested operations on each sequence in fasta file
if args.aligned:

	for seq in AlignIO.read(args.in_filename, args.in_form):

		align_seq = re.sub('-','',str(seq.seq))
		gc_count, codon_results = sliding_window(align_seq.upper(),args.wind_size,args.step)

		print('GC count')
		print(seq.id)
		print(gc_count)

		print('codon count')
		print(seq.id)
		print(codon_results)

else:

	for dna_record in SeqIO.parse(args.in_filename, args.in_form):

		dna_seq = dna_record.seq
		#gc_count,ind_count,substi_count,cod_count = sliding_window(dna_seq)
		gc_count, codon_results = sliding_window(dna_seq.upper(),args.wind_size,args.step)

		print('GC count')
		print(dna_record.id)
		print(gc_count)

		print('codon count')
		print(dna_record.id)
		print(codon_results)

if args.aligned:
	align = AlignIO.read(args.in_filename, args.in_form)
	indel_count, ind_pos, sub_count = indel_sub_count(align)

	print('INDEL count')
	print(indel_count)

	print('INDEL position')
	print(ind_pos)

	print('Substitution count')
	print(sub_count)

	results_df = return_results(gc_count, codon_results, ind_pos)

else:
	results_df = return_results(gc_count, codon_results)


results_df.to_csv(re.sub('.'+args.in_form,'',args.in_filename) + '_feat.tsv', sep='\t')
