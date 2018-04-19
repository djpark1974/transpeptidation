#!/usr/bin/env python

'''
transpeptidation.py

Tool to generate pairings of peptide fragments derived from perfect cleavge or accepting up to n cleavage failure events.

Author: Danny Park
Initiated: 30th November 2017
Inputs: peptide fasta file, residues to cleave adjacent to, max number of missed cleavage events
Outputs: peptide fasta file listing transpeptidation products

Key observation is that for increasing n, failed cleavage events must be contiguous to yield novel peptide species.
As such, if k is number of possible cleavge sites, number of peptides possible for a single input peptide scales approximately:

num_peptides = (k+1) + [(k-(i-1)) for i = 2 to n]

'''
import sys
import math
import argparse
import logging

parser = argparse.ArgumentParser(description='A primer design tool for amplicon sequencing')
parser.add_argument('--fasta', metavar='FASTA', type=str, required=True,
                    help='fasta format file containing peptides to undergo transpeptidation')
parser.add_argument('--sites', nargs='+', type=str, required=True,
                    help='Single letter amino acid cut sites - if multiple amino acids, separate by space - e.g. --sites R K')
parser.add_argument('--max_missed', type=int, required=True,
                    help='Maximum number of missed cleavage events')
parser.add_argument('--max_fragments', type=int, default=50000,
                    help='Maximum number of fragments to shuffle with eack other')
parser.add_argument('--out', type=str, required=True,
                    help='Output file name')
parser.add_argument('--log', type=str, required=False,
                    help='Logging file name')

def Fasta(x):
    fasta_dict, headers = {}, []
    for line in x:
        if fasta_dict == {} and '>' not in line:
            continue
        elif '>' in line:
            header = line.strip().split('>')[1]
            fasta_dict[header] = ''
            headers.append(header)
        else:
            fasta_dict[header] += line.strip()
    return fasta_dict, headers

def get_Peptides(name, peptide, options):
    peptide_fragments = []
    possible_cut_sites = []
    counter = 0
    for aa in peptide:
        counter += 1
        if aa in options.sites and counter != len(peptide):			# EXCLUDE C-TERMINAL CLEAVAGE POSITION
            possible_cut_sites.append(counter)
    logging.info("Possible cleavage sites: %s", possible_cut_sites)
    num_poss_sites = len(possible_cut_sites)
    if num_poss_sites == 0:							# IF NO CUT SITES, RETURN COMPLETE PEPTIDE
        peptide_fragments.append((name, peptide))
        return peptide_fragments
    else:
        for num_missed in range(options.max_missed + 1):
            if num_missed < num_poss_sites:			# DO NOT ALLOW MORE MISSED SITES THAN POSSIBLE
                logging.info("Number of missed cleavage events: %s", num_missed)
                for index in range(num_poss_sites - num_missed):
                    if (index == 0) and (index == num_poss_sites - num_missed - 1):
                        logging.info(name+'_'+str(possible_cut_sites[index + num_missed])+'N')
                        logging.info(peptide[:possible_cut_sites[index + num_missed]])
                        peptide_fragments.append((name+'_'+str(possible_cut_sites[index + num_missed])+'N', peptide[:possible_cut_sites[index + num_missed]]))
                        logging.info(name+'_'+str(possible_cut_sites[index])+'C')
                        logging.info(peptide[possible_cut_sites[index]:])
                        peptide_fragments.append((name+'_'+str(possible_cut_sites[index])+'C', peptide[possible_cut_sites[index]:]))
                    elif index == num_poss_sites - num_missed - 1:
                        logging.info(name+'_'+str(possible_cut_sites[index-1])+'_'+str(possible_cut_sites[index + num_missed]))
                        logging.info(peptide[possible_cut_sites[index-1]:possible_cut_sites[index + num_missed]])
                        peptide_fragments.append((name+'_'+str(possible_cut_sites[index-1])+'_'+str(possible_cut_sites[index + num_missed]), peptide[possible_cut_sites[index-1]:possible_cut_sites[index + num_missed]]))
                        logging.info(name+'_'+str(possible_cut_sites[index])+'C')
                        logging.info(peptide[possible_cut_sites[index]:])
                        peptide_fragments.append((name+'_'+str(possible_cut_sites[index])+'C', peptide[possible_cut_sites[index]:]))
                    elif 0 < index < num_poss_sites - num_missed - 1:
                        logging.info(name+'_'+str(possible_cut_sites[index-1])+'_'+str(possible_cut_sites[index + num_missed]))
                        logging.info(peptide[possible_cut_sites[index-1]:possible_cut_sites[index + num_missed]])
                        peptide_fragments.append((name+'_'+str(possible_cut_sites[index-1])+'_'+str(possible_cut_sites[index + num_missed]), peptide[possible_cut_sites[index-1]:possible_cut_sites[index + num_missed]]))
                    elif index == 0:
                        logging.info(name+'_'+str(possible_cut_sites[index + num_missed])+'N')
                        logging.info(peptide[:possible_cut_sites[index + num_missed]])
                        peptide_fragments.append((name+'_'+str(possible_cut_sites[index + num_missed])+'N', peptide[:possible_cut_sites[index + num_missed]]))
            elif (name, peptide) not in peptide_fragments:	# APPEND COMPLETE PEPTIDE IF TOO MANY MISSED SITES
                logging.info(name)
                logging.info(peptide)
                peptide_fragments.append((name, peptide))
    return peptide_fragments


def main():
    options = parser.parse_args()
    if options.log:
        logging.basicConfig(filename=options.log, level=logging.DEBUG, filemode='w', format='%(message)s',)
    logging.info('Command line: {}'.format(' '.join(sys.argv)))
    out_file = open(options.out, 'w')
    fasta_handle = open(options.fasta, 'r')
    fasta_dict, headers = Fasta(fasta_handle)
    all_peptides = []
    for name in headers:
        logging.info("\nPeptide name: %s \nPeptide sequence: %s", name, fasta_dict[name])
        sub_list = get_Peptides(name, fasta_dict[name], options)
        for member in sub_list:
            if member not in all_peptides:
                all_peptides.append(member)
    logging.info("Total number of fragments: %d", len(all_peptides))
    if len(all_peptides) > options.max_fragments:
        sys.exit("WARNING! The number of fragments to shuffle is {}.\nThe max threshold is set at {}.\
               \nPlease review the input file size and the max_fragments threshold before trying again.".format(len(all_peptides), options.max_fragments))
    for peptide in all_peptides:
        for other_peptide in all_peptides:
            fused_name = '>' + peptide[0] + '_' + other_peptide[0]
            fused_peptide = peptide[1] + other_peptide[1]
            out_file.write(fused_name + '\n')
            out_file.write(fused_peptide + '\n')
    fasta_handle.flush()
    fasta_handle.close()
    out_file.flush()
    out_file.close()


if __name__ == '__main__':
    main()
