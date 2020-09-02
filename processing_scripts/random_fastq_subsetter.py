#!/usr/bin/env python
#random_fastq_subsetter.py
##written 12-4-19
##Groves Dixon

#import modules
import argparse
from sys import exit
from Bio import SeqIO
import numpy as np

##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


#START RUN TIME CLOCK
  ##SET UP ARGUMENT PARSING
  Description = '''
  Description:
  Pull a random set of n sequences from a fastq file
  '''

  parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
  parser.add_argument('-i', required = True, dest = 'input_file', help = 'The the input fastq')
  parser.add_argument('-n', required = True, dest = 'n_seqs', help = 'The number of random sequences to pull')
  parser.add_argument('-o', required = True, dest = 'output_file', help = 'The the output file name')

  #--- PARSE ARGUMENTS ---#
  args = parser.parse_args()
  infile = args.input_file
  nsub = int(args.n_seqs)
  outName = args.output_file  

  #first count up the seqs
  print('counting reads in file...')
  fasSeqs = SeqIO.parse(open(infile), 'fastq')
  nseqs=0
  for seq in fasSeqs:
    nseqs += 1

  #assign random subset
  print('{} total reads'.format(nseqs))
  print('Sorting...')
  keep = np.random.choice(nseqs, nsub, replace=False)
  keep.sort()
  
  #select seqs to write
  print('selecting random read set...')
  fasSeqs = SeqIO.parse(open(infile), 'fastq')
  toWrite = []
  thisSeq = 0
  i=0
  match=keep[i]
  for seq in fasSeqs:
    thisSeq += 1
    if thisSeq==match:
      print(match)
      toWrite.append(seq)
      i+=1
      match=keep[i]
  
  #write out
  print('Writing out to {}...'.format(outName))
  SeqIO.write(toWrite, outName, "fastq")
  
    
  

