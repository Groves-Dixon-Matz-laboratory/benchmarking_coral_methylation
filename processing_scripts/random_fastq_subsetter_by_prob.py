#!/usr/bin/env python
#random_fastq_subsetter.py
##written 12-4-19
##Groves Dixon

#import modules
import argparse
from sys import exit
from Bio import SeqIO
import numpy as np
from decimal import Decimal

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
  parser.add_argument('-pct', required = True, dest = 'probability', help = 'The probability of pulling a sequence as a percentage')
  parser.add_argument('-o', required = True, dest = 'output_file', help = 'The the output file name')

  #--- PARSE ARGUMENTS ---#
  args = parser.parse_args()
  infile = args.input_file
  pct = args.probability
  prop=float(pct) / float(100)
  outName = args.output_file  


 
  #write out random seqs
  print('Writing random {}% of sequences from {} to {}...'.format(pct, infile, outName))
  with open(outName, 'w') as out:
    fasSeqs = SeqIO.parse(open(infile), 'fastq')
    written=0
    total = 0
    for seq in fasSeqs:
      total += 1
      if np.random.binomial(n=1, p=prop):
        outRec = seq.format("fastq")
        out.write(outRec)
        written += 1
  pctWritten = float(written) / float(total) * 100
  print('{} of {} total sequences written ({}%).'.format(written, total, pctWritten))
