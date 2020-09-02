#!/usr/bin/env python
#filter_methylGBS_reads.py
##written 12-4-19
##Groves Dixon

#import modules
import argparse
from sys import exit
from Bio import SeqIO
import numpy as np
from decimal import Decimal
import re


def check_fq(fqFile, regExp):
  """Function to check that a methylGBS fastq file has a given regular
  expression in it, as expected based on the library preparation, as 
  well as checking for PCR duplicates. Returns two lists, indexed to the 
  set of reads in the fastq, each indicating True if the read should be 
  kept, or False if it failed."""
  fastqIterator = SeqIO.parse(open(fqFile), 'fastq')
  checkList = []
  duplicateList = []
  duplicateDict = {}
  sameGenomicCount = 0
  pcrDupCount = 0
  for record in fastqIterator:
      recSeq = str(record.seq)
      nnrw = recSeq[0:31]
      dupPass = True
      if regExp.match(recSeq):
          checkList.append(True)
      else:
          checkList.append(False)
      try:
        duplicateDict[nnrw] += 1
        pcrDupCount += 1
        duplicateList.append(False)
      except KeyError:
        duplicateDict[nnrw] = 1
        duplicateList.append(True)
  totRight = sum(checkList)
  totRead = len(checkList)

  #print out results for correct product based on expected sequences
  pct = round(float(totRight) / float(totRead), 3)*100
  pct = '{}%'.format(pct)
  print('{} of {} reads ({}) found in {} had expected start'.format(totRight, totRead, pct, fqFile))

  #print out results for deduplication test
  totUniq = len(duplicateDict.keys())
  uPct = '{}%'.format(round(float(totUniq) / float(totRead), 3)*100)
  dupPct = '{}%'.format(round(float(pcrDupCount) / float(totRead), 3)*100)
  print('{} of {} reads ({}) found in {} were unqiue'.format(totUniq, totRead, uPct, fqFile))
  print('{} of {} reads ({}) found in {} were duplicates'.format(pcrDupCount, totRead, dupPct, fqFile))
  return(checkList, duplicateList)

def write_out(fqFile, outName, checkList):
  i=-1
  written=0
  fastqIterator = SeqIO.parse(open(fqFile), 'fastq')
  with open(outName, 'w') as out:
    for record in fastqIterator:
      i+=1
      check=checkList[i]
      if check:
        outRec = record.format("fastq")
        out.write(outRec)
        written += 1
      else:
        continue
  print('{} total reads written from {} to {}'.format(written, fqFile, outName))



##################################
############## MAIN ##############
##################################

if __name__ == '__main__':


  #START RUN TIME CLOCK
  ##SET UP ARGUMENT PARSING
  Description = '''
  Description:
  Filter methylGBS reads for correct starting sequences and PCR duplicates.
  Correct starting sequences are based on the library preparation.
  Forward reads should all start with NNRWCC
  Reverse reads (optional) should all start with ACAC
  To run in single end mode simply leave off the -r2 argument
  PCR duplicates are assigned with the first 20 bp of the reads are identical.
  This will include the degenerate NNRW sequence in the R1 reads, 
  which should be distinct for ~63/64 non-PCR duplicates cut from same recognition site.
  The 20bp length can be altered with the -lIdent argument.
  '''

  parser = argparse.ArgumentParser(description=Description) ##create argument parser that will automatically return help texts from global variables above
  parser.add_argument('-r1', required = True, dest = 'input1', help = 'The the input fastq')
  parser.add_argument('-r2', required = False, default = False, dest = 'input2', help = 'The the input fastq')
  parser.add_argument('-o1', required = True, dest = 'output1', help = 'The the output file name')
  parser.add_argument('-o2', required = False, default = False, dest = 'output2', help = 'The the output file name')
  parser.add_argument('-reg1', required = False, default = '^....CC', dest = 'R1_regular_expression', help = 'The regular expression expected at start of R1 reads. Default = "^....CC" ')
  parser.add_argument('-reg2', required = False, default = '^ACAC', dest = 'R2_regular_expression', help = 'The regular expression expected at start of R2 reads. Default = "^ACAC"')
  parser.add_argument('-lIdent', required = False, default = 20, dest = 'length_identical', help = 'The length of starting sequence from a read that must be identical to count as a duplicate. Default=20.')


  #--- PARSE ARGUMENTS ---#
  args = parser.parse_args()
  r1File = args.input1
  r2File = args.input2
  r1Out = args.output1
  r2Out = args.output2
  lIdentical = args.length_identical
  regExpR1 = re.compile('^....CC')
  regExpR2 = re.compile('^ACAC')


  #RUN THE CHECKING FUNCTION
  print('\nRunning filter_methylGBS_reads.py...')

  #single end mode
  if not r2File:
    print("No R2 input file given. Will run in single-end mode...")
    r1CheckList, r1DuplicateList = check_fq(r1File, regExpR1)

  #or paired end mode
  else:
    print('Both R1 and R2 input files given. Will run in paired-end mode...')
    if not r2Out:
      exit('Error. Please designate an output name with -o2')
    r1CheckList, r1DuplicateList = check_fq(r1File, regExpR1)
    r2CheckList, r2DuplicateList = check_fq(r2File, regExpR2)
    #double-check that list lengths are as they should be
    if len(r1CheckList) != len(r2CheckList):
      print('ERROR. Paired end files are of different length')
      exit()

  #double-check the shape and duplicate checkLists match in length
  if len(r1CheckList) != len(r1DuplicateList):
    print('ERROR. Check and deduplicate lists do not match. There must be a bug')
    exit()

  #BUILD COMBINED CHECKLISTS FOR WHETHER TO KEEP READS

  #set for single end mode
  if not r2Out:
    combinedCheckList = r1CheckList
    combinedDupList = r1DuplicateList

  #set for paired end mode
  else:
    #keep only reads that have right shape in forward and reverse direction
    print('Checking reads have correct shape in both forward and reverse reads...')
    combinedCheckList = []
    for i in range(len(r1CheckList)):
      combinedCheckList.append( (r1CheckList[i] and r2CheckList[i]))
    totRight = sum(combinedCheckList)
    totRead = len(combinedCheckList)
    pct = round(float(totRight) / float(totRead), 3)*100
    pct = '{}%'.format(pct)
    print('{} of {} reads ({}) had correct start for both pe reads'.format(totRight, totRead, pct))
    
    #remove only reads that look like duplicates in the forward and reverse directions
    print('Checking for duplication based on both forward and reverse reads...')
    combinedDupList = []
    for i in range(len(r1DuplicateList)):
      combinedDupList.append( (r1DuplicateList[i] or r2DuplicateList[i])) #both will be False to indicate a true duplicate
    totDup = totRead - sum(combinedDupList)
    DupPct = round(float(totDup) / float(totRead), 3)*100
    DupPct = '{}%'.format(DupPct)
    print('{} of {} reads ({}) were duplicates based on both R1 and R2 read'.format(totDup, totRead, DupPct))

  #combine the shape check and duplicate check lists
  finalCheckList = []
  for i in range(len(combinedCheckList)):
    finalCheckList.append( (combinedCheckList[i] and combinedDupList[i])) #both need to be True to keep read
  totPass = sum(finalCheckList)
  totRead = len(finalCheckList)
  pctPass = round(float(totPass) / float(totRead), 3)*100
  pctPass = '{}%'.format(pctPass)
  print('{} of {} reads ({}) had correct shape and were unique'.format(totPass, totRead, pctPass))



  #WRITE OUT

  #single end mode
  if not r2Out:
    write_out(r1File, r1Out, finalCheckList)
  else:
    write_out(r1File, r1Out, finalCheckList)
    write_out(r2File, r2Out, finalCheckList)
