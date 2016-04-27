#!/usr/bin/env python
__author__ = 'ranko'
# --------------------------------------------------------------------------------------------------
#    ToxClassify: classifies (toxin) sequences using selection of machine learning classifiers
#
#
#  workflow:
#  1) grab sequence (from fasta file)
#  2) vectorize them using <>
#  3) run R script that does actual classification, as follows:
#    3-A) loads appropriate ML models
#    3-B) loads vectors generated from input seqs
#    3-C) runs vectors against models
#    3-D) writes output file with simple stuff: sequence -> model -> classification [csv]
#
#  4) loads output file from 3)
#  5) write 'pretty' html output [ or wrap 5 in jsp ??? ]
#
# NOTES:
#  -> in html GUI, this script is started by HTML itself (hm... how?) [jsp?]
#
#
#
#
#
# --------------------------------------------------------------------------------------------------
# import os and stuff
import os
import csv
import argparse

# import bioinformatics stuff
from RBioTools import getFaHeader
from RBioTools import getUniProtID
from RBioTools import saveAllFastas
from RBioTools import loadAllFastas
from RBioTools import getFaSeq
# import ML stuff
from RBioML import biFreqVectorize
from RBioML import triBlastEnhVectorize
from RBioML import scoredToxBitsVectorize

# ------------------------------------------------------------------
# ----------------- PARSE COMMAND LINE ARGUMENTS -------------------
# ------------------------------------------------------------------

pars = argparse.ArgumentParser()
pars.add_argument('-I','--input',required=True,help='input fasta file (one or more seqs)')
pars.add_argument('-O','--output',default='_v',help='output prefix: ex: for _v, output will be _v_etb_vec.csv, _v_stb_vec.csv and _v_bif_vec.csv [def: _v]')
pars.add_argument('--triblast',default='1',help='if 1, will do triBlastEnhanced vectorization, save it as _etb_vec.csv')
pars.add_argument('--bifreq',default='1',help='if 1, will do biFreq vectorization, save it as _bif_vec.csv')
pars.add_argument('--scoredtox',default='1',help='if 1, will do scored ToxBits vectorization, save it as _stb_vec.csv')
pars.add_argument('--hmmer',default='../tools/',help='path to hmmer binaries folder')
pars.add_argument('--toxbits',default='../dbs/blocks_hmms',help='path to ToxBits')
pars.add_argument('--classifynotb',default='1',help='scoredToxBit vectorization: if not 1, will drop sequences in which no toxbits appear [def = 1]')
pars.add_argument('--blast',default='../tools/blastp',help='path to blastp')
pars.add_argument('--tbdbs',default='../dbs', help='path to triBlastEnhanced DBs')
pars.add_argument('--verbose',default='0',help='if 1, will output various diagnostic messages and status updates [def=0]')
pars.add_argument('--tmpfolder',default='',help='temporary files folder [def = '']')


args = pars.parse_args()
ver = False
if str(args.verbose) == '1':
    ver = True

allfa = loadAllFastas(args.input)
# vectorization 1: biFreq vectorizer
if args.bifreq == '1':
    vec = biFreqVectorize(allfa, verbose=ver)
    with open(args.output+'_bif_vec.csv','w') as oF:
        writ = csv.writer(oF,delimiter=',',quotechar='"')
        for v in vec:
            writ.writerow(v)

# vectorization 2: triBLASTenhanced vectorizer [codename ultimus]
if args.triblast == '1':
    vec = triBlastEnhVectorize(allfa,blastPath=args.blast, dbNeg1=args.tbdbs+'/db_neg1.fa',dbNeg2=args.tbdbs+'/db_neg2.fa',dbPos=args.tbdbs+'/db_pos.fa',verbose=ver,tmpfolder=args.tmpfolder)
    with open(args.output+'_tbe_vec.csv','w') as oF:
        writ = csv.writer(oF,delimiter=',',quotechar='"')
        for v in vec:
            writ.writerow(v)

# vectorization 3: scored Toxbits
if args.scoredtox == '1':
    vec = scoredToxBitsVectorize(allfa,toxBits=args.toxbits,hmmer=args.hmmer,classifyNoTB=args.classifynotb,verbose=ver,tmpfolder=args.tmpfolder)
    with open(args.output+'_stb_vec.csv','w') as oF:
        writ = csv.writer(oF,delimiter=',',quotechar='"')
        for v in vec:
            writ.writerow(v)
