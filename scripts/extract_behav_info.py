#!/usr/bin/env python

#####################
###     Input     ###
#####################
#   Subj_id: str
#       Subject ID
#       e.g.) 'S01', ..., 'R14'
#
#####################
###     output    ###
#####################
#	Behavioural Data: A tsv file
#		Basic information to construct GLM
#

import sys
from os.path import join
from os import makedirs
import numpy as np
import pandas as pd
import warnings

if ('-s' in sys.argv):
	idx = sys.argv.index('-s')
elif ('--subject' in sys.argv):
	idx = sys.argv.index('--subject')
subj = sys.argv[idx+1]

dir_root = '/Volumes/Diedrichsen_data$/data/SeqSpatialSupp_fMRI'
dir_behav = join(dir_root,'behavDir')
dir_work = join(dir_behav,'sub-%s'%subj)

raw_data = join(dir_work,'ssh__%s.dat'%subj)

df = pd.read_csv(raw_data, delimiter='\t')

df['cue'] = df['seqType'].map({0:'Letter',1: 'Spatial'})
df['response'] = df[['response%d'%d for d in range(5)]].apply(lambda row: int(''.join(map(str, row))), axis=1)
df = df.filter(items=['BN','TN','startTime','Cue','cueP','response','MT','RT','isError','prepTime','iti'])
df.rename(columns={'startTime':'onset', 'cueP':'sequence', 'iti':'ITI'}, inplace=True)

df.to_csv(
	join(dir_work,'behav_info.tsv'),
	sep='\t', index=False
)
