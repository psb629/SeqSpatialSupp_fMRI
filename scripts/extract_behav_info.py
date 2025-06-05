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

df['cue'] = df['seqType'].map({0:'L', 1:'S'})
df['sequence'] = df['cueP'].map({32451:1, 35124:2, 13254:3, 14523:4})
df['response'] = df[['response%d'%d for d in range(5)]].apply(lambda row: int(''.join(map(str, row))), axis=1)
df = df.filter(items=['BN','TN','startTime','PrepTime','cue','sequence','MT','RT','isError','iti','response'])
df.rename(columns={'startTime':'onset', 'PrepTime':'prepTime', 'iti':'ITI'}, inplace=True)

df.to_csv(
	join(dir_work,'behav_info.tsv'),
	sep='\t', index=False
)
