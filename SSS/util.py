import platform

from os.path import join
from os import getcwd

import numpy as np
import pandas as pd

def convert_pval_to_star(pvalue):
	res = "ns"
	if pvalue <= 0.0001:
		res = "****"
	elif pvalue <= 0.001:
		res = "***"
	elif pvalue <= 0.01:
		res = "**"
	elif pvalue <= 0.05:
		res = "*"
	
	return res

def get_dir_SSS():
	dir_current = getcwd().replace('\\','/')

	tmp = dir_current.split('/')
	idx = [ii for ii, s in enumerate(tmp) if s=='SeqSpatialSupp_fMRI'][0]

	dir_SSS = '/'.join(tmp[:idx+1])
	
	return dir_SSS

def get_dir_atlas():
	dir_SSS = get_dir_SSS()

	return join(dir_SSS,'atlas/fs_LR_32k')

def get_dir_root():
	system_name = platform.system()
	if system_name == "Linux":
	    dir_root = join('/mnt/f/SeqSpatialSupp_fMRI')
	elif system_name == "Darwin":
		dir_root = join('/Volumes/Diedrichsen_data$/data/SeqSpatialSupp_fMRI')
	elif system_name == "Windows":
		dir_root = join('F:/SeqSpatialSupp_fMRI')

	return dir_root

def get_dir_behav():
	dir_root = get_dir_root()

	return join(dir_root,'behavDir')

def get_dir_anat():
	dir_root = get_dir_root()

	return join(dir_root,'anatomicals')

def get_dir_epi():
	dir_root = get_dir_root()

	return join(dir_root,'imaging_data')

def get_dir_glm(glm):
	if isinstance(glm, int):
		glm_ = 'glm_%d'%glm
	else:
		glm_ = glm
	dir_root = get_dir_root()

	return join(dir_root,glm_)

def get_dir_surf():
	dir_root = get_dir_root()

	return join(dir_root,'surfaceWB')

def get_dir_roi():
	dir_root = get_dir_root()

	return join(dir_root,'ROI')

def get_dir_result():
	dir_root = get_dir_root()

	return join(dir_root,'results')

def get_list_sn():
	nn = [1, 2, 3, 5, 6, 8, 9, 10, 11, 12, 13, 14]

	return np.array(['%02d'%ii for ii in nn])

def get_S_id(subj):

	return subj.replace('R','S')

def get_list_cue():

	return np.array(['Letter', 'Spatial'])

def get_list_seq():

	return np.array([32451, 35124, 13254, 14523])

def get_reginfo(subj, glm):
	dir_glm = get_dir_glm(glm)
	reginfo = join(dir_glm, subj, 'reginfo.tsv')
	df = pd.read_csv(reginfo, delimiter='\t')

	return df
