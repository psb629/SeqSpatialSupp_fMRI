import sys
import platform
from os.path import join, exists, abspath
from os import getcwd
import numpy as np
import h5py

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

def get_dir_glm(glm):
	dir_root = get_dir_root()

	return join(dir_root,glm)

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

	return np.array(['%02d'%(ii+1) for ii in nn])

def get_list_roi():

	return np.array(['S1', 'M1', 'PMd', 'PMv', 'SMA', 'V1', 'SPLa', 'SPLp'])

def get_S_id(subj):

	return subj.replace('R','S')

def load_spm(fname):
	file = h5py.File(fname)
	SPM = file['SPM']

	return SPM

