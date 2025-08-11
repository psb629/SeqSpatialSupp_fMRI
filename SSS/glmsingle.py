from os.path import join
from glob import glob

import numpy as np
import pandas as pd
import h5py

from SSS import util as ut
from SSS import deal_spm

from nilearn import image
import nibabel as nb

import surfAnalysisPy as surf
# import Functional_Fusion.atlas_map as am
# import Functional_Fusion.reliability as rel

def get_dir_glmsingle(glm=None):
	dir_work = join(ut.get_dir_root(),'GLMsingle')
	if glm is not None:
		dir_work = join(dir_work,'glm_%d'%glm)

	return dir_work

def load_reginfo(subj, glm):
	dir_glm = get_dir_glmsingle(glm)
	reginfo = pd.read_csv(join(dir_glm,subj,'reginfo.tsv'), sep='\t')

	fnames = []
	for idx in reginfo.index.values:
		fname = 'beta_%04d.nii'%(idx+1)
		fnames.append(fname)
	reginfo['beta']=fnames

	return reginfo

def load_designinfo(subj, glm):
	dir_glm = get_dir_glmsingle(glm)
	designinfo = h5py.File(join(dir_glm,subj,'DESIGNINFO.mat'))

	return designinfo

def get_designSINGLE(subj, glm, run=1):
	info = load_designinfo(subj, glm)
	ref = info['designSINGLE'][run-1,0]
	Xs = info[ref][:]

	return Xs.T

def load_map_order(subj, glm, map):
	dir_surf = join(get_dir_glmsingle(glm),'surfaceWB')
	fname = join(dir_surf,subj,'%s.%s_orders.csv'%(subj,map))
	order = np.loadtxt(fname, delimiter='\t',dtype=str)

	return order

def get_y_raw(subj, glm, run=1, trim=True, as_nii=False):
	Xs = get_designSINGLE(subj, glm, run)
	T, P = Xs.shape

	dir_raw = join(ut.get_dir_root(),'imaging_data')
	img = nb.load(join(dir_raw,subj,'%s_run_%02d.nii'%(subj,run)))
	y_raw = img.get_fdata()
	
	if trim:
		y_raw = y_raw[..., -T:]

	if as_nii:
		y_raw = nb.Nifti1Image(y_raw, affine=img.affine, header=img.header)
		return y_raw
	else:
		return y_raw, img.affine, img.header

def calc_y_hat(subj, glm, run=1):

	return 1
