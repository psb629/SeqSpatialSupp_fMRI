from os.path import join
from glob import glob

import numpy as np
import pandas as pd
import h5py

from SSS import util as su
from SSS import deal_spm
from SSS import image as simage

from nilearn import image
import nibabel as nb

import surfAnalysisPy as surf
# import Functional_Fusion.atlas_map as am
# import Functional_Fusion.reliability as rel

def get_dir_glmsingle(glm=None):
	dir_work = join(su.get_dir_root(),'GLMsingle')
	if glm is not None:
		dir_work = join(dir_work,'glm_%d'%glm)

	return dir_work

 #def load_reginfo(subj, glm):
 #	dir_glm = get_dir_glmsingle(glm)
 #	reginfo = pd.read_csv(join(dir_glm,subj,'reginfo.tsv'), sep='\t')
 #
 #	fnames = []
 #	for idx in reginfo.index.values:
 #		fname = 'beta_%04d.nii'%(idx+1)
 #		fnames.append(fname)
 #	reginfo['beta']=fnames
 #
 #	return reginfo

def load_designinfo(subj, glm):
	dir_glm = get_dir_glmsingle(glm)
	designinfo = h5py.File(join(dir_glm,subj,'DESIGNINFO.mat'))

	return designinfo

def load_model(subj, glm, type=0):
	dir_glm = get_dir_glmsingle(glm)
	dir_work = join(dir_glm, subj)
	if (type == 0)|(type ==-4)|(type == 'a')|(type == 'A'):
		model = 'TYPEA_ONOFF.mat'
	elif (type == 1)|(type ==-3)|(type == 'b')|(type == 'B'):
		model = 'TYPEB_FITHRF.mat'
	elif (type == 2)|(type ==-2)|(type == 'c')|(type == 'C'):
		model = 'TYPEC_FITHRF_GLMDENOISE.mat'
	elif (type == 3)|(type ==-1)|(type == 'd')|(type == 'D'):
		model = 'TYPED_FITHRF_GLMDENOISE_RR.mat'

	info = h5py.File(join(dir_work,model))

	return info

def load_FIR(subj, glm):
	dir_glm = get_dir_glmsingle(glm)
	dir_work = join(dir_glm, subj)
	
	info = h5py.File(join(dir_work,'RUNWISEFIR.mat'))
	
	return info

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

def calc_y_hat(subj, glm):
	list_run = su.get_list_run()
	dir_work = join(get_dir_glmsingle(glm),subj)

	Xs = get_designSINGLE(subj=subj, glm=glm, run=1)  # (T,P)
	T, P = Xs.shape
	del Xs

	## load mask
	mask, affine, header = simage.load_mask(subj, glm)
	spatial_shape = mask.shape
	V = np.prod(spatial_shape)
	mask = mask.get_fdata()

	## load betas (V_3d, P)
	betas = np.ones((*spatial_shape, P)) * np.nan
	for ii, fname in enumerate(sorted(glob(join(dir_work,'beta_*.nii')))):
		beta = nb.load(fname).get_fdata()
		## masking
		betas[...,ii] = beta * mask

	## check the validation
	assert betas.shape[-1] == P, f'P mismatch: X has {P} cols, beta has {betas.shape[-1]} regressors'

	## Flattened in 2D for vectorized operations
	B_2d = betas.reshape(*spatial_shape,P).reshape(V,P)  # (V,P)
	del betas

	for rr, run in enumerate(list_run):
		## get y_hat
		Xs = get_designSINGLE(subj=subj, glm=glm, run=rr+1)
		Y_hat = Xs @ B_2d.T  # (T,V) = (T,P) * (P,V)

		## reshape
		yhat_4d = Y_hat.T.reshape(*spatial_shape,T)  # (V_3d,T)

		## save the y_hat
		img = nb.Nifti1Image(yhat_4d, affine=affine, header=header)
		output = join(dir_work,'%s.Yhat.%s.nii'%(subj,run))
		nb.save(img, output)

