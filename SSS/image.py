from os.path import join
from glob import glob

import numpy as np
import pandas as pd

from SSS import util as ut
from SSS import deal_spm

from nilearn import image
import nibabel as nb

def get_list_roi():

	return np.array(['S1', 'M1', 'PMd', 'PMv', 'SMA', 'V1', 'SPLa', 'SPLp'])

def get_underlay(path_surfAnalysisPY, hemi='L'):
	hemi_ = hemi.upper()
	path_underlay = join(path_surfAnalysisPY,'standard_mesh/fs_%s/fs_LR.32k.LR.sulc.dscalar.nii'%hemi)

	return path_underlay

def get_border(path_surfAnalysisPY, hemi='L'):
	hemi_ = hemi.upper()
	path_border = join(path_surfAnalysisPY,'standard_mesh/fs_%s/fs_LR.32k.%s.border'%(hemi_,hemi_))

	return path_border

def load_summed_roi(subj, list_roi):
	"""
	Output
		ROI_img: 3D nifti data
			ROI image with sequentially assigned natural number labels for each ROI.
	"""
	dir_roi = ut.get_dir_roi()
	S_id = ut.get_S_id(subj)
	for ii, roi in enumerate(list_roi):
		fname = join(dir_roi,S_id,'ROI.L.%s.%s.nii'%(S_id,roi))
		img_deform = nb.load(fname)
		if ii == 0:
			img_roi = img_deform
		else:
			img_roi = image.math_img(
				formula="img1 + (img1==0)*img2*%d"%(ii+1),
				img1=img_roi, img2=img_deform,
				copy_header_from="img1"
			)

	return img_roi

def load_yraw(subj,roi,hemi='L'):
	"""
	Output
		y: cifti data with (# total TRs) X (# voxels) shape
	"""
	dir_roi = ut.get_dir_roi()
	dir_work = join(dir_roi,subj)

	hemi_ = hemi.upper()

	fname = join(dir_work,'cifti.%s.%s.%s.y_raw.dtseries.nii'%(hemi_,subj,roi))

	return nb.load(fname)

def load_hrf_tune(subj, glm, roi, param=[6,16], hemi='L', map_='beta'):
	"""
	Output
		y: cifti data with (# total TRs) X (# voxels) shape
	"""
	dir_glm = ut.get_dir_glm(glm)
	glm_ = dir_glm.split('/')[-1]
	dir_work = join(dir_glm,subj,'hrf_tune')

	if np.array(param).dtype == 'int64':
		param_ = str(param).replace(' ','').replace(',',' ')
	else:
		param_ = param
	param_ = param_.replace(']','?')

	hemi_ = hemi.upper()
	fname = glob(join(dir_work,'cifti.%s.%s.%s.%s.%s.%s.*.nii'%(hemi_,glm_,param_,subj,roi,map_)))[0]
	cii = nb.load(fname)

	return cii

def get_df_y(subj, glm, roi, param=[6,16], hemi='L', show_yraw=False, melt=False):
	"""
	Output
		df: DataFrame
			len(df) = # total TRs
			len(df.melt) = nRuns*nTRs(=# TRs per Run)*2(y_adj/y_hat)
	"""
	y_hat = load_hrf_tune(subj=subj,glm=glm,roi=roi,param=param,map_='y_hat')
	y_res = load_hrf_tune(subj=subj,glm=glm,roi=roi,param=param,map_='y_res')

	nTRs = 410
	nRuns = 8

	df = pd.DataFrame(
		{
			'run':np.repeat(np.arange(1,nRuns+.1), nTRs).astype(int),
			'TR':np.tile(np.arange(nTRs),nRuns),
			'y_hat':y_hat.get_fdata().mean(axis=1),
			'y_res':y_res.get_fdata().mean(axis=1),
		}
	)
	df['y_adj'] = df.y_hat + df.y_res
	value_vars = ['y_hat','y_res','y_adj']
	if show_yraw:
		y_raw = load_yraw(subj=subj, roi=roi)
		df['y_raw'] = y_raw.get_fdata().mean(axis=1)
		value_vars = ['y_raw','y_hat','y_res','y_adj']

	if melt:
		df = df.melt(
			id_vars=['run','TR'],
			value_vars=value_vars, var_name='hue', value_name='y'
		)
	
	return df

def get_df_window_y(subj, glm, roi, param, pre=10, post=20, TR=1):
	"""
	Output
		df: DataFrame
			len(df) = nRuns*nTrials*nTRs(=pre+post+1)*2(y_adj/y_hat)
	"""
	## load onset times
	dir_glm = ut.get_dir_glm(glm)
	SPM = join(dir_glm,subj,'SPM.mat')
	df_onset = deal_spm.get_df_onset(SPM)

	## load y
	df_y = get_df_y(
		subj=subj, glm=glm, roi=roi, param=param,
		hemi='L', show_yraw=False, melt=False
	)

	nTRs = len(df_y.TR.unique())
	runs = df_onset.run.unique()

	## shape=(# runs, # trials)
	onset_idxs_by_run = []
	for run in runs:
		onset_idxs_by_run.append(
			np.round(
				np.sort(
					np.concatenate(df_onset[df_onset.run==run].onset.values) / TR
				)
			).astype(int)
		)

	lines = {
		'run':[], 'trial':[], 'TR':[], 'y':[], 'hue':[]
	}
	for rr, onset_idxs in enumerate(onset_idxs_by_run):
		run = rr+1
		for tt, idx in enumerate(onset_idxs):
			trial = tt+1
			start_idx = int(idx - pre)
			end_idx = int(idx + post + 1)

			valid_start = int(max(0, start_idx))
			valid_end = int(min(nTRs, end_idx))

			window_start = int(valid_start - start_idx)
			window_end = int(window_start + (valid_end - valid_start))

			for hue in ['y_adj','y_hat']:
				y = df_y[df_y.run==run][hue]
				window_y = np.full(end_idx - start_idx, np.nan)
				window_y[window_start:window_end] = y[valid_start:valid_end]
				for idx, y in enumerate(window_y):
					TR = idx - pre
					lines['run'].append(run)
					lines['trial'].append(trial)
					lines['TR'].append(TR)
					lines['hue'].append(hue)
					lines['y'].append(y)

	df_window_y = pd.DataFrame(lines)

	return df_window_y
