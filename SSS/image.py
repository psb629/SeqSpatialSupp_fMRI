from os.path import join
from os import getcwd, makedirs
from glob import glob

import numpy as np
import pandas as pd

from SSS import util as su
from SSS import deal_spm

from nilearn import image
import nibabel as nb

import surfAnalysisPy as surf
# import Functional_Fusion.atlas_map as am
# import Functional_Fusion.reliability as rel

def get_list_roi():

	return np.array(['S1', 'M1', 'PMd', 'PMv', 'SMA', 'V1', 'SPLa', 'SPLp'])

def get_underlay(path_surfAnalysisPY, hemi='L'):
	hemi_ = hemi.upper()
	path_underlay = join(path_surfAnalysisPY,'standard_mesh/fs_%s/fs_LR.32k.LR.sulc.dscalar.nii'%hemi)

	return path_underlay

def get_border(path_surfAnalysisPY, hemi='L'):
	hemi_ = hemi.upper()
	path_border = join(path_surfAnalysisPY,'standard_mesh/fs_%s/fs_LR.32k.%s.border'%(hemi_,hemi_))

	## Sulci: superior frontal sulcus (SFS), inferior frontal sulcus (IFS), precentral sulcus (PrCS), central sulcus (CS), postcentral sulcus (PoCS), intraparietal sulcus (IPS), parieto-occipital sulcus (POS), lateral occipital sulcus (LOS), lunate sulcus (LnS), superior temporal sulcus (STS), inferior temporal sulcus (ITS), collateral sulcus (CoS), sylvian fissure (SF)
	labels = {}
	labels['PrCS'] = [-20, 110]
	labels['CS'] = [20, 125]
	labels['PoCS'] = [50, 115]
	labels['SFS'] = [-70, 70]
	labels['IFS'] = [-70, 15]
	labels['IPS'] = [100, 90]
	labels['POS'] = [145, 125]
	labels['STS'] = [95, 30]

	return path_border, labels

def load_mask(subj, glm=1, as_nii=True):
	"""
	Return
		mask image: nifti
	"""
	dir_glm = su.get_dir_glm(glm)
	mask = join(dir_glm,subj,'mask.nii')
	img = nb.load(mask)

	if as_nii:
		return img
	else:
		return img.get_fdata(), img.affine, img.header

def load_summed_roi(subj, list_roi):
	"""
	Return
		ROI_img: 3D nifti data
			ROI image with sequentially assigned natural number labels for each ROI.
	"""
	dir_roi = su.get_dir_roi()
	S_id = su.get_S_id(subj)
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

def load_yraw(subj, run=None, roi=None, hemi='L'):
	"""
	Return
		y: cifti data with (# total TRs) X (# voxels) shape
	"""
	if roi is None:
		assert run is not None, 'The input, RUN, is necessary.'
		dir_raw = join(su.get_dir_root(),'imaging_data')
		## load the y_raw
		fname = join(dir_raw,subj,'%s_run_%02d.nii'%(subj,run))
	else:
		assert run is None, 'The input, RUN, is not necessary.'
		dir_roi = su.get_dir_roi()
		hemi_ = hemi.upper()
		fname = join(dir_roi,subj,'cifti.%s.%s.%s.y_raw.dtseries.nii'%(hemi_,subj,roi))

	return nb.load(fname)

def trim_ydata(ydata, glm, subj, as_nii=True):
	"""
	Return
		ydata: nifti 
			time series, nVolumns x # nTRs
		glm: int
	"""
	## RUNs
	list_run = su.get_list_run()

	## Load the design matrix X of the 1st Run
	X = deal_spm.get_SPM_X(SPM=deal_spm.fname_SPM(subj=subj, glm=glm))
	T,K = X.shape

	## Load the mask image
	_, affine, header = load_mask(subj=subj, glm=glm, as_nii=False)

	## Get the values from ydata
	if not (isinstance(ydata, np.memmap))|(isinstance(ydata, np.ndarray)):
		ydata = ydata.get_fdata()
	_,P = ydata.shape

	## Trim the ydata
	ydata = ydata.reshape(len(list_run),-1,P)[:,-T:,:].reshape(-1,P)
	
	if as_nii:
		return nb.Nifti1Image(ydata, affine=affine, header=header)
	else:
		return ydata

def masking_data(data, mask):
	"""
	Return
		data: nifti 
			Volumn X K
		mask: nifti
			Volumn
	"""
	affine = mask.affine
	header = mask.header
	assert (data.affine == affine).all(), "The data and mask are not matched"
	assert data.shape[:3] == mask.shape, "The data and mask are not matched"

	data = data.get_fdata()
	mask = mask.get_fdata()
	#mask[mask==0] = np.nan

	res = np.ones(data.shape) * np.nan
	if len(data.shape)>3:
		for t in np.arange(data.shape[-1]):
			res[...,t] = data[...,t] * mask
	else:
		res = data * mask

	return nb.Nifti1Image(res, affine=affine, header=header)

def load_hrf_tune(subj, glm, roi, param=[6,16], hemi='L', map='beta'):
	"""
	Return
		y: cifti 
			(# total TRs) X (# voxels) for map='y_hat' or 'y_res'
		beta: cifti
			(# runs * # interests) X (# voxels) for map='beta'
	"""
	dir_glm = su.get_dir_glm(glm)
	glm_ = dir_glm.split('/')[-1]
	dir_work = join(dir_glm,subj,'hrf_tune')

	param_ = su.convert_param_to_hrf(params=param, type='str')

	hemi_ = hemi.upper()
	fname = join(dir_work,'cifti.%s.%s.%s.%s.%s.%s.*.nii'%(hemi_,glm_,param_.replace('[','?'),subj,roi,map)) 
	fnames = glob(fname)
	if len(fnames)>0:
		cii = nb.load(fnames[0])

	return cii

def get_df_y(subj, glm, roi, param=[6,16], hemi='L', show_yraw=False, melt=False):
	"""
	Return
		df: DataFrame
			len(df) = # total TRs
			len(df.melt) = nRuns*nTRs(=# TRs per Run)*2(y_adj/y_hat)
	"""
	y_hat = load_hrf_tune(subj=subj,glm=glm,roi=roi,param=param,map='y_hat')
	y_res = load_hrf_tune(subj=subj,glm=glm,roi=roi,param=param,map='y_res')

	nTRnRUN, _ = y_hat.shape
	nRuns = 8
	nTRs = int(nTRnRUN/nRuns)

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
		df['y_raw'] = y_raw.get_fdata().mean(axis=1)[-nTRnRUN:]
		value_vars = ['y_raw','y_hat','y_res','y_adj']

	if melt:
		df = df.melt(
			id_vars=['run','TR'],
			value_vars=value_vars, var_name='hue', value_name='y'
		)
	
	return df

def get_df_window_y(subj, glm, roi, param, pre=10, post=20, gap=0, TR=1):
	"""
	Return
		df: DataFrame
			len(df) = nRuns*nTrials*nTRs(=pre+post+1)*2(y_adj/y_hat)
	"""
	## load onset times
	dir_glm = su.get_dir_glm(glm)
	SPM = join(dir_glm,subj,'SPM.mat')

	if not gap:
		onsets_by_run = deal_spm.get_concat_onset(SPM)
	else:
		onsets = np.array(deal_spm.get_concat_onset(SPM))
		idxs = np.diff(onsets,axis=1) > gap
		tmp = np.ones((len(onsets),1), dtype=bool)
		idxs = np.concatenate([idxs, tmp], axis=1)

		onsets_by_run = onsets[idxs].reshape((len(onsets),-1))
	
	## load y
	df_y = get_df_y(
		subj=subj, glm=glm, roi=roi, param=param,
		hemi='L', show_yraw=False, melt=False
	)

	nTRs = len(df_y.TR.unique())
	# runs = df_onset.run.unique()

	## shape=(# runs, # trials)
	onset_idxs_by_run = []
	for onsets in onsets_by_run:
		onset_idxs_by_run.append(
			np.round(onsets/TR).astype(int)
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

def get_WPM(subj, glm, hem='L'):
	"""
	Return
		cerebral images: gifti, gifti, and nifti
	"""
	dir_surf = su.get_dir_surf()
	dir_glm = su.get_dir_glm(glm)
	S_id = su.get_S_id(subj)

	white = join(dir_surf,S_id,'%s.%s.white.32k.surf.gii'%(S_id,hem))
	pial = join(dir_surf,S_id,'%s.%s.pial.32k.surf.gii'%(S_id,hem))

	mask = join(dir_glm,subj,'mask.nii')

	return white, pial, mask

def calc_sigma(subj, glm, roi, hemi, param):
	"""
	Return
		Sigma: 2-D numpy array
			Covariance matrix of Y_res for the given ROI
	"""

	y_res = load_hrf_tune(subj=subj, glm=glm, roi=roi, param=param, hemi=hemi, map='y_res')
	T, P = y_res.shape

	sigma = (y_res.get_fdata().T @ y_res.get_fdata()) / T
	
	return sigma

def get_prewhitened_beta(subj, glm, region='2D', param=[5,15], hemi='L'):
	"""
	Return
		beta_whiten : 1-D or 2-D numpy array
			cifti data with (# runs * # interests) X (# voxels) for map='beta'
	"""
	if region=='2D':
		dir_surf = su.get_dir_surf()
		betas = nb.load(join(dir_surf,'glm_%d/%s.%s.glm_%d.beta.func.gii'%(glm,subj,hemi,glm))).darrays
		res = nb.load(join(dir_surf,'glm_%d/%s.%s.glm_%d.ResMS.func.gii'%(glm,subj,hemi,glm))).darrays[0].data

	else:
		pp = su.convert_param_to_hrf(params=param, type='str')
		betas = load_hrf_tune(subj=subj, glm=glm, roi=region, param=param, hemi=hemi, map='beta').get_fdata()
		## Sometimes 'res' contains NaN values—for example, in the M1 region of R11—but the reason is unknown.
		# res = nb.load(join(
		# 	su.get_dir_roi(),'glm_%d'%glm,'cifti.%s.%s.%s.ResMS.dscalar.nii'%(hemi,subj,region)
		# )).get_fdata().reshape(-1)
		sigma = calc_sigma(subj=subj, glm=glm, roi=region, hemi=hemi, param=param)
		res = np.diagonal(sigma)
		
	beta_whiten = []
	for beta in betas:
		if isinstance(beta, nb.gifti.gifti.GiftiDataArray):
			beta = beta.data
		beta_whiten.append(
			beta/(np.sqrt(res)+1.e-14)
		)

	return np.nan_to_num(beta_whiten, nan=0.0)

def get_optimal_hrf(subj, roi, r2_score=None):
	"""
	Select the HRF parameter with the highest R2 score in ROI.

	Return
		param : string
			HRF parameters
	"""
	df_r2 = pd.read_csv(r2_score, sep='\t', header=0)

	df_tmp = df_r2.groupby(['subj','roi','param'], as_index=False).mean(['r2'])
	df_param = df_tmp[df_tmp.r2==df_tmp.groupby(['subj','roi'])['r2'].transform('max')]
	df_param.sort_values(by=['subj','roi'], ascending=[True,True])
	# df_param['subj'] = df_param.subj.astype(str).str.zfill(2)

	if isinstance(subj,int):
		sidx = subj
	else:
		if len(subj)==3:
			sidx = int(subj[1:])
		elif len(subj)==2:
			sidx = int(subj)

	param = df_param[(df_param.subj==sidx)&(df_param.roi==roi)].param.values[0]
	
	return param

def load_contrast_order(subj, glm, map='con'):
	dir_work = join(su.get_dir_surf(),'glm_%d'%glm)

	order = np.genfromtxt(
		join(dir_work,'%s.%s_orders.csv'%(subj,map)),
		delimiter='\t', dtype=str
	)

	return order

def save_surf2cifti(data, label_axis, dir_output, prefix='p'):
	"""
	Save the data as a CifTi file.
	"""
	makedirs(dir_output, exist_ok=True)
	bm_axis = nb.cifti2.BrainModelAxis.from_surface(
		vertices=np.arange(32492), nvertex=32492, name='CortexLeft'
	)
	scalar_axis = nb.cifti2.ScalarAxis(label_axis)
	header = nb.Cifti2Header.from_axes((scalar_axis, bm_axis))

	cii = nb.Cifti2Image(dataobj=data, header=header)
	nb.save(cii, join(dir_output, '%s.dscalar.nii'%prefix))
