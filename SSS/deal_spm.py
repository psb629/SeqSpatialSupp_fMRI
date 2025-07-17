from os.path import join
from glob import glob

import numpy as np
import pandas as pd

import re

import h5py

from SSS import util as ut

def load_SPM(SPM):
	"""
	Output
		SPM: h5py data
			Data loaded from SPM.mat, the GLM result file processed using SPM.
	"""
	if isinstance(SPM, str):
		file = h5py.File(SPM)
		SPM = file['SPM']

	return SPM

def get_trials_per_sess(SPM):
	"""
	Output
		tps: N-D numpy array
			Output the number of trials for each session (run).
	"""
	SPM = load_SPM(SPM)

	tmp = SPM['xsDes/Trials_per_session'][:]
	#if tmp.dtype == np.uint16:
	tmp = np.reshape(tmp,-1)
	tmp = ''.join(map(chr, tmp))
	## len(tps): # of runs, tps[run]: # of regressors
	tps = np.array(re.findall(r'\d+', tmp)).astype(int)

	return tps

def get_TR(SPM):
	"""
	Output
		tr: float
			The number of TRs
		unit: string
			The unit of TR
	"""
	SPM = load_SPM(SPM)

	tmp = SPM['xsDes/Interscan_interval'][:]
	tmp = np.reshape(tmp,-1)
	tmp = ''.join(map(chr, tmp))
	## TR
	tr = np.array(re.findall(r'\d+\.\d+',tmp)).astype(float)
	## Unit
	match = re.search(r'\{(.*?)\}', tmp)

	if match:
	    unit = match.group(1)

	return tr, unit

def get_column_head(SPM):
	"""
	Output
		column_head: 1-D numpy array
			Labels of the GLM betas of interest.
	"""
	SPM = load_SPM(SPM)

	U = SPM['Sess/U']
	nruns = len(U)
	
	column_head = []
	for run in np.arange(nruns):
		ref = U[run,0]
		name = SPM[ref]['name']
		ntrial = len(name)
		for trial in np.arange(ntrial):
			ref = SPM[name[trial,0]][0,0]
			reg = SPM[ref][:]
			reg = np.reshape(reg,-1)
			reg = ''.join(map(chr, reg))
			column_head.append(reg)

	return np.array(column_head)

def get_df_onset(SPM):
	"""
	Output
		df_onset: pandas dataframe
			Data storing the onset times of each trial.
	"""
	SPM = load_SPM(SPM)

	column_head = get_column_head(SPM)

	U = SPM['Sess/U']
	nruns = len(U)
	
	df = {'run':[],'reg':[],'onset':[]}
	for run in np.arange(nruns):
		ref = U[run,0]
		ons = SPM[ref]['ons']
		ntrial = len(ons)
		for trial in np.arange(ntrial):
			ref = ons[trial,0]
			onset = SPM[ref][:].reshape(-1)
			df['run'].append(run+1)
			idx = trial + run*ntrial
			df['reg'].append(column_head[idx])
			df['onset'].append(onset)

	return pd.DataFrame(df)

def get_concat_onset(SPM):
	df_onset = get_df_onset(SPM)

	runs = df_onset.run.unique()
	onsets_by_run = []
	for run in runs:
		onsets = np.sort(
			np.concatenate(
				df_onset[df_onset.run==run].onset.values
			)
		)
		onsets_by_run.append(onsets)
	
	return onsets_by_run

def get_df_vec(SPM):
	"""
	A regressor in 'label (SPM/Sess/U/name)' that contains the string 'Non' is assigned condition = 0, and the rest are numbered in ascending order starting from 1.

	Input
		SPM: str / h5py
			It is either the filename of the SPM.mat file or the data itself loaded using h5py.
	
	Output
		df: Pandas DataFrame
			A dataframe storing the condition vector and partition vector.
	"""
	SPM = load_SPM(SPM)
	
	df = get_df_onset(SPM).filter(items=['run','reg'])
	nRuns = len(df.run.unique())

	df['cond_vec'] = df.apply(lambda row: 0 if 'Non' in row['reg'] else None, axis=1)

	for r in np.arange(nRuns)+1:
		df_tmp = df[df.run==r]
		cnt = 1
		for row in df_tmp[np.isnan(df_tmp.cond_vec)].index:
			df.loc[row,'cond_vec'] = cnt
			cnt += 1
	df.cond_vec = df.cond_vec.astype(int)
	
	df.rename(columns={'run':'part_vec'}, inplace=True)
	df = df.filter(items=['part_vec','cond_vec'])

	column_head = get_column_head(SPM)
	df['column_head'] = column_head

	return df

def get_SPM_X(SPM, run=1):
	"""
	Output
		X: 2-D numpy array
			Design matrix for the regressors of interest for each run.
	"""
	SPM = load_SPM(SPM)

	nRuns = len(SPM['xX/K/row'])

	## row
	rr = run-1
	ref = SPM['xX/K/row'][rr,0]
	idx_r = SPM[ref][:].reshape(-1).astype(int)-1

	## column
	iC = SPM['xX/iC'][:].reshape(-1).astype(int)-1
	idx_c = iC.reshape(nRuns,-1)[rr,:]

	X = SPM['xX/X'][:].T

	return X[idx_r,:][:,idx_c]

def get_df_X(SPM, mean=False):
	"""
	Output
		df: pandas dataframe
			Design matrix for the regressors of interest for all runs.
	"""
	nrun = 8
	head = get_column_head(SPM).reshape(nrun,-1)[0]

	df_list = []
	for rr in range(nrun):
		run = rr+1
		X = get_SPM_X(SPM, run=run)
		df = pd.DataFrame(X)
		nrows, ncols = df.shape
		df['Run']=run
		df['TR']=np.arange(nrows)
		df_list.append(df)
	df = pd.concat(df_list, ignore_index=True)
	df = df.filter(['Run','TR',*np.arange(ncols)])

	columns = {}
	for ii, reg in enumerate(head):
		columns[ii]=reg
		df.rename(columns=columns, inplace=True)

	if mean:
		df['mean'] = df.drop(columns=['Run','TR']).mean(axis=1)
		df = df.filter(['Run','TR','mean'])

	return df

def get_xBF_params(xBF):
	"""
	Output
		xBF: dictionary
			Data containing HRF information.
	"""
	if not isinstance(xBF, dict):
		if isinstance(xBF, str):
			fname = glob(xBF.replace('[','?'))[0]
			file = h5py.File(fname)
			xBF = file['xBF']

    	# xBF = SPM['xBF']
		xBF_ = {}
		for key, feature in xBF.items():
			tmp = feature[:].copy()
			if np.ndim(tmp) > 1:
				tmp = tmp.flatten()
			if len(tmp) == 1:
				tmp = tmp[0]
			if tmp.dtype == np.uint16:
				tmp = ''.join(map(chr, tmp))
			xBF_[key] = tmp
	else:
		xBF_ = xBF

	return xBF_

def load_reginfo(subj, glm):
	dir_glm = ut.get_dir_glm(glm)
	reginfo = pd.read_csv(
		join(dir_glm,subj,'reginfo.tsv'),
		sep='\t', header=0
	)

	return reginfo
