import platform
from os.path import join
from os import getcwd
import numpy as np
import h5py
import re
import SSS.util as ut
import pandas as pd


def get_trials_per_sess(SPM):
	SPM = ut.load_SPM(SPM)

	tmp = SPM['xsDes/Trials_per_session'][:]
	#if tmp.dtype == np.uint16:
	tmp = np.reshape(tmp,-1)
	tmp = ''.join(map(chr, tmp))
	## len(tps): # of runs, tps[run]: # of regressors
	tps = np.array(re.findall(r'\d+', tmp)).astype(int)

	return tps

def get_TR(SPM):
	SPM = ut.load_SPM(SPM)

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
	SPM = ut.load_SPM(SPM)

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
	SPM = ut.load_SPM(SPM)

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

def get_df_vec(SPM):
	SPM = ut.load_SPM(SPM)
	
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

	return df.filter(items=['part_vec','cond_vec'])

def get_SPM_X(SPM, run=1):
	SPM = ut.load_SPM(SPM)

	nRuns = len(SPM['xX/K/row'])

	## row
	rr = run-1
	ref = SPM['xX/K/row'][rr,0]
	idx_r = SPM[ref][:].reshape(-1).astype(int)

	## column
	iC = SPM['xX/iC'][:].reshape(-1).astype(int)
	

	X = SPM['xX/X'][:].T

	return X[idx_r,:][:,idx_c]
