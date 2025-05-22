import platform
from os.path import join
from os import getcwd
import numpy as np
import h5py
import re
import SSS.util as ut
import pandas as pd

def get_trials_per_sess(SPM):
	if isinstance(SPM, str):
		SPM = ut.load_spm(SPM)

	tmp = SPM['xsDes/Trials_per_session'][:]
	#if tmp.dtype == np.uint16:
	tmp = np.reshape(tmp,-1)
	tmp = ''.join(map(chr, tmp))
	## len(tps): # of runs, tps[run]: # of regressors
	tps = np.array(re.findall(r'\d+', tmp)).astype(int)

	return tps

def get_TR(SPM):
	if isinstance(SPM, str):
		SPM = ut.load_spm(SPM)

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
	if isinstance(SPM, str):
		SPM = ut.load_spm(SPM)

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

def get_onset(SPM):
	if isinstance(SPM, str):
		SPM = ut.load_spm(SPM)

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
