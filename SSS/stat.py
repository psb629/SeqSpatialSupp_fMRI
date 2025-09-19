import platform

from os.path import join
from os import getcwd

import numpy as np
import pandas as pd
import re
import scipy

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

def convert_alpha_to_tval(alpha=0.05, df=None, alternative='two-sided'):
	if alternative == 'one-sided':
		alpha = alpha * 2
	thresh = scipy.stats.t.ppf(1-0.5*alpha, df=df)
	
	return thresh

def compute_t_stat(betas, conditions):
	"""
	# Input:
		betas
			2D array (n_trials, n_channels)
		conditions
			list of trial indices for each condition (0-based indexing)
	
	# Return: 
		t_stat
			2D array (n_conditions, n_channels)
	"""
	conditions = np.array(conditions)
	sz = betas.shape
	n_conditions = len(np.unique(conditions))
	t_stat = np.zeros((n_conditions, sz[-1]))

	for ii, cond in enumerate(np.unique(conditions)):
		idx = [True if c==cond else False for c in conditions]
		# 조건 trial들의 beta 값 추출
		betas_cond = betas[idx,:]
		# 평균 (nan 무시)
		mean_beta = np.nanmean(betas_cond, axis=0)
		# 표준편차 (nan 무시)
		std_beta = np.nanstd(betas_cond, axis=0, ddof=1)  # ddof=1 → 표본 표준편차
		# 표준오차
		se_beta = std_beta / np.sqrt(len(idx)) + 1.e-14
		# t-statistic
		t_stat[ii] = mean_beta / se_beta

	return t_stat

