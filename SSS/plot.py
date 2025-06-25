from os.path import join

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import colormaps, cm
import matplotlib.colors as mcolors

from SSS import util as ut
from SSS import deal_spm

import nibabel as nb

def gifti_to_cmap(label_img):
	if isinstance(label_img, str):
		label_img = nb.load(label_img)
	if not isinstance(label_img, nb.gifti.gifti.GiftiImage):
		raise TypeError("Expected 'gifti' file")

	labels = label_img.labeltable.labels
	rgba = np.zeros((len(labels),4))
	for i, label in enumerate(labels[1:]):
		rgba[i+1,:] = labels[i].rgba
		
	cmap = mcolors.ListedColormap(rgba, N=len(rgba))

	return cmap

def cmap_for_cbar(label_list, cmap):
	N = len(label_list)

	cmap_ = colormaps.get_cmap(cmap).resampled(N)
	colors = [cmap_(i) for i in range(N)]
	colors.insert(0,(0,0,0,1))
	cmap_ = mcolors.LinearSegmentedColormap.from_list('modified_%s'%cmap, colors, N+1)

	norm = mcolors.BoundaryNorm(boundaries=np.arange(-0.5,N+1,1), ncolors=N+1)

	fig, ax = plt.subplots(figsize=(8,1))
	fig.subplots_adjust(bottom=0.5)

	sm = plt.cm.ScalarMappable(cmap=cmap_, norm=norm)
	sm.set_array([])

	cbar = fig.colorbar(sm, cax=ax, orientation='horizontal', ticks=np.arange(N+1))

	labels = np.concatenate((['0'], label_list))
	cbar.ax.set_xticklabels(labels)

	plt.show()

	return cmap_

def plot_SPM_X(subj, glm, run=1):
	dir_glm = join(ut.get_dir_root(),glm)
	SPM = join(dir_glm,subj,'SPM.mat')

	df_onset = deal_spm.get_df_onset(SPM)

	X = get_SPM_X(SPM, run=run)

	df = pd.DataFrame()
	for ii, reg in enumerate(df_onset.reg.unique()):
		plt.plot(X[:,ii], label=reg)
		df['onset'] = df_onset[(df_onset.run==run)&(df_onset.reg==reg)]
		df['stim'] = 0.1

def plot_BF(xBF):
	xBF_ = deal_spm.get_xBF_params(xBF)

	# x = np.arange(xBF_['T0'],xBF_['length']+xBF_['T0'],xBF_['dt'])
	t0 = xBF_['T0']
	dt = xBF_['dt']
	x = np.arange(t0,xBF_['length']+t0,dt)
	y = xBF_['bf']

	fig, ax = plt.subplots()

	ax.plot(x, y)
	ax.grid(axis='y')
	ax.set_xlim(-1,34)

	y = 0.007
	dy = 0.001

	# x = xBF_['T0']
	x = xBF_['params'][5]
	ax.axvline(x=x, color='gray', linestyle='--', linewidth=2)
	ax.text(x=x, y=y+dy, s='onset', ha='center', va='center')

	x = xBF_['length']+t0
	ax.axvline(x=x, color='gray', linestyle='--', linewidth=2)
	ax.text(x=x, y=y+dy, s='offset', ha='center', va='center')

	x = xBF_['params'][0]
	ax.axvline(x=x, color='blue', linestyle='--', linewidth=2)
	ax.text(x=x, y=y-dy, s='response', ha='center', va='center')

	x = xBF_['params'][1]
	ax.axvline(x=x, color='red', linestyle='--', linewidth=2)
	ax.text(x=x, y=y+dy, s='undershoot', ha='center', va='center')

	ax.set_title(xBF_['params'])
	plt.show()

