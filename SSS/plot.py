from os.path import join
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps, cm
import matplotlib.colors as mcolors
import seaborn as sns
import scipy
from nilearn import plotting, image
import nibabel as nb
import SSS.util as ut

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

def cmap_to_cbar(label_list, cmap):
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

def get_WPM(subj, glm):
	dir_surf = ut.get_dir_surf()
	dir_glm = join(ut.get_dir_root(),glm)
	S_id = ut.get_S_id(subj)

	white = join(dir_surf,S_id,'%s.L.white.32k.surf.gii'%S_id)
	pial = join(dir_surf,S_id,'%s.L.pial.32k.surf.gii'%S_id)

	mask = join(dir_glm,subj,'mask.nii')

	return white, pial, mask
