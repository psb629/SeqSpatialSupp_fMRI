from os.path import join
from glob import glob

import numpy as np
import pandas as pd

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

def load_map_order(subj, glm, map):
	dir_surf = join(get_dir_glmsingle(glm),'surfaceWB')
	fname = join(dir_surf,subj,'%s.%s_orders.csv'%(subj,map))
	order = np.loadtxt(fname, delimiter='\t',dtype=str)

	return order
