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

def get_dir_glmsingle(glm):

	return join(ut.get_dir_root(),'GLMsingle','glm_%d'%glm)

def load_reginfo(subj, glm):
	dir_glm = get_dir_glmsingle(glm)
	reginfo = pd.read_csv(join(dir_glm,subj,'reginfo.tsv'), sep='\t')

	fnames = []
	for idx in reginfo.index.values:
		fname = 'beta_%04d.nii'%(idx+1)
		fnames.append(fname)
	reginfo['beta']=fnames

	return reginfo
