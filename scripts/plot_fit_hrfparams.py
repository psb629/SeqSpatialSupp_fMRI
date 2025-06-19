#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import platform
from os.path import join, exists, abspath, dirname
from os import getcwd, makedirs
from glob import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps, cm, colors
import scipy
import h5py
import seaborn as sns

from tqdm import tqdm


# In[2]:


import nibabel as nb
from nilearn import plotting, image
from nipype.interfaces import fsl


# In[3]:


dir_current = getcwd().replace('\\','/')

tmp = dir_current.split('/')
idx = [ii for ii, s in enumerate(tmp) if s=='github'][0]

dir_git = '/'.join(tmp[:idx+1])
dir_git


# In[4]:


dname = join(dir_git,'nitools')
sys.path.append(dname)
import nitools as nt


# In[5]:


dname = join(dir_git,'SUITPy')
sys.path.append(dname)
import SUITPy as suit


# In[6]:


dname = join(dir_git)
sys.path.append(dname)
import surfAnalysisPy as surf


# In[7]:


dname = join(dir_git,'SeqSpatialSupp_fMRI')
sys.path.append(dname)
from SSS import deal_spm
from SSS import util as su
from SSS import plot as splt
from SSS import image as simage


# ---

# In[38]:

glm = 1
dir_glm = su.get_dir_glm(glm)


# In[39]:


dir_result = su.get_dir_result()
dir_work = join(dir_result,'mean_y_across_run')
makedirs(dir_work, exist_ok=True)


# In[40]:


nrows, ncols = 6, 1

for ss in ['S','R']:
    for nn in su.get_list_sn():
        subj = ss+nn
        SPM = join(dir_glm,subj,'SPM.mat')
        df_onset = deal_spm.get_df_onset(SPM)
        onsets_by_run = []
        for rr in range(8):
            run = rr+1
            onsets_by_run.append(
                np.sort(
                    np.concatenate(df_onset[df_onset.run==run].onset.values)
                ).astype(int)
            )
        for roi in simage.get_list_roi():
            fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(30*ncols,4*nrows))
            plt.close()
            for ii, param in enumerate([[4,14],[5,15],[6,16],[7,17],[8,18],[9,19]]):
                ax = axs[ii]
                df_y = simage.get_df_y(subj=subj,glm=glm,roi=roi,param=param,hemi='L',show_yraw=False,melt=True)
                g = sns.lineplot(
                    data=df_y[df_y.hue!='y_res'],
                    x='TR', y='y', hue='hue',
                    ax=ax
                )
                handles, labels = g.get_legend_handles_labels()
                g.legend(handles, [r'$%s_{%s}$'%(s.split('_')[0],s.split('_')[1]) for s in labels], loc='upper left', fontsize=12)
                g.grid(axis='x', linestyle='--', color='gray')
                g.set_ylabel(r'mean $y$ across run', fontsize=16)
                g.set_xlabel('TR', fontsize=16)
                xticks = g.get_xticks()[1:-1]
                g.set_xticks(xticks)
                g.set_xticklabels(['%d'%x for x in xticks], fontsize=14)
                yticks = g.get_yticks()[1:-1]
                g.set_yticks(yticks)
                g.set_yticklabels(['%.2f'%y for y in yticks], fontsize=14)
                g.set_title('%s (%s, %s)'%(deal_spm.convert_param_to_hrf(params=param,type='str'), subj, roi))

                for onset in onsets_by_run[ii]:
                    g.axvline(x=onset, color='red', linestyle='--')

            fig.tight_layout()
            # plt.show()
            fig.savefig(
                join(dir_work,'y_mean.%s.%s.png'%(subj,roi)),
                dpi=300, facecolor=[1,1,1,1],
                bbox_inches='tight'
            )


# #### 4. Time series time lock to stimulus onset

# In[41]:


glm = 1


# In[42]:


dir_result = su.get_dir_result()
dir_work = join(dir_result,'y_window')
makedirs(dir_work, exist_ok=True)


# In[41]:


nrows, ncols = 4, 2

for ss in ['S','R']:
    for nn in su.get_list_sn():
        subj = ss+nn
        for roi in simage.get_list_roi():
            for param in [[4,14],[5,15],[6,16],[7,17],[8,18],[9,19]]:
                df_window_y = simage.get_df_window_y(subj=subj,glm=glm,roi=roi,param=param)

                fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6*ncols,4*nrows))
                axs = np.concatenate(axs)
                plt.close()
                for rr in range(8):
                    ax = axs[rr]
                    run = rr+1
                    g = sns.lineplot(
                        data=df_window_y[df_window_y.run==run],
                        x='TR', y='y', hue='hue',
                        ax=ax
                    )
                    handles, labels = g.get_legend_handles_labels()
                    g.legend(handles, [r'$%s_{%s}$'%(s.split('_')[0],s.split('_')[1]) for s in labels], loc='upper left', fontsize=12)
                    g.grid(axis='x', linestyle='--', color='gray')
                    g.set_ylabel(r'mean $y$ across window', fontsize=16)
                    g.set_xlabel('TR', fontsize=16)
                    xticks = g.get_xticks()[1:-1]
                    g.set_xticks(xticks)
                    g.set_xticklabels(['%d'%x for x in xticks], fontsize=14)
                    yticks = g.get_yticks()[1:-1]
                    g.set_yticks(yticks)
                    g.set_yticklabels(['%.2f'%y for y in yticks], fontsize=14)
                    g.set_title('run%02d (%s, %s, %s)'%(run, subj, roi, str(param)))
                    g.axvline(x=0, color='red', linestyle='-')

                fig.tight_layout()
                # plt.show()
                fig.savefig(
                    join(dir_work,'y_window.%s.%s.%s.png'%(subj,roi,deal_spm.convert_param_to_hrf(params=param,type='str'))),
                    dpi=300, facecolor=[1,1,1,1],
                    bbox_inches='tight'
                )


# ---

# ## Load $\beta$

