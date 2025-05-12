# SeqSpatialSupp_fMRI

## 1. Preprocessing

### 1-1. Anatomical Image

```
sss_imana('PREP:ANAT-all','sn',<subject number>)
```

#### 1-1-i. Pre-freesurfer

#### 1-1-ii. freesurfer

i) recon-all
```
scripts/reconall.zsh -s <subject id>
```

ii) reslice
```
scripts/reslice.zsh -s <subject id> -a <surface name>
```

### 1-2. Functional Image

i) Automatic run (pre):
```
sss_imana('PREP:FUNC-pre','sn',<subject number>)
```

ii) Manual run:
- Open `fsleyes`
- Add anatomical image and b*mean.nii* (bias corrected mean) image to overlay
- click on the bias corrected mean image in the ‘Overlay list' in the bottom left of the fsleyes window. list to highlight it.
- Open `tools` -> `Nudge`
- Manually adjust *bmean.nii* image to the anatomical by changing the 6 parameters (translation xyz and rotation xyz). Do not change the **scales**!
- When done, click apply and close the tool tab. Then to save the changes, click on the save icon next to the mean image name in the ‘Overlay list’ and save the new image by adding ‘r’ in the beginning of the name: r*bmean.nii*. If you don’t set the format to be .nii, fsleyes automatically saves it as a .nii.gz so either set it or gunzip afterwards to make it compatible with SPM.

iii) Automatic run (post)
```
sss_imana('PREP:FUNC-post','sn',<subject number>)
```

#### 1-2-i. field map (optional)

#### 1-2-ii. EPI

---

## 2. GLM 

```
sss_GLM('GLM:all','sn',<subject number>,'glm',<GLM number>)
```

### GLM number

Index (i,j)
- Sequence
	- i=0: 32451
	- i=1: 35124
	- i=2: 13254
	- i=3: 14523
- Cue
	- j=0: Letter
	- j=1: Spatial

i) GLM = 2: Repetition

|  | trial $t-1$ | trial $t$ |
|---------|---------|---------|
| Both-Rep| $(i,j)$ | $(i,j)$ |
| Cue-Rep | $(i,\neg j)$ | $(i,j)$ |
| Seq-Rep | $(\neg i,j)$ | $(i,j)$ |
| NRep    | $(\neg i,\neg j)$ | $(i,j)$ |

|     |(0,0)|(0,1)|(1,0)|(1,1)|(2,0)|(2,1)|(3,0)|(3,1)| 
|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|(0,0)|  B  | S00 | C01 | N02 | C03 | N04 | C05 | N06 |
|(0,1)|     |  B  | N07 | C08 | N09 | C10 | N11 | C12 |
|(1,0)|     |     |  B  | S13 | C14 | N15 | C16 | N17 |
|(1,1)|     |     |     |  B  | N18 | C19 | N20 | C21 |
|(2,0)|     |     |     |     |  B  | S22 | C23 | N24 |
|(2,1)|     |     |     |     |     |  B  | N25 | C26 |
|(3,0)|     |     |     |     |     |     |  B  | S27 |
|(3,1)|     |     |     |     |     |     |     |  B  |

ii) GLM = 3: Trial State (i,j)
- 1: (0,0)
- 2: (0,1)
- 3: (1,0)
- 4: (1,1)
- 5: (2,0)
- 6: (2,1)
- 7: (3,0)
- 8: (3,1)
