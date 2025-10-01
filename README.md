# SeqSpatialSupp_fMRI

## 1. Preprocessing

### 1-0. Behavioural Data
Create a 'behav_info.tsv' file containing the necessary information for the GLM in each subject's raw directory.
```
scripts/extract_behav_info.py --subject <subj_id>
```

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

Index (s,c)
- Sequence
	- s=0: 32451
	- s=1: 35124
	- s=2: 13254
	- s=3: 14523
- Cue
	- c=0: Letter
	- c=1: Spatial

#### i) GLM = 1: Trial State (s,c)
i-1) onset: trial onset + prep (1s)

i-2) duration: 2s

- 1: (0,0)
- 2: (0,1)
- 3: (1,0)
- 4: (1,1)
- 5: (2,0)
- 6: (2,1)
- 7: (3,0)
- 8: (3,1)

#### ii) GLM = 2: Repetition
ii-1) onset: trial onset + prep (1s)

ii-2) duration: 2s

 | | symbol | trial$_{t-1}$ | trial $_{t}$ |
 |---------|---------|---------|---------|
 | Both-Rep$_{c}$ | $B_{c}$| $(s,c)$ | $(s,c)$ |
 | Cue-Rep$_{c}$ | $C_{c}$| $(\neg s,c)$ | $(s,c)$ |
 | Seq-Rep$_{c}$ | $S_{c}$| $(s,\neg c)$ | $(s,c)$ |
 | Non-Rep$_{c}$ | $N_{c}$| $(\neg s,\neg c)$ | $(s,c)$ |

within Cue RS effect (=$B_{c} - S_{c}$): seq가 반복되는 맥락에서의 cue 반복 효과의 단순효과
- seq 반복 여부를 고정한 상태에서, cue가 반복되었을 때 추가로 생기는 RS(반복 억제) 크기를 추정한다.
- 해석은 “이미 같은 seq를 수행 중일 때 cue가 같으면 얼마나 더 억제되는가?”에 해당한다.

across Cue RS effect (=$C_{c} - N_{c}$): seq가 반복되지 않는 맥락에서의 cue 반복 효과의 단순효과
- seq 비반복을 고정한 상태에서, cue 반복의 RS 크기를 추정한다.
- 해석은 “seq가 바뀔 때라도 cue가 같으면 얼마나 억제되는가?”에 해당한다.

두 대비를 비교할 때 드러나는 것:
- 상호작용 검정: $(B-S)\neq(C-N)$ 이면 cue 반복 효과가 seq 반복 맥락에 의해 변조됨을 뜻한다. 즉, cue×seq 상호작용의 존재를 시사한다.
- 차이의 차이(difference-in-differences): 상호작용 대비는 $(B-S)-(C-N)$로 정리되며, 이 값이 0에서 유의하게 벗어나면 두 요인의 비가법성이 통계적으로 확인된다.

시각적 cue 종류(letter vs. spatial)를 고려한 추가 팁:
- 모달리티 요인 추가 (letter/spatial): cue가 두 종류(문자, 공간)이므로 3요인 설계(모달리티 × cue반복 × seq반복)로 각 모달리티에서의 단순효과와 상호작용을 분리 추정하는 것이 좋다.
- 모달리티별 RS 차이: $(B-S)\neq(C-N)$을 모달리티별로 계산 후, 모달리티 간 차이를 비교하면 “어떤 cue 유형에서 맥락 변조가 더 강한지”를 확인할 수 있다.
- 균형과 전이효과: 각 조합의 시행 수 균형, 블록·러닝에 따른 전이(exposure)와 피로 누적을 공변량으로 통제하면 RS 추정이 안정적이다.

RS effect 의 양(+)의 값과 음(-)의 값: RS effect 는 곧 어떤 자극에 대한 절대적인 활동이 줄어들었다는 의미. 즉, $|A_{t}|>|A_{t+1}|$
- 애초에 어떤 자극에 대한 활동이 + 였으면 ($A_{t}>0$), RS effect 는 - 로 나타남. cf) $A_{t+1}-A_{t}<0$
- 그러나 어떤 자극에 대한 활동이 - 였으면 ($A_{t}<0$), RS effect 는 + 로 나타남. cf) $A_{t+1}-A_{t}>0$

|(s,c)|(0,0)|(0,1)|(1,0)|(1,1)|(2,0)|(2,1)|(3,0)|(3,1)| 
|-----|-----|-----|-----|-----|-----|-----|-----|-----|
|(0,0)|  B  | S00 | C01 | N02 | C03 | N04 | C05 | N06 |
|(0,1)|     |  B  | N07 | C08 | N09 | C10 | N11 | C12 |
|(1,0)|     |     |  B  | S13 | C14 | N15 | C16 | N17 |
|(1,1)|     |     |     |  B  | N18 | C19 | N20 | C21 |
|(2,0)|     |     |     |     |  B  | S22 | C23 | N24 |
|(2,1)|     |     |     |     |     |  B  | N25 | C26 |
|(3,0)|     |     |     |     |     |     |  B  | S27 |
|(3,1)|     |     |     |     |     |     |     |  B  |

#### iii) GLM = 3: GLM1 + Preparation
iii-1) onset: trial onset

iii-2) duration: 0.1s

- 1: (0,0)
- 2: (0,1)
- 3: (1,0)
- 4: (1,1)
- 5: (2,0)
- 6: (2,1)
- 7: (3,0)
- 8: (3,1)
- 1+8: p(0,0)
- 2+8: p(0,1)
- 3+8: p(1,0)
- 4+8: p(1,1)
- 5+8: p(2,0)
- 6+8: p(2,1)
- 7+8: p(3,0)
- 8+8: p(3,1)

