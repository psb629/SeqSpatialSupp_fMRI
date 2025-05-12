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
reconall -s <subject id>
```

ii) reslice
```
reslice -s <subject id> -a <surface name>
```

### 1-2. Functional Image

```
sss_imana('PREP:FUNC-all','sn',<subject number>)
```

#### 1-2-i. field map (optional)

#### 1-2-ii. EPI

---

## 2. GLM 

```
sss_GLM('GLM:all','sn',<subject number>,'glm',<GLM number>)
```

### Note, GLM numbers

i) GLM = 2: Repetition

|  | trial ${}_{t-1}$ | trial ${}_{t}$ |
|---------|---------|---------|
| Both-Rep| $(i,j)$ | $(i,j)$ |
| Cue-Rep | $(i,\neg j)$ | $(i,j)$ |
| Seq-Rep | $(\neg i,j)$ | $(i,j)$ |
| NRep    | $(\neg i,\neg j)$ | $(i,j)$ |

trial state: $ (i,j) $

- i=0: 32451
- i=1: 35124
- i=2: 13254
- i=3: 14523$
- j=0: Letter
- j=1: Spatial

ii) GLM = 3: Trial State

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
