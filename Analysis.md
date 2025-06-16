## `GLM:all`

### (1) `GLM:design`: `Diedrichsen_data$/data/SeqSpatialSupp_fMRI/glm_<num>/<subj_id>`
$\leftarrow$ `behav_info.tsv`

`behav_info.tsv`를 load하여 *SPM.mat*과 *reginfo.tsv*를 생성함 

### (2) `GLM:estimate`: `Diedrichsen_data$/data/SeqSpatialSupp_fMRI/glm_<num>/<subj_id>`
$\leftarrow$ `GLM:design`

Whole brain에 대하여 default HRF parameter 인 [5,15,1,1,6,0,32]로 GLM을 실행하고, $\beta$ 값을 계산함.

### (optional) `WB:vol2surf`: `Diedrichsen_data$/data/SeqSpatialSupp_fMRI/surfaceWB/glm_<num>`
$\leftarrow$ `GLM:estimate`

계산한 $\beta$*.nii* 및 *ResMS.nii*를 2D surface에 mapping 함.

### (5) `GLM:HRF_tuner`: `Diedrichsen_data$/data/SeqSpatialSupp_fMRI/glm_<num>/<subj_id>/hrf_tune`
$\leftarrow$ `GLM:estimate`, `ROI:make_cifti.y_raw`

SPM.mat 파일을 불러들여 `SPM = spmj_glm_convolve(SPM)`를 통해 새로운 HRF parameter 를 적용한 후, ROI 별 $y_{raw}$ CIFTI 파일을 불러들여 국소적인 GLM 을 계산하고 (`[beta, Yhat, Yres] = spmj_glm_fit(SPM,Yraw)`) 결과물을 CIFTI 포멧으로 저장한다.

$\rightarrow$ 각 HRF parameter 별, GLM 모델인 $y_{hat}$ 이 잘 동작하는지 평가할 $R^{2}=1-\frac{RSS}{TSS}$ 값을 계산하고 적절한 HRF parameter 를 선택한다.

$\rightarrow$ 새롭게 얻은 ROI 별 beta는 (univariate) spatial prewhitening 을 위해 ResMS.nii 의 값이 필요하므로, ResMS.nii 역시 ROI 별로 `Diedrichsen_data$/data/SeqSpatialSupp_fMRI/ROI/glm_<num>`에 저장할 필요가 있다.

---

## (3) `ROI:init`

### `ROI:calc_region`: `Diedrichsen_data$/data/SeqSpatialSupp_fMRI/ROI/<subj_id>`
$\leftarrow$ `GLM:estimate`

먼저, `FreeSurfer`의 결과인 *white.surf.gii*와 *pial.surf.gii* 파일을 바탕으로 2D template surface *ROI.32k.L.label.gii* 파일과 3D original volumn space *mask.nii*을 연결하는 *<subj_id>.Task_regions.mat* 파일을 생성한다.

그 후 quality 체크를 위해, 생성한 <subj_id>.Task_regions.mat를 토대로 mask.nii space에 맞춘 *ROI.nii*를 각 부위별로 생성한다. 

$\rightarrow$ 생성한 ROI.nii 들은 anatomical.nii 와 함께 불러와서 align 을 체크해야한다!

cf) **anatomical.nii** 가 아닌 **mask.nii** 를 쓰는 이유는 전처리 과정에서 EPI image 이미지가 anatomical image에 align이 되기도 했고, $y_{raw}$나 $\beta$와 같은 EPI-related 데이터를 추출하는게 주 목적이기 때문이다.

### `ROI:make_cifti.y_raw`: `Diedrichsen_data$/data/SeqSpatialSupp_fMRI/ROI/<subj_id>`
$\leftarrow$ `ROI:calc_region`, `GLM:estimate`

GLM 결과인 *SPM.mat* 에 저장된 (3+1)D whole brain $y_{raw}$ (=*SPM.xY.VY*) 를 <subj_id>.Task_regions.mat 에 저장된 ROI 정보를 토대로 2D surface로 추출하여 *cifti.<hemi>.<subj_id>.<roi>.y_raw.dtseries.nii* 꼴(CIFTI)로 저장한다.

## (4) `ROI:glm`

### `ROI:make_cifti.ResMS`: `Diedrichsen_data$/data/SeqSpatialSupp_fMRI/ROI/glm_<num>`
$\leftarrow$ `ROI:calc_region`, `GLM:estimate`

GLM 결과인 residual variance image *ResMS.nii* 를 <subj_id>.Task_regions.mat 에 저장된 ROI 정보를 토대로 2D surface로 추출하여 *cifti.<hemi>.<subj_id>.<roi>.ResMS.dscalar.nii* 꼴(CIFTI)로 저장한다.

---
