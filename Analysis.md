## `GLM:all`

### `GLM:design`: `SeqSpatialSupp_fMRI/glm_<num>/<subj_id>`
$\leftarrow$ `behav_info.tsv`

`behav_info.tsv`를 load하여 *SPM.mat*과 reginfo.tsv*를 생성함 

### `GLM:estimate`: `SeqSpatialSupp_fMRI/glm_<num>/<subj_id>`
$\leftarrow$ `GLM:design`

Whole brain에 대하여 default HRF parameter 인 **[5,15,1,1,6,0,32]**로 GLM을 실행하고, $\beta$ 값을 계산함.

### `WB:vol2surf`: `SeqSpatialSupp_fMRI/surfaceWB/glm_<num>`
$\leftarrow$ `GLM:estimate`

계산한 *$\beta$.nii* 및 *ResMS.nii*를 2D surface에 mapping 함.

### `GLM:HRF_tuner`: `SeqSpatialSupp_fMRI/glm_<num>/<subj_id>/hrf_tune`
$\leftarrow$ `GLM:estimate`, `ROI:make_cifti.y_raw`

SPM.mat 파일을 불러들여 `SPM = spmj_glm_convolve(SPM)`를 통해 새로운 HRF parameter 를 적용한 후, ROI별 $y_{raw}$ CIFTI 파일을 불러들여 국소적인 GLM 을 계산하고 (`[beta, Yhat, Yres] = spmj_glm_fit(SPM,Yraw)`) 저장한다.

---

## `ROI:all`

### `ROI:calc_region`: `SeqSpatialSupp_fMRI/ROI/<subj_id>`
$\leftarrow$ `GLM:estimate`

먼저, `FreeSurfer`의 결과인 *white.surf.gii*와 *pial.surf.gii* 파일을 바탕으로 2D template surface *ROI.32k.L.label.gii* 파일과 3D original volumn space *mask.nii*을 연결하는 *<subj_id>.Task_regions.mat* 파일을 생성한다.

그 후 quality 체크를 위해, 생성한 <subj_id>.Task_regions.mat를 토대로 mask.nii space에 맞춘 *ROI.nii*를 각 부위별로 생성한다. 

$\rightarrow$ 생성한 ROI.nii 들은 anatomical.nii 와 함께 불러와서 align 을 체크해야한다!

cf) **anatomical.nii** 가 아닌 **mask.nii** 를 쓰는 이유는 전처리 과정에서 EPI image 이미지가 anatomical image에 align이 되기도 했고, $y_{raw}$나 $\beta$와 같은 데이터를 추출하는게 주 목적이기 때문이다.

### `ROI:make_cifti.y_raw`: `SeqSpatialSupp_fMRI/ROI/<subj_id>`
$\leftarrow$ `ROI:calc_region`, `GLM:estimate`

GLM 결과인 *SPM.mat* 에 저장된 3+1D whole brain $y_{raw}$ (=*SPM.xY.VY*) 를 <subj_id>.Task_regions.mat 에 저장된 ROI 정보를 토대로 2D surface로 추출하여 *cifti.<hemisphere>.<subj_id>.<ROI>.y_raw.dtseries.nii* 꼴(CIFTI)로 저장한다.

---
