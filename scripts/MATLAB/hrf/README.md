# Region: Lightweight and simple Region toolbox for Matlab/SPM

## Puropose

1. Find appropriate HRF parameters for each ROI.
2. Utilize CifTi

## Function

1. `region`: creates region structure
2. `region_calcregion`: calculates the locations in regions and stores them as R.data
3. `region_getdata`: gets the values from a series of Image files for a series of regions
4. `region_getts`: gets Raw, predicted, adjusted and residual time series from series of regions
```
for rr = [1:8]
	idx = SPM.xX.K(rr).row;
	k = SPM.xX.K(rr).X0;
	w = SPM.xX.W(idx,idx);
	y_filt(idx,1) = w*y_raw(idx,1) - k*k'*w*y_raw(idx,1);
end

B = SPM.xX.pKX*y_filt;
y_hat = SPM.xX.xKXs.X(:,idx_interest)*B(idx_interest,:);
y_res = spm_sp('r', SPM.xX.xKXs, y_filt);
y_adj = y_hat + y_res;
```
5. `region_getirf`: extracts the ts (using `getts`) and gets evoked response for all events
6. `region_saveasimg`: saves a certain region as an image
7. `region_deformation`: deforms regions into individual space over a non-linear transformation

## Example of usage

1. make `ROIs.gii` file
```
sss_hrf('ROI:findall')
```

2. make `SubjID_Task_regions.glm_??.mat` file (using `region_calcregion` function)
```
sss_hrf('ROI:redefine')
```

3. (Using `region_getts` function)
```
sss_hrf('HRF:ROI_hrf_get')
```

4. (Using `region_getdata` and `spmj_fit_hrfparams` function)
```
sss_hrf('HRF:fit')
```
