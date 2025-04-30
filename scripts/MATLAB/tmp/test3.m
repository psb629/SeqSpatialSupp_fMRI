workdir='/srv/diedrichsen/data/SeqSpatialSupp_fMRI';
baseDir         = (sprintf('%s/',workdir));     % Base directory of the project
BIDSDir        = 'BIDS';                       % Raw data post AutoBids conversion
behavDir        = 'behavDir';           % Timing data from the scanner
imagingRawDir   = 'imaging_data_raw';           % Temporary directory for raw functional data
imagingDir      = 'imaging_data';               % Preprocesses functional data
anatomicalDir   = 'anatomicals';                % Preprocessed anatomicalcentr data (LPI + center AC + segemnt)
fmapDir         = 'fieldmaps';                  % Fieldmap dir after moving from BIDS and SPM make fieldmap
suitDir         = 'suit';
regDir          = 'RegionOfInterest';
freesurferDir    = 'freesurf';
wbDir = 'surfaceWB';  %% standard surface?
glmDir = '/glm_%d';
roiDir = 'ROI';

sn = [1 2 3 5 6];
glm = 1;
for s=sn % for each subj
    fprintf('\nSubject: s%02d\n', s) % output to user

    % change directory to subject glm
    cd(fullfile(baseDir,sprintf(glmDir,glm),sprintf('S%02d',s)))
    %temp = dir('psc*nii');
    O = {};           
    O = {'psc_Motor2-1.nii','psc_Cue2-1.nii','psc_BothLetter2-1.nii',...
        'psc_BothSpatial2-1.nii','psc_CueLetter2-1.nii','psc_CueSpatial2-1.nii',...
        'psc_NRepMotor2-1.nii','psc_NRepCue2-1.nii','psc_NRep2-1.nii',...
        'psc_Letter.nii','psc_Spatial.nii'};
    O = {'pcs_Motor2-1.nii'};
    % load ROI
    load(fullfile(baseDir,roiDir,sprintf('%s_%s_regions.mat',sprintf('S%02d',s),'SSS')));

    %cond = [1 2 3 1 2 3]';
    cond = [1 2 3 4 5 6 7 8 9 10 11]';

    V=spm_vol(char(O));            
    % get raw data for voxels in region
    for r=1:length(R) % for each region

        Y=region_getdata(V,R{r});  
        S.psc=mean(Y,2);

        S.hemi=repmat(hemi,length(O),1);
        S.roi=repmat(r,length(O),1);
        S.SN=repmat(s,length(O),1);
        S.cond=cond;

        T=addstruct(T,S);
        fprintf('%d.',r)
    end
end