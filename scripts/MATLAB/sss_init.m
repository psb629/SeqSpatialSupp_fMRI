if ispc
    rootDir='F:\SeqSpatialSupp_fMRI';
    dir_git = '\\wsl.localhost/ubuntu-24.04/home/sungbeenpark/github';
elseif ismac
    rootDir = '/Volumes/Diedrichsen_data$/data/SeqSpatialSupp_fMRI';
    dir_git = '/Users/sungbeenpark/github';
    addpath(genpath('/Users/sungbeenpark/SPM')); % nanmean (in surf_vol2surf)
    % SPM
    addpath('/Users/sungbeenpark/spm');
    % FreeSurfer
    setenv('FREESURFER_HOME', '/Applications/freesurfer/8.1.0');
    setenv('PATH', [getenv('FREESURFER_HOME') '/bin:' getenv('PATH')]);
    % WorkBench
    setenv('PATH', ['/Users/sungbeenpark/workbench/bin_macosxub:' getenv('PATH')]);
else
    fprintf('Workdir not found. Mount or connect to server and try again.');
end

baseDir         = rootDir;     % Base directory of the project
BIDSDir         = 'BIDS';                       % Raw data post AutoBids conversion
behavDir        = 'behavDir';           % Timing data from the scanner
imagingRawDir   = 'imaging_data_raw';           % Temporary directory for raw functional data
imagingDir      = 'imaging_data';               % Preprocesses functional data
anatomicalDir   = 'anatomicals';                % Preprocessed anatomicalcentr data (LPI + center AC + segemnt)
fmapDir         = 'fieldmaps';                  % Fieldmap dir after moving from BIDS and SPM make fieldmap
suitDir         = 'suit';
fsDir           = 'Freesurfer';
wbDir           = 'surfaceWB';
roiDir          = 'ROI';

if exist(dir_git, 'dir')
    addpath(genpath(dir_git));
end
dir_SSS = fullfile(dir_git,'SeqSpatialSupp_fMRI');
atlasDir = fullfile(dir_SSS,'atlas');