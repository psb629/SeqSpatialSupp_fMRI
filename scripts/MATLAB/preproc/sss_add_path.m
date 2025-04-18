% After mounting the diedrichsen datashare on a mac computer.
if ismac
    workdir='/Volumes/Diedrichsen_data$/data/SeqSpatialSupp_fMRI';
    atlasDir = '/Volumes/Diedrichsen_data$/data/Atlas_templates';
% After mounting the diedrichsen datashare on the CBS server.
elseif isunix
    workdir='/srv/diedrichsen/data/SeqSpatialSupp_fMRI';
    atlasDir = '/srv/diedrichsen/data/Atlas_templates';
    standardmeshDir = '/home/ROBARTS/skim2764/imaging_tools/surfAnalysis/standard_mesh';
elseif ispc
    workdir='D:/milli/diedrichsenlab/SeqSpatialSupp_fMRI';
    atlasDir='D:/mobaxterm/sungbeenpark/github/diedrichsenlab/atlas';
    standardmeshDir = 'D:/mobaxterm/sungbeenpark/github/surfAnalysis/standard_mesh';
    
else
    fprintf('Workdir not found. Mount or connect to server and try again.');
end

baseDir         = (sprintf('%s/',workdir));     % Base directory of the project
BIDSDir         = 'BIDS';                       % Raw data post AutoBids conversion
behavDir        = 'behavDir';           % Timing data from the scanner
imagingRawDir   = 'imaging_data_raw';           % Temporary directory for raw functional data
imagingDir      = 'imaging_data';               % Preprocesses functional data
anatomicalDir   = 'anatomicals';                % Preprocessed anatomicalcentr data (LPI + center AC + segemnt)
fmapDir         = 'fieldmaps';                  % Fieldmap dir after moving from BIDS and SPM make fieldmap
suitDir         = 'suit';
regDir          = 'RegionOfInterest';
freesurferDir   = 'freesurf';
wbDir = 'surfaceWB';  %% standard surface?
glmDir = '/glm_%d';
roiDir = 'ROI';

% pathToSave = fullfile(baseDir,'Figures');
% 
% addpath(genpath('D:/mobaxterm/sungbeenpark/github'));
% addpath(genpath('/home/ROBARTS/skim2764/Documents/MATLAB/scripts'));
% addpath(genpath(fullfile(workdir,behavDir)));
% addpath(genpath(fullfile(workdir,BIDSDir)));
% addpath(genpath('/home/ROBARTS/skim2764/imaging_tools'));
% addpath(genpath('/srv/diedrichsen/matlab/imaging/surfing'));
% addpath(genpath('/srv/diedrichsen/matlab/imaging/freesurfer'));
% addpath(genpath('/srv/diedrichsen/matlab/spm12'));