function varargout = sss_imana(what,varargin)
%%
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
    workdir='F:/SeqSpatialSupp_fMRI';
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
addpath(genpath('D:/mobaxterm/sungbeenpark/github'));
% addpath(genpath('/home/ROBARTS/skim2764/Documents/MATLAB/scripts'));
% addpath(genpath(fullfile(workdir,behavDir)));
% addpath(genpath(fullfile(workdir,BIDSDir)));
% addpath(genpath('/home/ROBARTS/skim2764/imaging_tools'));
% addpath(genpath('/srv/diedrichsen/matlab/imaging/surfing'));
% addpath(genpath('/srv/diedrichsen/matlab/imaging/freesurfer'));
% addpath(genpath('/srv/diedrichsen/matlab/spm12'));
%% subject info

% Read info from participants .tsv file 
pinfo = dload('D:/mobaxterm/sungbeenpark/github/diedrichsenlab/SeqSpatialSupp_fMRI/participants.tsv');


% ROI
hemi = {'lh','rh'}; % left & right hemi folder names/prefixes
hem = {'L', 'R'}; % hemisphere: 1=LH 2=RH
hname = {'CortexLeft', 'CortexRight'}; % 'CortexLeft', 'CortexRight', 'Cerebellum'


% Subject
subj_name = {
    'S01', 'S02', 's03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10', 'S11', ...
    'S12', 'S13', 'S14', 'S15', 'S16', 'S17', ... 
    };

ns = numel(subj_name);
subj_vec = zeros(1,ns);
for i = 1:ns; subj_vec(1,i) = str2double(subj_name{i}(2:3)); end


% Scanner
TR = 1000;
numDummys  = 8;  % dummy images at the start of each run (these are discarded)
numTRs = 410; %% Adjust this for experiment
endTR = 410; % S.Park: prevent the error in sss_imana('FUNC:realign_unwarp','sn',s,'rtm',0)

% nTRs = [410*ones(1,10) 401 406 410 404 410 410 385]; % For S11
% nTRs = [410*ones(1,16) 385]; % for S09

fsl=16;
fs=13;
fontname='Arial';

%% MAIN OPERATION 
switch(what)
    case 'PREP:step0'           
        % All preprocessing steps by just one go (AC coordinates (loc_AC) are prerequisite)
        % handling input args:
        sn = [];
        vararginoptions(varargin,{'sn'})
        if isempty(sn)
            error('PREP:step0 -> sn must be inputted to this function.')
        end
        
        % loop on subjects and preprocessing the data:
        for s = sn
            % BIDS transfers:
            sss_imana('BIDS:move_unzip_raw_anat','sn',s);
            sss_imana('BIDS:move_unzip_raw_func','sn',s);
            sss_imana('BIDS:move_unzip_raw_fmap','sn',s);
            % ANAT functions:
            sss_imana('ANAT:reslice_LPI','sn',s);
            %  after this step, find the AC coordinate and enter it into
            %  participants.tsv
        end
        fprintf('Manually align the mean bias corrected image and the anatomical\n')
    case 'PREP:step1'           
        % All preprocessing steps by just one go (AC coordinates (loc_AC) are prerequisite)
        % handling input args:
        sn = [];
        vararginoptions(varargin,{'sn'})
        if isempty(sn)
            error('PREP:step1 -> sn must be inputted to this function.')
        end
        
        % loop on subjects and preprocessing the data:
        for s = sn
            % ANAT functions:
            sss_imana('ANAT:centre_AC','sn',s);
            sss_imana('ANAT:segmentation','sn',s);
            
            % FUNC functions:
            sss_imana('FUNC:make_fmap','sn',s);
            sss_imana('FUNC:realign_unwarp','sn',s,'rtm',0);
            sss_imana('FUNC:move_realigned_images','sn',s,'prefix','u','rtm',0);
            sss_imana('FUNC:meanimage_bias_correction','sn',s,'prefix','u','rtm',0);
            % after this step, manual aligment of mean bias corrected image
            % and the anatomical is required.
        end
        fprintf('Manually align the mean bias corrected image and the anatomical\n')

    case 'PREP:step2'
        % handling input args:
        sn = [];
        vararginoptions(varargin,{'sn'})
        if isempty(sn)
            error('PREP:step2 -> sn must be inputted to this function.')
        end
        
        % loop on subjects and preprocessing the data:
        for s = sn
            % FUNC:
            sss_imana('FUNC:coreg','sn',s,'prefix','u','rtm',0);
            sss_imana('FUNC:make_samealign','sn',s);
            sss_imana('FUNC:make_maskImage','sn',s);
        end


    case 'BIDS:move_unzip_raw_anat'
        % Moves, unzips and renames anatomical images from BIDS directory
        % to anatomicalDir
        
        % handling input args:
        sn = [];
        vararginoptions(varargin,{'sn'})
        if isempty(sn)
            error('BIDS:move_unzip_raw_anat -> ''sn'' must be passed to this function.')
        end

        % path to the subj anat data:
        % anat_raw_path = fullfile(baseDir,BIDSDir,sprintf('sub-%s',sid),'ses-01','anat',[char(pinfo.AnatRawName(pinfo.sn==sn)) '.nii.gz']);
        % anat_raw_path = fullfile(baseDir,BIDSDir,sprintf('sub-S%.02d',sn),'anat',[char(pinfo.AnatRawName(pinfo.sn==sn)) '.nii.gz']);
        anat_raw_path = fullfile(baseDir,BIDSDir,['sub-' char(pinfo.subj_id(pinfo.sn==sn))],'anat',[char(pinfo.AnatRawName(pinfo.sn==sn)) '.nii.gz']);
        % fprintf(anat_raw_path);
        % destination path:
        output_folder = fullfile(baseDir,anatomicalDir,char(pinfo.subj_id(pinfo.sn==sn)));
        output_file = fullfile(output_folder,[char(pinfo.subj_id(pinfo.sn==sn)) '_anatomical_raw.nii.gz']);

        if ~exist(output_folder,"dir")
            mkdir(output_folder);
        end

        % copy file to destination:
        [status,msg] = copyfile(anat_raw_path,output_file);
        if ~status  
            error('ANAT:move_anatomical -> subj %d raw anatmoical was not moved from BIDS to the destenation:\n%s',sn,msg)
        end

        % unzip the .gz files to make usable for SPM:
        gunzip(output_file);

        % delete the compressed file:
        delete(output_file);

    case 'BIDS:move_unzip_raw_func'
        % Moves, unzips and renames raw functional (BOLD) images from BIDS
        % directory
        
        % handling input args:
        sn = [];
        vararginoptions(varargin,{'sn'})
        if isempty(sn)
            error('BIDS:move_unzip_raw_func -> ''sn'' must be passed to this function.')
        end
        
        % loop on sessions:
        % for sess = 1:pinfo.numSess(pinfo.sn==sn)
            
        % pull list of runs from the participant.tsv:
        % run_list = eval(['pinfo.runsSess',num2str(sess),'(pinfo.sn==',num2str(sn),')']);
        run_list = str2double(split(pinfo.runlist{sn},'.'));

        % loop on runs of sess:
        for i = 1:length(run_list)
            
            % pull functional raw name from the participant.tsv:
            FuncRawName_tmp = eval(['pinfo.FuncRawNameSess(pinfo.sn==',num2str(sn),')']);
            FuncRawName_tmp = [char(FuncRawName_tmp) '.nii.gz'];

            % add run number to the name of the file:
            FuncRawName_tmp = replace(FuncRawName_tmp,'XX',sprintf('%.02d',i));

            % path to the subj func data:
            func_raw_path = fullfile(baseDir,BIDSDir,sprintf('sub-%s',pinfo.subj_id{pinfo.sn==sn}),'func',FuncRawName_tmp);
    
            % destination path:
            output_folder = fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn)));
            output_file = fullfile(output_folder,[char(pinfo.subj_id(pinfo.sn==sn)) sprintf('_run_%.02d.nii.gz',run_list(i))]);
            
            if ~exist(output_folder,"dir")
                mkdir(output_folder);
            end
            
            % copy file to destination:
            [status,msg] = copyfile(func_raw_path,output_file);
            if ~status  
                error('FUNC:move_unzip_raw_func -> subj %d raw functional (BOLD) was not moved from BIDS to the destenation:\n%s',sn,msg)
            end
    
            % unzip the .gz files to make usable for SPM:
            gunzip(output_file);
    
            % delete the compressed file:
            delete(output_file);
        end
        % end    

    case 'BIDS:move_unzip_raw_fmap'
        % Moves, unzips and renames raw fmap images from BIDS
        % directory
        
        % handling input args:
        sn = [];
        vararginoptions(varargin,{'sn'})
        if isempty(sn)
            error('BIDS:move_unzip_raw_fmap -> ''sn'' must be passed to this function.')
        end
        
        % loop on sessions:
        for sess = 1:pinfo.numSess(pinfo.sn==sn)
            % pull fmap raw names from the participant.tsv:
            fmapMagnitudeName_tmp = eval(['pinfo.fmapMagnitudeNameSess(pinfo.sn==',num2str(sn),')']);
            magnitude = [char(fmapMagnitudeName_tmp) '.nii.gz'];
            
            fmapPhaseName_tmp = eval(['pinfo.fmapPhaseNameSess(pinfo.sn==',num2str(sn),')']);
            phase = [char(fmapPhaseName_tmp) '.nii.gz'];

            % path to the subj fmap data:
            magnitude_path = fullfile(baseDir,BIDSDir,sprintf('sub-%s',pinfo.subj_id{pinfo.sn==sn}),'fmap',magnitude);
            phase_path = fullfile(baseDir,BIDSDir,sprintf('sub-%s',pinfo.subj_id{pinfo.sn==sn}),'fmap',phase);
    
            % destination path:
            output_folder = fullfile(baseDir,fmapDir,char(pinfo.subj_id(pinfo.sn==sn)));
            output_magnitude = fullfile(output_folder,[char(pinfo.subj_id(pinfo.sn==sn)) '_magnitude.nii.gz']);
            output_phase = fullfile(output_folder,[char(pinfo.subj_id(pinfo.sn==sn)) '_phase.nii.gz']);
            
            if ~exist(output_folder,"dir")
                mkdir(output_folder);
            end
            
            % copy magnitude to destination:
            [status,msg] = copyfile(magnitude_path,output_magnitude);
            if ~status  
                error('BIDS:move_unzip_raw_fmap -> subj %d, fmap magnitude was not moved from BIDS to the destenation:\n%s',sn,msg)
            end
            % unzip the .gz files to make usable for SPM:
            gunzip(output_magnitude);
    
            % delete the compressed file:
            delete(output_magnitude);

            % copy phase to destination:
            [status,msg] = copyfile(phase_path,output_phase);
            if ~status  
                error('BIDS:move_unzip_raw_fmap -> subj %d, fmap phase was not moved from BIDS to the destenation:\n%s',sn,msg)
            end
            % unzip the .gz files to make usable for SPM:
            gunzip(output_phase);
    
            % delete the compressed file:
            delete(output_phase);
        end 
        
    case 'ANAT:reslice_LPI'          
        % Reslice anatomical image within LPI coordinate systems
        % handling input args:
        sn = [];
        vararginoptions(varargin,{'sn'})
        if isempty(sn)
            error('ANAT:reslice_LPI -> ''sn'' must be passed to this function.')
        end
        
        % (1) Reslice anatomical image to set it within LPI co-ordinate frames
        source = fullfile(baseDir,anatomicalDir,char(pinfo.subj_id(pinfo.sn==sn)),[char(pinfo.subj_id(pinfo.sn==sn)) '_anatomical_raw.nii']);
        dest = fullfile(baseDir,anatomicalDir,char(pinfo.subj_id(pinfo.sn==sn)),[char(pinfo.subj_id(pinfo.sn==sn)) '_anatomical.nii']);
        spmj_reslice_LPI(source,'name', dest);
        
        % (2) In the resliced image, set translation to zero (why?)
        V               = spm_vol(dest);
        dat             = spm_read_vols(V);
        V.mat(1:3,4)    = [0 0 0];
        spm_write_vol(V,dat);


    case 'ANAT:centre_AC'            
        % Description:
        % Recenters the anatomical data to the Anterior Commissure
        % coordiantes. Doing that, the [0,0,0] coordiante of subject's
        % anatomical image will be the Anterior Commissure.

        % You should manually find the voxel coordinates of AC 
        % for each from their anatomical scans and add it to the
        % participants.tsv file under the loc_ACx loc_ACy loc_ACz columns.

        % This function runs for all subjects and sessions.

        % handling input args:
        sn = [];
        vararginoptions(varargin,{'sn'})
        if isempty(sn)
            error('ANAT:centre_AC -> ''sn'' must be passed to this function.')
        end
        
        % path to the raw anatomical:
        anat_raw_file = fullfile(baseDir,anatomicalDir,char(pinfo.subj_id(pinfo.sn==sn)),[char(pinfo.subj_id(pinfo.sn==sn)) '_anatomical.nii']);
        if ~exist(anat_raw_file,"file")
            error('ANAT:centre_AC -> file %s was not found.',anat_raw_file)
        end
        
        % Get header info for the image:
        V = spm_vol(anat_raw_file);
        % Read the volume:
        dat = spm_read_vols(V);
        
        % changing the transform matrix translations to put AC near [0,0,0]
        % coordinates:
        % R = V.mat(1:3,1:3);
        % t = -1 * R * [pinfo.locACx(pinfo.sn==sn),pinfo.locACy(pinfo.sn==sn),pinfo.locACz(pinfo.sn==sn)]';
        % V.mat(1:3,4) = t;
        V.mat(1:3,4) = -[pinfo.locACx(pinfo.sn==sn),pinfo.locACy(pinfo.sn==sn),pinfo.locACz(pinfo.sn==sn)];
        % writing the image with the changed header:
        spm_write_vol(V,dat);
        fprintf(1,'Done\n');

    case 'ANAT:segmentation'
        % Segmentation + Normalization
        % manually check results when done

        % handling input args:
        sn = [];
        vararginoptions(varargin,{'sn'})
        if isempty(sn)
            error('ANAT:segmentation -> ''sn'' must be passed to this function.')
        end
        
        SPMhome=fileparts(which('spm.m'));
        J=[];
        for s=sn
            J.channel.vols = {fullfile(baseDir,anatomicalDir,char(pinfo.subj_id(pinfo.sn==sn)),[char(pinfo.subj_id(pinfo.sn==sn)),'_anatomical.nii,1'])};
            J.channel.biasreg = 0.001;
            J.channel.biasfwhm = 60;
            J.channel.write = [0 0];
            J.tissue(1).tpm = {fullfile(SPMhome,'tpm/TPM.nii,1')};
            J.tissue(1).ngaus = 1;
            J.tissue(1).native = [1 0];
            J.tissue(1).warped = [0 0];
            J.tissue(2).tpm = {fullfile(SPMhome,'tpm/TPM.nii,2')};
            J.tissue(2).ngaus = 1;
            J.tissue(2).native = [1 0];
            J.tissue(2).warped = [0 0];
            J.tissue(3).tpm = {fullfile(SPMhome,'tpm/TPM.nii,3')};
            J.tissue(3).ngaus = 2;
            J.tissue(3).native = [1 0];
            J.tissue(3).warped = [0 0];
            J.tissue(4).tpm = {fullfile(SPMhome,'tpm/TPM.nii,4')};
            J.tissue(4).ngaus = 3;
            J.tissue(4).native = [1 0];
            J.tissue(4).warped = [0 0];
            J.tissue(5).tpm = {fullfile(SPMhome,'tpm/TPM.nii,5')};
            J.tissue(5).ngaus = 4;
            J.tissue(5).native = [1 0];
            J.tissue(5).warped = [0 0];
            J.tissue(6).tpm = {fullfile(SPMhome,'tpm/TPM.nii,6')};
            J.tissue(6).ngaus = 2;
            J.tissue(6).native = [0 0];
            J.tissue(6).warped = [0 0];
            J.warp.mrf = 1;
            J.warp.cleanup = 1;
            J.warp.reg = [0 0.001 0.5 0.05 0.2];
            J.warp.affreg = 'mni';
            J.warp.fwhm = 0;
            J.warp.samp = 3;
            J.warp.write = [1 1]; %
            matlabbatch{1}.spm.spatial.preproc=J;
            spm_jobman('run',matlabbatch);
        end

    case 'FUNC:make_fmap'                
        % Description:
        % Generates VDM files from the presubtracted phase & magnitude
        % images acquired from the field map sequence. Also, just as a
        % quality control this function creates unwarped EPIs from the
        % functional data with the prefix 'u' for each run.
        
        % handling input args:
        sn = [];
        vararginoptions(varargin,{'sn'})
        if isempty(sn)
            error('FUNC:make_fmap -> ''sn'' must be passed to this function.')
        end
        
        % loop on sessions:
        % for sess = 1:pinfo.numSess(pinfo.sn==sn)
            % Prefix of the functional files:
        prefixepi  = '';
        % read this one
        % https://lcni.uoregon.edu/kb-articles/kb-0003

        % echo times of the gradient eho sequence:
        % total EPI readout time = = echo spacing (in ms) * base resolution 
        % (also knows as number of echos). If you use GRAPPA acceleration, 
        % you need to divide the total number of echos by two:
        [et1, et2, tert] = spmj_et1_et2_tert(workdir, char(pinfo.subj_id(pinfo.sn==sn)));




        % pull list of runs from the participant.tsv:
        % run_list = pinfo.runlist;

        % subfolderFieldmap = sprintf('sess%d',sess);
        subfolderFieldmap = '';
        % function to create the makefieldmap job and passing it to the SPM
        % job manager:
        run_list_cell = cellfun(@(x) sprintf('%.02d',str2double(x)), split(pinfo.runlist{sn},'.'), 'UniformOutput', false);
        
        spmj_makefieldmap(baseDir,char(pinfo.subj_id(pinfo.sn==sn)),run_list_cell, ...
                          'et1', et1, ...
                          'et2', et2, ...
                          'tert', tert, ...
                          'image', numDummys+1, ...
                          'prefix',prefixepi, ...
                          'rawdataDir',fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn))), ...
                          'subfolderFieldmap',subfolderFieldmap);
        % end

    case 'FUNC:realign_unwarp'      
        % Do spm_realign_unwarp
        
        % handling input args:
        sn = [];
        rtm = 0; % 0: first volume, 1: mean volume
        vararginoptions(varargin,{'sn','rtm','endTR'})
        if isempty(sn)
            error('FUNC:realign_unwarp -> ''sn'' must be passed to this function.')
        end

        % loop on sessions:
        % for sess = 1:pinfo.numSess(pinfo.sn==sn)
        % Prefix of the functional files:
        prefixepi  = '';
        %startTR = numDummys+1;
        startTR = 1;
        % pull list of runs from the participant.tsv:
        run_list_cell = cellfun(@(x) sprintf('%.02d',str2double(x)), split(pinfo.runlist{sn},'.'), 'UniformOutput', false);
        
        subfolderFieldmap = '';
%         endTR = [405*ones(1,16) 380];  %% For S11
        % below inf -> enTR
        spmj_realign_unwarp(baseDir,char(pinfo.subj_id(pinfo.sn==sn)),run_list_cell, startTR, endTR, ...
                            'prefix',prefixepi,...
                            'rawdataDir',fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn))),...
                            'subfolderFieldmap',subfolderFieldmap,...
                            'rtm',rtm);
        % end

    case 'FUNC:move_realigned_images'          
        % Move images created by realign(+unwarp) into imaging_data
        
        % handling input args:
        sn = [];
        prefix = 'u';   % prefix of the 4D images after realign(+unwarp)
        rtm = 0;        % realign_unwarp registered to the first volume (0) or mean image (1).
        vararginoptions(varargin,{'sn','prefix','rtm'})
        if isempty(sn)
            error('FUNC:move_realigned_images -> ''sn'' must be passed to this function.')
        end

        % loop on sessions:
        % for sess = 1:pinfo.numSess(pinfo.sn==sn)
        % pull list of runs from the participant.tsv:
        run_list_cell = cellfun(@(x) sprintf('%.02d',str2double(x)), split(pinfo.runlist{sn},'.'), 'UniformOutput', false);
        
        % loop on runs of the session:
        for r = 1:length(run_list_cell)
            % realigned (and unwarped) images names:
            file_name = [prefix, char(pinfo.subj_id(pinfo.sn==sn)), '_run_', run_list_cell{r}, '.nii'];
            source = fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn)),file_name);
            dest = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)));
            if ~exist(dest,'dir')
                mkdir(dest)
            end
            dest = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),file_name);
            % move to destination:
            [status,msg] = movefile(source,dest);
            if ~status  
                error('BIDS:move_realigned_images -> %s',msg)
            end

            % realign parameters names:
            source = fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn)),['rp_', char(pinfo.subj_id(pinfo.sn==sn)), '_run_', run_list_cell{r}, '.txt']);
            dest = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),['rp_', char(pinfo.subj_id(pinfo.sn==sn)), '_run_', run_list_cell{r}, '.txt']);
            % move to destination:
            [status,msg] = movefile(source,dest);
            if ~status  
                error('BIDS:move_realigned_images -> %s',msg)
            end
        end
        
        % mean epi name - the generated file name will be different for
        % rtm=0 and rtm=1. Extra note: rtm is an option in
        % realign_unwarp function. Refer to spmj_realign_unwarp().
        if rtm==0   % if registered to first volume of each run:
            source = fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn)),['mean', prefix, char(pinfo.subj_id(pinfo.sn==sn)), '_run_', run_list_cell{1}, '.nii']);
            dest = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),['mean', prefix, char(pinfo.subj_id(pinfo.sn==sn)), '_run_', run_list_cell{1}, '.nii']);
        else        % if registered to mean image of each run:
            source = fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn)),[prefix, 'meanepi_', char(pinfo.subj_id(pinfo.sn==sn)), '.nii']);
            dest = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),[prefix, 'meanepi_', char(pinfo.subj_id(pinfo.sn==sn)), '.nii']);
        end
        % move to destination:
        [status,msg] = movefile(source,dest);
        if ~status  
            error('BIDS:move_realigned_images -> %s',msg)
        end
        % end
    
    case 'FUNC:meanimage_bias_correction'                                         
        % correct bias for mean image of 1st session. If the realignment
        % was done with respect to the first volume of each run of each 
        % session, the mean image will be calculated on the first run of
        % each session and will be called 'meanu*_run_01.nii' ('mean' 
        % indicates the image is average of the volumes and 'u' indicates
        % it's unwarped). Therefore, we do the bias correction on this 
        % file.
        % But if you do the realignment to the mean epi of every run, the
        % generated mean file will be named 'umeanepi_*' and we do the bias
        % correction on this file.

        % handling input args:
        sn = [];
        prefix = 'u';   % prefix of the 4D images after realign(+unwarp)
        rtm = 0;        % realign_unwarp registered to the first volume (0) or mean image (1).
        vararginoptions(varargin,{'sn','prefix','rtm'})
        if isempty(sn)
            error('FUNC:meanimage_bias_correction -> ''sn'' must be passed to this function.')
        end

        % loop on sessions:
        % for sess = 1:pinfo.numSess(pinfo.sn==sn)
        % pull list of runs from the participant.tsv:
        run_list_cell = cellfun(@(x) sprintf('%.02d',str2double(x)), split(pinfo.runlist{sn},'.'), 'UniformOutput', false);


        if rtm==0   % if registered to first volume of each run:
            P{1} = fullfile(baseDir, imagingDir, char(pinfo.subj_id(pinfo.sn==sn)),['mean', prefix,  char(pinfo.subj_id(pinfo.sn==sn)), '_run_', run_list_cell{1}, '.nii']);
        else        % if registered to mean image of each run:
            P{1} = fullfile(baseDir, imagingDir, char(pinfo.subj_id(pinfo.sn==sn)),[prefix, 'meanepi_', char(pinfo.subj_id(pinfo.sn==sn)), '.nii']);
        end
        spmj_bias_correct(P);
        % end


    case 'FUNC:coreg'                                                      
        % coregister rbumean image to anatomical image for each session
        
        % handling input args:
        sn = [];
        prefix = 'u';   % prefix of the 4D images after realign(+unwarp)
        rtm = 0;        % realign_unwarp registered to the first volume (0) or mean image (1).
        vararginoptions(varargin,{'sn','prefix','rtm'})
        if isempty(sn)
            error('FUNC:coreg -> ''sn'' must be passed to this function.')
        end
        
        % loop on sessions:
        % for sess = 1:pinfo.numSess(pinfo.sn==sn)
        % (1) Manually seed the functional/anatomical registration
        % - Open fsleyes
        % - Add anatomical image and b*mean*.nii (bias corrected mean) image to overlay
        % - click on the bias corrected mean image in the 'Overlay
        %   list' in the bottom left of the fsleyes window.
        %   list to highlight it.
        % - Open tools -> Nudge
        % - Manually adjust b*mean*.nii image to the anatomical by 
        %   changing the 6 paramters (tranlation xyz and rotation xyz) 
        %   and Do not change the scales! 
        % - When done, click apply and close the tool tab. Then to save
        %   the changes, click on the save icon next to the mean image 
        %   name in the 'Overlay list' and save the new image by adding
        %   'r' in the beginning of the name: rb*mean*.nii. If you don't
        %   set the format to be .nii, fsleyes automatically saves it as
        %   a .nii.gz so either set it or gunzip afterwards to make it
        %   compatible with SPM.
        
        % (2) Run automated co-registration to register bias-corrected meanimage to anatomical image
        
        % pull list of runs from the participant.tsv:
        run_list_cell = cellfun(@(x) sprintf('%.02d',str2double(x)), split(pinfo.runlist{sn},'.'), 'UniformOutput', false);


        if rtm==0   % if registered to first volume
            % when you run realign and unwarp separately
            mean_file_name = sprintf('rbmean%s%s_run_%s.nii', prefix, char(pinfo.subj_id(pinfo.sn==sn)), run_list_cell{1}); 
            % when you run realign and unwarp simultaneously (after FUNC:realign_unwarp)
%           mean_fie_name = sprintf('bmean%s%s_run_%s.nii', prefix, char(pinfo.subj_id(pinfo.sn==sn)), run_list_cell{1});

        else    % if registered to the mean image
            % when you run realign and unwarp separately
            mean_file_name = sprintf('rb%smeanepi_%s.nii', prefix, char(pinfo.subj_id(pinfo.sn==sn)));
            % when you run realign and unwarp simultaneously (after FUNC:realign_unwarp)
%             mean_file_name = sprintf('b%smeanepi_%s.nii', prefix, char(pinfo.subj_id(pinfo.sn==sn)));

        end
        J.source = {fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),mean_file_name)}; 
        J.ref = {fullfile(baseDir,anatomicalDir,char(pinfo.subj_id(pinfo.sn==sn)),[char(pinfo.subj_id(pinfo.sn==sn)), '_anatomical','.nii'])};
        J.other = {''};
        % J.eoptions.cost_fun = 'ncc';
        J.eoptions.cost_fun = 'nmi';

        J.eoptions.sep = [4 2];
        J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        J.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estimate=J;
        spm_jobman('run',matlabbatch);
        
        % (3) Check alignment manually by using fsleyes similar to step
        % one.
        % end
    case 'FUNC:make_samealign'
        % align to registered bias corrected mean image of each session
        % (rb*mean*.nii). Alignment happens only by changing the transform
        % matrix in the header files of the functional 4D .nii files to the
        % transform matrix that aligns them to anatomical. The reason that
        % it works is: 1) in the realignment (+unwarping) process, we have
        % registered every single volume of every single run to the first
        % volume of the first run of the session. 2) In the same step, for
        % each session, a mean functional image (meanepi*.nii or meanu*.nii
        % based on the rtm option) was generated. This mean image is alread
        % in the space of all the functional volumes. Later we coregister
        % this image to the anatomical space. Therefore, if we change the
        % transformation matrices of all the functional volumes to the
        % transform matrix of the coregistered image, they will all
        % tranform into the anatomical coordinates space.

        % handling input args:
        sn = [];
        prefix = 'u';   % prefix of the 4D images after realign(+unwarp)
        rtm = 0;        % realign_unwarp registered to the first volume (0) or mean image (1).
        vararginoptions(varargin,{'sn','prefix','rtm'})
        if isempty(sn)
            error('FUNC:make_samealign -> ''sn'' must be passed to this function.')
        end

        % loop on sessions:
        % for sess = 1:pinfo.numSess(pinfo.sn==sn)
        % pull list of runs from the participant.tsv:
         run_list_cell = cellfun(@(x) sprintf('%.02d',str2double(x)), split(pinfo.runlist{sn},'.'), 'UniformOutput', false);

        
        % select the reference image:
        if rtm==0
            P{1} = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),['rb' 'mean' prefix char(pinfo.subj_id(pinfo.sn==sn)) '_run_' run_list_cell{1} '.nii']);
        else
            P{1} = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),['rb' prefix 'meanepi_' char(pinfo.subj_id(pinfo.sn==sn)) '.nii']);
        end

        % select images to be realigned:
        Q = {};
        for r = 1:length(run_list_cell)
            for i = 1:numTRs
                 Q{end+1} = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),[prefix char(pinfo.subj_id(pinfo.sn==sn)) '_run_' run_list_cell{r} '.nii,' num2str(i)]);
            end
        end

        spmj_makesamealign_nifti(char(P),char(Q));
        % end
    
    case 'FUNC:make_maskImage'       
        % Make mask images (noskull and gray_only) for 1st level glm
        
        % handling input args:
        sn = [];
        prefix = 'u';   % prefix of the 4D images after realign(+unwarp)
        rtm = 0;        % realign_unwarp registered to the first volume (0) or mean image (1).
        vararginoptions(varargin,{'sn','prefix','rtm'})
        if isempty(sn)
            error('FUNC:check_samealign -> ''sn'' must be passed to this function.')
        end

        % loop on sessions:
        % for sess = 1:pinfo.numSess(pinfo.sn==sn)
        run_list_cell = cellfun(@(x) sprintf('%.02d',str2double(x)), split(pinfo.runlist{sn},'.'), 'UniformOutput', false);

        subj_id = char(pinfo.subj_id(pinfo.sn==sn));
        % bias corrected mean epi image:
        if rtm==0
            nam{1} = fullfile(baseDir,imagingDir,subj_id,['rb' 'mean' prefix subj_id '_run_' run_list_cell{1} '.nii']);
        else
            nam{1} = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),['rb' prefix 'meanepi_' char(pinfo.subj_id(pinfo.sn==sn)) '.nii']);
        end
        nam{2}  = fullfile(baseDir,anatomicalDir, subj_id, ['c1' subj_id, '_anatomical.nii']);
        nam{3}  = fullfile(baseDir,anatomicalDir, subj_id, ['c2' subj_id, '_anatomical.nii']);
        nam{4}  = fullfile(baseDir,anatomicalDir, subj_id, ['c3' subj_id, '_anatomical.nii']);
        spm_imcalc(nam, fullfile(baseDir,imagingDir,subj_id,'rmask_noskull.nii'), 'i1>1 & (i2+i3+i4)>0.2')
        
        source = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)), 'rmask_noskull.nii');
        dest = fullfile(baseDir,anatomicalDir,char(pinfo.subj_id(pinfo.sn==sn)),'rmask_noskull.nii');
        movefile(source,dest);
        
        % gray matter mask for covariance estimation
        % ------------------------------------------
        nam={};
        nam{1}  = fullfile(baseDir, imagingDir,subj_id, ['rbmean' prefix subj_id '_run_' run_list_cell{1} '.nii']);
        nam{2}  = fullfile(baseDir, anatomicalDir, subj_id, ['c1' subj_id, '_anatomical.nii']);
        spm_imcalc(nam, fullfile(baseDir,imagingDir,subj_id,'rmask_gray.nii'), 'i1>1 & i2>0.4')
        
        source = fullfile(baseDir,imagingDir,subj_id, 'rmask_gray.nii');
        dest = fullfile(baseDir, anatomicalDir,subj_id,'rmask_gray.nii');
        movefile(source,dest);
        % end      

    case 'GLM:design'
        % handling input args:
        sn = [];
        prefix = 'u';
        % TODOfreesurferDir
        hrf_params = [6 16];
        hrf_cutoff = 128;
        % cvi_type = 'wls';
        cvi_type = 'fast';
        % nTR = 410; %% modify when necessary
        vararginoptions(varargin,{'sn','glm','nTR','hrf_params'});   
        % def_params = spm_get_defaults('stats.fmri.hrf');
        % def_params([1:2]) = hrf_params; hrf_params=def_params;
        % run_list = str2double(split(pinfo.runlist(pinfo.sn==sn),'.'));
        run_list = [1:8];
        % Save the aux. information file (SPM_info.mat).
        % This file contains user-friendly information about the glm
        % model, regressor types, condition names, etc.
        T = struct();
        % initialize 
        J = struct();
        % load preprocessed behavioral data
        %if glm~=0
        % R = construct_dsgmat(sn, glm); %% Important: index of R indicates a trial number
        %end
        % R1 = construct_dsgmat(sprintf('S%02d',sn),glm);
        % R2 = construct_dsgmat(sprintf('R%02d',sn),glm);
        % R = combine_behavdata(R1,R2);
        subj_id = char(pinfo.subj_id(pinfo.sn==sn));
        behav_data = fullfile(workdir,behavDir,sprintf('sub-%s/ssh__%s.dat',subj_id,subj_id));
        R = construct_dsgmat(behav_data,glm);
        
        J.dir = {fullfile(baseDir, sprintf(glmDir, glm), subj_id)};
        J.timing.units = 'secs'; % timing unit that all timing in model will be
        J.timing.RT = 1; % TR (in seconds, as per 'J.timing.units')
        J.timing.fmri_t = 16;
        J.timing.fmri_t0 = 1;
        % make GLM dir
        if (~exist(J.dir{1},'dir'))
            mkdir(J.dir{1});
        end
        % loop through runs within the current sessions
        % itaskUni = 0;
        if glm==0  %% for 9th run, slow-event design
           run_list = 17; % use run 17
           % run_list = 9; % use run 9
           % nTR = 385; % number of TRs for run 09
           R.cond = ones(1,nTR);
           
        end
        n_cond = length(unique(R.cond));
        if sn==9 % for S09
            nTR = [410*ones(1,8) 385];
        elseif sn==11 % for S11
            nTR = [410*ones(1,8) 385];
        elseif sn==25 % for R11
            nTR = [410*ones(1,2) 401 406 410 404 410 410];
        else
            nTR = [410*ones(1,8)];
        end
        
        % if length(run_list)==9
        %    nTR = [nTR 385];
        % end
        for r=run_list
            % get functional runs
            for i=1:nTR(r)
                N{i} = fullfile(baseDir,imagingDir,subj_id,[prefix sprintf('%s_run_%02d',subj_id,run_list(r)) '.nii,',num2str(i)]);
            end                
            J.sess(r).scans= N;
            if glm==0
                cond_name = {'All trials'};  %% only for 9th run
            elseif glm==1
                for c=1:n_cond cond_name{c} = sprintf('Trial %d',c);end
            elseif glm==2 
                cond_name = {'MotorOnly-L','MotorOnly-S','CueOnly-L','CueOnly-S',...
                            'BothRep-L','BothRep-S','NonRep-L','NonRep-S','Non-Interest'};
            elseif glm==3
                for c=1:n_cond cond_name{c} = sprintf('Trial-State %d',c);end;
                cond_name{n_cond+1} = 'Non-Interest';
            elseif glm==4
                cond_name = {'Letter','Spatial','Non-Interest'};
            end

            for c=1:n_cond  %% c : condition index

                % filling in "reginfo"
                % Transition ID, TT.task_name: 1~68, transition: (seqID, seqType), (0,0)->(0,0):1, (0,0)->(0,1):2,
                % (0,1)->(0,0):3, (0,1)->(0,1):4,...(3,1)->(3,1):64, 
                % 1st trial:65, 18th: 66, 35th: 67, 52nd: 68 
%                 % itaskUni = itaskUni + 1;
%                 TT.sn        = sn;
%                 TT.sess      = ses;
%                 TT.run       = r;
%                 TT.tS = R.tS(r,t);
%                 TT.transID(r,t) = R.transID(r,t);
%                 TT.isRepMotor(r,t) = R.isRepMotor(r,t);
%                 TT.isRepCue(r,t) = R.isRepCue(r,t);
%                 TT.isRepBoth(r,t) = R.isRepBoth(r,t);
%                 TT.isNrep(r,t) = R.isNrep(r,t);
%                 TT.isNint(r,t) = R.isNint(r,t);

%                 TT.task_name = sprintf('Transition I/D: %d',c);
%                 if glm==1
%                     TT.transId      = R.cond(r,t);  %% Transition ID
%                 else 
%                     TT.repIdx = R.cond(r,t); % 1: rep, -1: non-rep, 0: non-interest
%                 end
                % TT.taskUni   = itaskUni;
                % 
                % J.sess(r).cond(c).name = sprintf('transition %d',c);
                % J.sess(r).cond(c).onset = R.onset(r,find(R.cond(r,:)==c));
                % J.sess(r).cond(c).duration = R.dur(r,find(R.cond(r,:)==c));
                if glm==0
                    J.sess(1).cond(1).name = 'All trials';
                    J.sess(1).cond(1).onset = R.onset'+1;  %% added 1
                    J.sess(1).cond(1).duration = 0.001;
                elseif glm==1
                    J.sess(r).cond(c).name = sprintf('trial %d',c);
                    J.sess(r).cond(c).onset = R.onset(r,c);
                    J.sess(r).cond(c).duration = 0.001; % used fixed time, 2 secon
                    % J.sess(r).cond(c).duration = R.dur(r,c);
                else
                    J.sess(r).cond(c).name = cond_name{c};
                    J.sess(r).cond(c).onset = R.onset(r,find(R.cond(r,:)==c));
                    J.sess(r).cond(c).duration = 0.001; % used fixed time, 2 secon
%                     J.sess(r).cond(c).duration = R.dur(r,find(R.cond(r,:)==c));
                end

                J.sess(r).cond(c).tmod=0;
                J.sess(r).cond(c).orth=0;
                J.sess(r).cond(c).pmod=struct('name',{},'param',{},'poly',{});

               % add the condition info to the reginfo structure

            end
            % add any additional regressors here.
            J.sess(r).multi={''};
            J.sess(r).regress=struct('name',{},'val',{});
            J.sess(r).multi_reg={''};
            
            % Define high pass filter cutoff (in seconds): see glm cases.
            J.sess(r).hpf=hrf_cutoff;
        end

        %T = addstruct(T, TT);

        J.fact = struct('name', {}, 'levels', {});
        J.bases.hrf.derivs = [0 0];
        J.bases.hrf.params = hrf_params;
        % defaults.stats.fmri.hrf([1:2])=J.bases.hrf.params; 
        J.volt = 1;
        J.global = 'None';
        J.mask = {fullfile(baseDir,anatomicalDir, subj_id, 'rmask_noskull.nii,1')};
        J.mthresh = 0.05;
        J.cvi_mask = {fullfile(baseDir, anatomicalDir, subj_id, 'rmask_gray.nii')};
        J.cvi = cvi_type;
        
        % Save the GLM file for this subject.
        spm_rwls_run_fmri_spec(J);
        tmp = load(fullfile(J.dir{1},'SPM.mat'));
        delete(fullfile(J.dir{1},'SPM.mat'));
        SPM = tmp.SPM;
        save(fullfile(J.dir{1},'SPM.mat'),'SPM','-v7.3');
        save(fullfile(J.dir{1},'SPM_info.mat'),'R');
    
    case 'GLM:estimate' % estimate beta coefficient
        % Estimate the GLM from the appropriate SPM.mat file.
        % Make GLM files with case 'GLM_make'.
        sn = [];
        fig = 0;
        vararginoptions(varargin, {'sn', 'glm','fig'});
        subj_id = char(pinfo.subj_id(pinfo.sn==sn));       
        load(fullfile(baseDir,sprintf(glmDir, glm), subj_id, 'SPM.mat'));
        SPM.swd = fullfile(baseDir,sprintf(glmDir, glm), subj_id);
        % Run the GLM.
        if fig==0
            spm_rwls_spm(SPM);
        else
            dm = SPM.xX.xKXs.X(SPM.Sess(1).row,SPM.Sess(1).col); % design matrix for one run
            figure; imagesc(dm); axis square; colorbar;
            title('Design matrix'); xlabel('Regressors'); ylabel('Volumes'); set(gca, 'fontsize', fs)
            cv = cov(dm);                   % covariance matrix for one run
            varE = nanmean(diag(inv(cv)));  % mean variance of the regression estimates
            varX = nanmean(1./diag(cv));    % mean variance of the regressors in the design matrix
            vif = varE ./ varX;             % variance inflation factor
            figure;
            subplot(2,2,1)
            imagesc(cv); axis square; title('Covariance matrix'); colorbar;
            subplot(2,2,2)
            imagesc(inv(cv)); axis square; title('Inverse of the covariance'); colorbar;
            subplot(2,2,3)
            imagesc(1./diag(cv)); axis square; title('Mean variance per regressor'); colorbar;
            subplot(2,2,4)
            imagesc(inv(cv)./(1./cv)); axis square; title('Covariance inflation factor'); colorbar;
            fprintf(1, '\nVariance estimates: %2.3f\nVariance regressors: %2.3f\nVariance inflation factor: %2.3f (the closer to 1, the better)\n', varE, varX, vif);
        end
        dir_output = fullfile(baseDir, sprintf(glmDir, glm), subj_id);
        tmp = load(fullfile(dir_output,'SPM.mat'));
        delete(fullfile(dir_output,'SPM.mat'));
        SPM = tmp.SPM;
        save(fullfile(dir_output,'SPM.mat'),'SPM','-v7.3');

    case 'GLM:tcontrast'                % 1ST-LEVEL GLM 3: Make t-contrast
        % Condition# 1:anodal 2:sham 3: cathodal
        % 1: overall movement

        vararginoptions(varargin, {'sn','glm','opt'}); %% 
        Names   = {};
        Nums    = [];
        C       = []; % contrast matrix
        % optimize= 1;
        subj_id = char(pinfo.subj_id(pinfo.sn==sn));

        % vararginoptions(varargin(3:end),{'optimize'});
        if glm~=5
            cd(fullfile(baseDir, sprintf(glmDir, glm), subj_id));
            load SPM;
            SPM     = rmfield(SPM,'xCon');  %% temporarily commented out
        else
            load(fullfile(baseDir, sprintf(glmDir, 1), subj_id,'SPM.mat')); % load glm=1
            cd(fullfile(baseDir, sprintf(glmDir, glm), subj_id));            
        end
%         load('SPM_info.mat');  %% load R
%         R = construct_dsgmat(sn,glm);
        %        con = calc_contrasts(R);  %% original contrast
        nRun = length(SPM.nscan);
        nTr = 68;
%         nTr = size(R.transID,2);
        % end
        if glm==0
            contrast_names{1} = 'All_trials';
            optC = [1 0];
        elseif glm==1
            if opt==0
                X = SPM.xX.X(:,1:end-length(SPM.nscan));
                R1 = construct_dsgmat(sprintf('S%02d',sn),glm);
                R2 = construct_dsgmat(sprintf('R%02d',sn),glm);
                R = combine_behavdata(R1,R2);
                for i=1:64 
                    contrast_names{i} = sprintf('Transition_%02d',i);
                    temp = zeros(nRun, nTr);
                    temp(R.transID == i) = 1;
                    optC(i,:) = calc_optWeights2(temp.*R.isValidRep,X);
                end
%                 for i=1:nRun*nTr contrast_names{i} = sprintf('Trial_%02d',i);end
%                 optC = calc_LSSbetasWeights(SPM);   
            else
%     %           tem = load(fullfile('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/glm_1',  sprintf('S%02d', sn), 'SPM.mat'));
%                 contrast_names = {'Motor1','Motor2','Motor2-1','Cue1','Cue2','Cue2-1', 'Both1','Both2','Both2-1'...
%                     'BothLetter1','BothLetter2','BothLetter2-1','BothSpatial1','BothSpatial2','BothSpatial2-1', ...
%                     'CueLetter1','CueLetter2','CueLetter2-1','CueSpatial1','CueSpatial2','CueSpatial2-1',...
%                     'NRepMotor1','NRepMotor2','NRepMotor2-1','NRepCue1','NRepCue2','NRepCue2-1','NRep1','NRep2','NRep2-1',...
%                     'Letter', 'Spatial','Letter-Spatial', 'Letter+Spatial'};
%                 
%                 X = SPM.xX.X(:,1:end-length(SPM.nscan));
%                 
%                 cond_list = {R.isRepMotor, R.isRepCue, R.isRepBoth, R.isRepBothLetter, ...
%                     R.isRepBothSpatial, R.isRepCueLetter, R.isRepCueSpatial, R.isNRepMotor, R.isNRepCue, R.isNRep};
%                     
%                 for i=1:length(cond_list)
%                     optC(3*(i-1)+1:3*i,:) = calc_optWeights(cond_list{i},X);
%                 end
%                 C = zeros(nTr*nRun,2); 
%                 C(find(R.isLetter'==1),1) = 1;
%                 C(:,2) = ones(nTr*nRun,1)-C(:,1);            
%                 optW = inv(C'*X'*X*C)*C'*X'*X;
%                 optC(3*length(cond_list)+1,:) = [optW(1,:)  zeros(1,nRun)];
%                 C = zeros(nTr*nRun,2); 
%                 C(find(R.isSpatial'==1),1) = 1;
%                 C(:,2) = ones(nTr*nRun,1)-C(:,1);            
%                 optW = inv(C'*X'*X*C)*C'*X'*X;
%                 optC(3*length(cond_list)+2,:) = [optW(1,:)  zeros(1,nRun)];
%                 optC(3*length(cond_list)+3,:) = optC(3*length(cond_list)+1,:) - optC(3*length(cond_list)+2,:);
%                 optC(3*length(cond_list)+4,:) = (optC(3*length(cond_list)+1,:) + optC(3*length(cond_list)+2,:))/2;  
%                 contrast_names = {'MotorR','MotorN','MotorR-N','CueR','CueN','CueR-N','BothR','BothN','BothR-N',...
%                               'MotorOnly','CueOnly','MotorOnly-CueOnly',...
%                                'Letter','Spatial','Letter-Spatial',...
%                                 'BothRep-CueOnly','MotorOnly-BothN'
%                                 % 'BothRep-CueOnly-L','BothRep-CueOnly-S','BothRep-CueOnly-L-S',...
%                               % 'MotorOnly-L','MotorOnly-S','CueOnlyRep-Letter','CueOnlyRep-Spatial',...
%                               };


%                 optC(1,:) = calc_optWeights2(R.isRepMotor, X);
%                 optC(2,:) = calc_optWeights2(R.isNRepMotor, X);
%                 optC(3,:) = optC(1,:)-optC(2,:);
%                 optC(4,:) = calc_optWeights2(R.isRepCue, X);
%                 optC(5,:) = calc_optWeights2(R.isNRepCue, X);
%                 optC(6,:) = optC(4,:)-optC(5,:);
%                 optC(7,:) = calc_optWeights2(R.isRepBoth, X);
%                 optC(8,:) = calc_optWeights2(R.isNRep, X);
%                 optC(9,:) = optC(7,:)-optC(8,:);
%                 optC(10,:) = calc_optWeights2(~R.isRepBoth.*R.isRepMotor, X);                
%                 optC(11,:) = calc_optWeights2(~R.isRepBoth.*R.isRepCue, X);
%                 optC(12,:) = optC(10,:)-optC(11,:);
%                 optC(13,:) = calc_optWeights2(R.isLetter, X);
%                 optC(14,:) = calc_optWeights2(R.isSpatial, X);
%                 optC(15,:) = optC(13,:)-optC(14,:);
%                 optC(16,:) = optC(7,:)-optC(11,:);
%                 optC(17,:) = optC(10,:)-optC(8,:);


                % optC(12,:) = calc_optWeights2(R.isRepBoth.*R.isLetter,X)-calc_optWeights2(~R.isRepBoth.*R.isRepMotor.*R.isLetter, X);
                % optC(13,:) = calc_optWeights2(R.isRepBoth.*R.isSpatial,X)-calc_optWeights2(~R.isRepBoth.*R.isRepMotor.*R.isSpatial, X);
                % optC(14,:) = calc_optWeights2(R.isRepBoth.*R.isLetter,X)-calc_optWeights2(~R.isRepBoth.*R.isRepCue.*R.isLetter, X);
                % optC(15,:) = calc_optWeights2(R.isRepBoth.*R.isSpatial,X)-calc_optWeights2(~R.isRepBoth.*R.isRepCue.*R.isSpatial, X);
            end
        elseif glm==2 
              % cond_name = {'MotorOnly-L','MotorOnly-S','CueOnly-L','CueOnly-S',...
              %              'BothRep-L','BothRep-S','NonRep-L','NonRep-S', 'Non-Interest'};
%             contrast_names = {'MotorR','MotorN','MotorR-MotorN','CueR','CueN','CueR-CueN','BothR','BothN','BothR-BothN','MotorR-BothR', 'CueR-BothR'};
            %contrast_names = {'MotorR-MotorN', 'CueR-CueN','BothR-BothN'}; 

              contrast_names = {'BothRep-L','CueRep-L','MotorRep-L','NRep-L','BothRep-S','CueRep-S','MotorRep-S','NRep-S',...
                  'wRS-L','wRS-S','wRS-L-S','acRS-L','acRS-S','acRS-L-S','Letter','Spatial','Letter-Spatial'};
              optC(1,:)  = [repmat([0 0 0 0 1 0 0 0 0],1,nRun) zeros(1,nRun)]/nRun;
              optC(2,:)  = [repmat([0 0 1 0 0 0 0 0 0],1,nRun) zeros(1,nRun)]/nRun;
              optC(3,:)  = [repmat([1 0 0 0 0 0 0 0 0],1,nRun) zeros(1,nRun)]/nRun;
              optC(4,:)  = [repmat([0 0 0 0 0 0 1 0 0],1,nRun) zeros(1,nRun)]/nRun;
              optC(5,:)  = [repmat([0 0 0 0 0 1 0 0 0],1,nRun) zeros(1,nRun)]/nRun;
              optC(6,:)  = [repmat([0 0 0 1 0 0 0 0 0],1,nRun) zeros(1,nRun)]/nRun;
              optC(7,:)  = [repmat([0 1 0 0 0 0 0 0 0],1,nRun) zeros(1,nRun)]/nRun;
              optC(8,:)  = [repmat([0 0 0 0 0 0 0 1 0],1,nRun) zeros(1,nRun)]/nRun;
              optC(9,:)  = [repmat([0 0 -1 0 1 0 0 0 0],1,nRun) zeros(1,nRun)]/nRun;
              optC(10,:) = [repmat([0 0 0 -1 0 1 0 0 0],1,nRun) zeros(1,nRun)]/nRun;
              optC(11,:) = [repmat([0 0 -1 -1 1 1 0 0 0],1,nRun) zeros(1,nRun)]/(2*nRun);
              optC(12,:)  = [repmat([1 0 0 0 0 0 -1 0 0],1,nRun) zeros(1,nRun)]/nRun;
              optC(13,:)  = [repmat([0 1 0 0 0 0 0 -1 0],1,nRun) zeros(1,nRun)]/nRun;
              optC(14,:)  = [repmat([1 1 0 0 0 0 -1 -1 0],1,nRun) zeros(1,nRun)]/(2*nRun);
              optC(15,:)  = (optC(1,:)+optC(2,:)+optC(3,:)+optC(4,:))/4;
              optC(16,:)  = (optC(5,:)+optC(6,:)+optC(7,:)+optC(8,:))/4;
              optC(17,:)  = optC(15,:)-optC(16,:);
              
%             contrast_names = {'MotorR','MotorN','MotorR-N',...
%                                 'CueR','CueN','CueR-N','BothR','BothN','BothR-N',...
%                                 'MotorOnly','CueOnly','MotorOnly-CueOnly',...
%                                 'MotorR-L','MotorR-S','MotorR-L-S',...
%                                 'CueR-L','CueR-S','CueR-L-S',...
%                                 'Letter','Spatial','Letter-Spatial',...
%                                 'BothRep-L','BothRep-S','BothRep-L-S',...
%                                 'CueOnly-L','CueOnly-S','CueOnly-L-S',...
%                                 'BothRep-CueOnly-L','BothRep-CueOnly-S','BothRep-CueOnly-L-S',...
%                                 'BothRep-CueOnly','MotorOnly-BothN'};
%             optC(1,:) = [repmat([1 1 0 0 1 1 0 0 0],1, nRun) zeros(1,nRun)]/(4*nRun);
%             optC(2,:) = [repmat([0 0 1 1 0 0 1 1 0],1, nRun) zeros(1,nRun)]/(4*nRun);
%             optC(3,:) = optC(1,:)-optC(2,:);
%             optC(4,:) = [repmat([0 0 1 1 1 1 0 0 0],1, nRun) zeros(1,nRun)]/(4*nRun);
%             optC(5,:) = [repmat([1 1 0 0 0 0 1 1 0],1, nRun) zeros(1,nRun)]/(4*nRun);
%             optC(6,:) = optC(4,:)-optC(5,:);
%             optC(7,:) = [repmat([0 0 0 0 1 1 0 0 0],1, nRun) zeros(1,nRun)]/(2*nRun);
%             optC(8,:) = [repmat([0 0 0 0 0 0 1 1 0],1, nRun) zeros(1,nRun)]/(2*nRun);
%             optC(9,:) = optC(7,:)-optC(8,:);
%             optC(10,:) = [repmat([1 1 0 0 0 0 0 0 0],1, nRun) zeros(1,nRun)]/(2*nRun);
%             optC(11,:) = [repmat([0 0 1 1 0 0 0 0 0],1, nRun) zeros(1,nRun)]/(2*nRun);
%             optC(12,:) = optC(10,:)-optC(11,:);    
% 
%             optC(13,:) = [repmat([1 0 0 0 1 0 0 0 0],1, nRun) zeros(1,nRun)]/(2*nRun);
%             optC(14,:) = [repmat([0 1 0 0 0 1 0 0 0],1, nRun) zeros(1,nRun)]/(2*nRun);
%             optC(15,:) = optC(13,:)-optC(14,:);
%             optC(16,:) = [repmat([0 0 1 0 1 0 0 0 0],1, nRun) zeros(1,nRun)]/(2*nRun);
%             optC(17,:) = [repmat([0 0 0 1 0 1 0 0 0],1, nRun) zeros(1,nRun)]/(2*nRun);
%             optC(18,:) = optC(16,:)-optC(17,:);
%             optC(19,:) = [repmat([1 0 1 0 1 0 1 0 0],1, nRun) zeros(1,nRun)]/(4*nRun);
%             optC(20,:) = [repmat([0 1 0 1 0 1 0 1 0],1, nRun) zeros(1,nRun)]/(4*nRun);
%             optC(21,:) = optC(19,:)-optC(20,:);
%             optC(22,:) = [repmat([0 0 0 0 1 0 0 0 0],1, nRun) zeros(1,nRun)]/(nRun);
%             optC(23,:) = [repmat([0 0 0 0 0 1 0 0 0],1, nRun) zeros(1,nRun)]/(nRun);
%             optC(24,:) = optC(22,:)-optC(23,:);
% 
%             optC(25,:) = [repmat([0 0 1 0 0 0 0 0 0],1, nRun) zeros(1,nRun)]/(nRun);
%             optC(26,:) = [repmat([0 0 0 1 0 0 0 0 0],1, nRun) zeros(1,nRun)]/(nRun);
%             optC(27,:) = optC(25,:)-optC(26,:);
% 
%             optC(28,:) = [repmat([0 0 -1 0 1 0 0 0 0],1, nRun) zeros(1,nRun)]/(nRun);
%             optC(29,:) = [repmat([0 0 0 -1 0 1 0 0 0],1, nRun) zeros(1,nRun)]/(nRun);
%             optC(30,:) = optC(28,:)-optC(29,:);
%             optC(31,:) = [repmat([0 0 -1 -1 1 1 0 0 0],1, nRun) zeros(1,nRun)]/(2*nRun);
%             optC(32,:) = optC(10,:)-optC(8,:);

        elseif glm==3  %% encoding trial-state, 9th state is for non-interest
            load('SPM_info.mat');
            n_cond = length(unique(R.cond));
            for c=1:n_cond 
                contrast_names{c} = sprintf('Trial-State %d',c);
                temp = zeros(1,n_cond);temp(c) = 1;
                optC(c,:) = [repmat(temp, 1, nRun) zeros(1,nRun)];
                optC = optC/nRun;
            end
        elseif glm==4
            contrast_names = {'Letter', 'Spatial','Letter-Spatial', 'Letter+Spatial'};
            if length(unique(R.cond))==3
                optC(1,:) = [repmat([1 0 0],1, nRun) zeros(1,nRun)];
                optC(2,:) = [repmat([0 1 0],1, nRun) zeros(1,nRun)];
                optC(3,:) = [repmat([1 -1 0],1, nRun) zeros(1,nRun)];
                optC(4,:) = [repmat([0.5 0.5 0], 1, nRun) zeros(1,nRun)];
                optC = optC/nRun;
            else    %% two conditions
                optC(1,:) = [repmat([1 0],1, nRun) zeros(1,nRun)];
                optC(2,:) = [repmat([0 1],1, nRun) zeros(1,nRun)];
                optC(3,:) = [repmat([1 -1],1, nRun) zeros(1,nRun)];
                optC(4,:) = [repmat([0.5 0.5], 1, nRun) zeros(1,nRun)];
                optC = optC/nRun;
            end
        elseif glm==5
%                 contrast_names = {'BothRep-L','CueRep-L','MotorRep-L','NRep-L','BothRep-S','CueRep-S','MotorRep-S','NRep-S'};
               contrast_names = {'BothRep-L','CueRep-L','MotorRep-L','NRep-L','BothRep-S','CueRep-S','MotorRep-S','NRep-S',...
                  'wRS-L','wRS-S','wRS-L-S','acRS-L','acRS-S','acRS-L-S','Letter','Spatial','Letter-Spatial'};
                R1 = construct_dsgmat(sprintf('S%02d',sn),2); % Use glm=2 for the constrast
                R2 = construct_dsgmat(sprintf('R%02d',sn),2);  % Use glm=2 for the contrast
                R = combine_behavdata(R1,R2);
                X = SPM.xX.X(:,1:end-length(SPM.nscan));
                for c=1:8 temp = zeros(nRun, nTr); temp(R.cond==c)=1; optC(c,:) = calc_optWeights2(temp,X);end
                optC(9,:) = optC(1,:)-optC(2,:);
                optC(10,:) = optC(5,:)-optC(6,:);
                optC(11,:) = optC(9,:)-optC(10,:);
                optC(12,:) = optC(3,:)-optC(4,:);
                optC(13,:) = optC(7,:)-optC(8,:);
                optC(14,:) = optC(12,:)-optC(13,:);
                optC(15,:)  = (optC(1,:)+optC(2,:)+optC(3,:)+optC(4,:))/4;
                optC(16,:)  = (optC(5,:)+optC(6,:)+optC(7,:)+optC(8,:))/4;
                optC(17,:)  = optC(15,:)-optC(16,:);        
        end
        
        for i=1:length(contrast_names)
            SPM.xCon(i) = spm_FcUtil('Set', contrast_names{i}, 'T', 'c', optC(i,:)', SPM.xX.xKXs);
        end

        SPM = spm_contrasts(SPM,1:length(SPM.xCon));
        save('SPM.mat', 'SPM','-v7.3');
        SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
        save(fullfile(baseDir, sprintf(glmDir, glm), subj_id, 'SPM_light.mat'), 'SPM')
       % rename contrast images and spmT images
        conName = {'con','spmT'};
        for i = 1:length(SPM.xCon)
            for n=1
%             for n = 1:numel(conName)
                oldName = fullfile(baseDir, sprintf(glmDir, glm), subj_id, sprintf('%s_%2.4d.nii',conName{n},i));
                newName = fullfile(baseDir, sprintf(glmDir, glm), subj_id, sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                movefile(oldName, newName);
            end % conditions (n, conName: con and spmT)
        end % i (contrasts)
        
    
%     case 'GLM:beta_estimate'                % Estimation of beta-value from glm=1 results encoding individual trials
%         glm = 1;  %% Using trial-wise beta estimats from glm=1
%         cd(fullfile(baseDir, sprintf(glmDir, glm), sprintf('S%02d', sn)));
%         load('SPM.mat');
%         optC = calc_LSSbetasWeights(SPM);

    case 'GLM:psc' % calculate percent signal change for selected contrasts
        vararginoptions(varargin, {'sn', 'glm'});
        
        % Loop over the subjects
        subj_id = char(pinfo.subj_id(pinfo.sn==sn));
        % Go to subject's directory and load SPM info
        cd(fullfile(baseDir,sprintf(glmDir, glm), subj_id));
        if glm~=5
            load SPM;
        else
            load(fullfile(baseDir,sprintf(glmDir, 1), subj_id,'SPM.mat'));
        end
%        load('SPM_info.mat');
        %load('SPM_info.mat');

        X = (SPM.xX.X(:,SPM.xX.iC));      % Design matrix - raw
        h = median(max(X));               % Height of response;
        numB = length(SPM.xX.iB);         % Partitions - runs
        P = cell(1,numB+1);
        r = 0;
        for p=SPM.xX.iB
            r=r+1;
            P{r}=sprintf('beta_%4.4d.nii',p);       % get the intercepts and use them to calculate the baseline (mean images)
        end

        for con=1:length(SPM.xCon)    % all contrasts
            P{numB+1}=sprintf('con_%s.nii',SPM.xCon(con).name);           % replace with T.con_name{con} in case
            outname=sprintf('psc_%s.nii',SPM.xCon(con).name);
            formula = sprintf(['100.*%f.*' pscFormula(numB)],h);
            spm_imcalc_ui(P,outname,formula,{0,[],spm_type(16),[]});       % Calculate percent signal change
            fprintf('Contrast: %s\n',SPM.xCon(con).name);
        end

        fprintf('Subject %s - Done\n', subj_id);

    
    case 'MNI:norm_write'
        vararginoptions(varargin, {'sn','glm'}); %% 
        groupDir = fullfile(baseDir,sprintf(glmDir,glm),'group');
        % for s = sn
            
        disp(sn)
    
        P = fullfile(baseDir, anatomicalDir,sprintf('S%02d',sn),sprintf('y_S%02d_anatomical.nii',sn));
        [Def,mat] = spmdefs_get_def(P);

        % apply deformation to contrasts
        intrp = 4;
        fnames = {};
        ofnames = {};
        load(fullfile(baseDir,sprintf(glmDir,glm),sprintf('S%02d',sn),'SPM_light.mat'));
        for c=1:length(SPM.xCon)
            fnames{c}=fullfile(baseDir, sprintf(glmDir,glm),sprintf('S%02d',sn),sprintf('con_%s.nii',SPM.xCon(c).name));
            ofnames{c} = fullfile(groupDir,sprintf('con_%s_S%02d.nii',SPM.xCon(c).name,sn));
            spmdefs_apply_def(Def,mat,fnames{c},intrp,ofnames{c})
        end


    
        % end

    case 'SURF:freesurfer'              % SURFACE PREPROCESS 1: Run recon-all on freesurfer
        vararginoptions(varargin, {'sn'}); %% 
        subj_id = char(pinfo.subj_id(pinfo.sn==sn));
        S_id = strrep(subj_id,'R','S');
        freesurf_dir = fullfile(baseDir, freesurferDir);
        if (~exist(freesurf_dir,'dir'))
            mkdir(freesurf_dir);
        end
        %mysetEnv(freesurferDir);
        freesurfer_reconall(freesurf_dir,S_id,fullfile(baseDir,anatomicalDir,S_id,[S_id '_anatomical.nii']));  
    case 'WB:surf_resample' % reslice indiv surfaces into fs_lr standard mesh
        % sn = subj_vec; cwd = pwd;
        surf  = '32k'; % 164k or 32k vertices
        vararginoptions(varargin, {'sn'});
        subj_id = char(pinfo.subj_id(pinfo.sn==sn));
        S_id = strrep(subj_id,'R','S');
        % for s = sn
        fprintf('reslicing %s\n', subj_id);
        freesurferDir = fullfile(baseDir, freesurferDir);
        wbDir = fullfile(baseDir, wbDir);
        surf_resliceFS2WB(S_id, freesurferDir, wbDir, 'resolution', surf);
        % end

    case 'WB:vol2surf_indiv' % map indiv vol contrasts (.nii) onto surface (.gifti)
            % sss_imana('GLM:psc','sn',s) ->
            % glm = 3;
            surf = '32k';
            hemis = [1 2];
            hem = {'L','R'};
            % map = 'psc'; % 't'; 'con'; 'dist'; 'psc';
            % vararginoptions(varargin,{'sn', 'glm', 'hemis', 'map', 'surf'});
            vararginoptions(varargin,{'sn', 'glm' 'map'});
            subj_id = char(pinfo.subj_id(pinfo.sn==sn));
            S_id = strrep(subj_id,'R','S');
            for h = hemis
                surfDir = fullfile(baseDir, wbDir, S_id);
                white = fullfile(surfDir, sprintf('%s.%s.white.%s.surf.gii', S_id, hem{h}, surf));
                pial = fullfile(surfDir, sprintf('%s.%s.pial.%s.surf.gii', S_id, hem{h}, surf));
                C1 = gifti(white);
                C2 = gifti(pial);
                if glm~=5
                    subjGLM = fullfile(baseDir,sprintf(glmDir, glm), subj_id);
                else
                    subjGLM =  fullfile(baseDir,sprintf(glmDir, 1), sprintf('S%02d', sn));
                end
                load(fullfile(subjGLM, 'SPM.mat')); 
                switch map
                    case 't' % t-values maps (univariate GLM)
                        fnames      = cell(1,numel(SPM.xCon));
                        con_name    = cell(1,numel(SPM.xCon));
                        for j=1:numel(fnames)
                            fnames{j}   = fullfile(subjGLM, sprintf('spmT_%s.nii', SPM.xCon(j).name));
                            con_name{j} = SPM.xCon(j).name;
                        end
                    case 'con' % contrast beta maps (univariate GLM)
                        fnames      = cell(1,numel(SPM.xCon));
                        con_name    = cell(1,numel(SPM.xCon));
                        for j=1:numel(fnames)
                            fnames{j}   = fullfile(subjGLM, sprintf('con_%s.nii', SPM.xCon(j).name));
                            con_name{j} = SPM.xCon(j).name;
                        end    
                    case 'psc' % percent signal change maps (univariate GLM)
                        fnames      = cell(1,numel(SPM.xCon));
                        con_name    = cell(1,numel(SPM.xCon));
                        for j=1:numel(fnames)
                            fnames{j}   = fullfile(subjGLM, sprintf('psc_%s.nii', SPM.xCon(j).name));
                            con_name{j} = SPM.xCon(j).name;
                        end
                    case 'rdm'
                        fname{1} = fullfile(baseDir, 'patterns', sprintf('S%02d', sn),'rdm.nii');
                        con_name{1} = 'rdm';
                    case 'res' % residual
                        fnames{1} = fullfile(subjGLM, 'ResMS.nii');
                        con_name{1} = 'ResMS';
                end
                outfile = fullfile(surfDir, sprintf('%s.%s.glm%d.%s.func.gii', subj_id, hem{h}, glm, map));
                G = surf_vol2surf(C1.vertices, C2.vertices, fnames, 'column_names', con_name, 'anatomicalStruct', hem{h}, 'exclude_thres', 0.75, 'faces', C1.faces, 'ignore_zeros', 0);
                save(G, outfile);
                fprintf('mapped %s %s glm%d \n', subj_id, hem{h}, glm);
            end

    case 'WB:vol2surf_group' % map group contrasts on surface (.gifti)
        %sn = subj_vec;
        %glm = 1;
        hemis = [1 2];
        map = 'psc'; % 't'; 'con'; 'dist'; 'psc';
        surf = '32k';
        vararginoptions(varargin,{'sn', 'glm','map','prefix'});
        %vararginoptions(varargin,{'sn'});
        groupDir = fullfile(baseDir, wbDir, sprintf('glm%d',glm), map, sprintf('group%s', surf));
%         dircheck(groupDir);
        if ~exist(groupDir, 'dir')
            mkdir(groupDir);
        end
        subjGLM = fullfile(baseDir,sprintf(glmDir, glm), 'S01');
        load(fullfile(subjGLM, 'SPM.mat'));
        for h = hemis
            fprintf(1, '%s ... ', hemi{h});
            inputFiles = {};
            columnName = {};
            for s = sn
                subj_id = char(pinfo.subj_id(pinfo.sn==s));
                S_id = strrep(subj_id,'R','S');
                inputFiles{end+1} = fullfile(baseDir, wbDir, S_id, sprintf('%s.%s.glm%d.%s.func.gii', subj_id, hem{h}, glm, map));
                columnName{end+1} = subj_id;
            end
            groupfiles = cell(1);
            switch map
                case 't' % t-values maps (univariate GLM)
                    con = SPM.xCon;
                    nc  = numel(con);
                    for ic = 1:nc
                        groupfiles{ic}  = fullfile(groupDir, sprintf('group.%s.%s.%s.glm%d.%s.func.gii', prefix, map, hem{h}, glm, con(ic).name));
                    end
                case 'con' % contrast beta maps (univariate GLM)
                    con = SPM.xCon;
                    for ic = 1:numel(con)
                        groupfiles{ic}  = fullfile(groupDir, sprintf('group.%s.%s.%s.glm%d.%s.func.gii', prefix, map, hem{h}, glm, con(ic).name));
                    end
                case 'psc' % percent signal change maps (univariate GLM)
                    con = SPM.xCon;
                    for ic = 1:numel(con)
                        groupfiles{ic}  = fullfile(groupDir, sprintf('group.%s.%s.%s.glm%d.%s.func.gii', prefix, map, hem{h}, glm, con(ic).name));
                    end
            end
            surf_groupGiftis(inputFiles, 'outcolnames', columnName, 'outfilenames', groupfiles);
            fprintf(1, 'Done.\n');
        end

    case 'WB:searchlight'
        S = rsa.rea;
        % EXAMPLE 1:
%   % Define a surface-based searchlight for the left hemisphere only
%   S = rsa.readSurf({'lh.white.surf.gii'},{'lh.pial.surf.gii'});
%   M = rsa.readMask('mask.img');
%   L = rsa.defineSearchlight_surface(S,M);
%
% EXAMPLE 2:
%   % Define a surface-based searchlight for both hemispheres, over a
%   80 voxel searchlight with a radius of 20mm
%   white   = {'lh.white.surf.gii', 'rh.white.surf.gii'};
%   pial    = {'lh.pial.surf.gii' , 'rh.pial.surf.gii'};
%   S       = rsa.readSurf(white,pial);
%   M       = rsa.readMask('mask.img');
%   L       = rsa.defineSearchlight(S,'mask',M,'sphere',[20 80]);

    case 'WB:vol2surf_stats' % perform group stats on surface (.gifti)
        % sn = subj_vec;
%         glm   = 1;
        hemis = [1 2];
        map   = 'psc'; % 't'; 'con'; 'psc';
        sm    = 3;     % smoothing kernel in mm (optional)
        surf  = '32k';  % 164k or 32k vertices
        vert  = str2double(regexp(surf, '\d+', 'match'));
        % vararginoptions(varargin,{'sn', 'glm', 'hemis', 'map', 'sm', 'surf'});
        vararginoptions(varargin,{'glm','map','prefix'});
        groupDir = fullfile(baseDir, wbDir, sprintf('glm%d',glm), map, sprintf('group%s',surf));
        subjGLM = fullfile(baseDir,sprintf(glmDir, glm), 'S01');
        load(fullfile(subjGLM, 'SPM.mat'));
        % Loop over the metric files and calculate the cSPM of each
        for h = hemis
            fprintf(1, '%s ... \n', hemi{h});
            groupfiles = cell(1);
            switch map
                case 't' % t-values maps (univariate GLM)
                    con = SPM.xCon;
                    nc  = numel(con);
                    con_name = cell(1);
                    for ic = 1:nc
                        con_name{ic} = con(ic).name;
                    end
                case 'con' % contrast beta maps (univariate GLM)
                    con = SPM.xCon;
                    nc  = numel(con);
                    con_name = cell(1);
                    for ic = 1:nc
                        con_name{ic} = con(ic).name;
                    end                    
                case 'psc' % percent signal change maps (univariate GLM)
                    con = SPM.xCon;
                    nc  = numel(con);
                    con_name = cell(1);
                    for ic = 1:nc
                        con_name{ic} = con(ic).name;
                    end
            end
            % Perform stats
            for ic = 1:nc
                groupfiles{ic} = fullfile(groupDir,sprintf('group.%s.%s.%s.glm%d.%s.func.gii', prefix, map, hem{h}, glm, con_name{ic}));
%                 metric = gifti(groupfiles{ic});
%                 data = double(metric.cdata);
%                 
%                 % Do the uw noise normalization
%                 subj_count=0;
%                 for s = sn
%                     subj_count=subj_count+1;
%                     surfDir = fullfile(wbDir, sprintf('s%02d', s));
%                     res = gifti(fullfile(surfDir, sprintf('%s.%s.glm%d.res.func.gii', sprintf('s%02d', s), hem{h}, glm)));
%                     res = double(res.cdata);
%                     data(res>0,subj_count) = data(res>0,subj_count)./sqrt(res(res>0)); 
%                 end
%                 
%                 % save the new group map for this particular contrast
%                 G = surf_makeFuncGifti(data,'anatomicalStruct',hname{h},'columnNames',surf_getGiftiColumnNames(metric));
%                 newname = fullfile(groupDir, sprintf('uwgroup.%s.%s.glm%d.%s.func.gii', map, hem{h}, glm, con_name{ic}));
%                 save(G,newname);
                
                
                % Perform smoothing (optional)
                if sm>0
                    surface = fullfile(atlasDir, sprintf('fS_LR_%s/fs_LR.%s.%s.flat.surf.gii', strrep(surf,'k',''), surf, hem{h}));
                    groupfiles{ic}  = surf_smooth(groupfiles{ic}, 'prefix', sprintf('smooth.%s.',prefix), 'surf', surface, 'kernel', vert);
                    %groupfiles{ic}  = surf_smooth(newname, 'surf', surface, 'kernel', sm);
                end
                metric          = gifti(groupfiles{ic});
                
                
                
                cSPM            = surf_getcSPM('onesample_t', 'data', metric.cdata); %, 'maskthreshold',0.5); % set maskthreshold to 0.5 = calculate stats at location if 50% of subjects have data at this point
                C.data(:, ic)   = cSPM.con.con; % mean
                C.c_name{ic}    = ['mean_' con_name{ic}];
                C.data(:,ic+nc) = cSPM.con.Z;   % t (confusing)
                C.c_name{ic+nc} = ['t_' con_name{ic}];
            end
            % Save output
            O = surf_makeFuncGifti(C.data, 'columnNames', C.c_name, 'anatomicalStruct', hname{h});
            %summaryfile = fullfile(groupDir, sprintf('uwsummary.%s.glm%d.%s.sm%d.func.gii', hem{h}, glm, map, sm));
            summaryfile = fullfile(groupDir, sprintf('summary.%s.%s.glm%d.%s.sm%d.func.gii', prefix, hem{h}, glm, map, sm));
            %uwsummary.L.glm3.psc.sm8.func.gii
            save(O, summaryfile);
            fprintf(1, 'Done.\n');
        end
    case 'HRF:fit' % finding optimal parameters for hrf
        % we have a instruction and a movement trial type and we want to
        % find a set of parameter that works for both. The important thing
        % is that the duration of underlying neural event is the only thing
        % that is different between these two trial types.

        sn = subj_vec;
        glm = 0;
        regN = [1:7]; % SMA, PMv, PMd, M1, S1, aSPL, pSPL
        duration = 1;
        onsetshift = 0;
        pre = 5;
        post = 25;
        fig = 1;
%         vararginoptions(varargin, {'sn', 'glm', 'regN', 'duration', 'onsetshift', 'pre', 'post', 'fig'});
        vararginoptions(varargin, {'sn','regN','fitMethod'});
%         cwd = pwd;
        subj_id = pinfo.subj_id{pinfo.sn==sn};
        
        % loop over the subject - TODO: group fitting
        for s = sn
            % initialize
            %T = [];
            %FIT = [];
            
            fprintf('subject %.2d\n', s);
            
            cd(fullfile(baseDir,sprintf('glm_%d',glm),subj_id)); % cd to subject's GLM dir
            load SPM;
            
            % % updating the filenames - since I ran the code on server
            % if ~strcmp(fullfile(baseDir,sprintf(glmDir, glm), sprintf('S%02d', s)), SPM.swd) % need to rename SPM
            %     SPM = spmj_move_rawdata(SPM, fullfile(imagingDir, sprintf('S%02d', s)));
            % endS
            
            load(fullfile(baseDir, roiDir, sprintf('%s_SSS_regions.mat', sprintf('S%02d', s)))); % Load R
            %load(fullfile(roiDir, sprintf('%s_Wang_regions.mat', sprintf('s%02d', s))));
            R = R(regN); % only use the specified regions
            
            Data = region_getdata(SPM.xY.VY, R);
            
            reg=[];
            data=[];
            for i = 1:length(Data)
                reg = [reg ones(1,size(Data{i},2))*regN(i)];
                data = [data Data{i}];
            end
            clear Data
            
            Y = spm_filter(SPM.xX.K, SPM.xX.W*data); % filter out low-frequence trends in Y
            Yres = spm_sp('r', SPM.xX.xKXs, Y);
            err_before = sum(sum(Yres.^2))/numel(Yres);
            
%             for r = 1:length(SPM.nscan)
%                 for u=1:length(SPM.Sess(r).U)
%                     SPM.Sess(r).U(u).dur = ones(size(SPM.Sess(r).U(u).dur))*duration;
%                     SPM.Sess(r).U(u).ons = SPM.Sess(r).U(u).ons+onsetshift;
%                 end
%                 SPM.Sess(r).U=spm_get_ons(SPM,r);
%             end
            
            % Fit a common hrf for specified regions
%             [SPMf, Yhat, Yres] = fit_hrf(SPM, data);  % 'fit',[1,2]'  
            [SPMf, Yhat, Yres, p_opt] = spmj_fit_hrfparams(SPM, data,fitMethod);  % 'fit',[1,2]'  
            % Check Error after
            err_after = sum(sum(Yres.^2))/numel(Yres);
            F.bf_before = SPM.xBF.bf;
            F.bf_after = SPMf.xBF.bf;
            F.params_before = SPM.xBF.params;
            F.params_after = p_opt;
           
            F.err_before = err_before;
            F.err_after = err_after;
            [err_before err_after p_opt]
            g = 0;
            figure;
            set(gcf,'color','w');
            for r=1:8
                subplot(2,4,r);
                sss_imana('HRF:ROI_hrf_get','sn',s,'glm',g,'post',20);
                sss_imana('HRF:ROI_hrf_plot','sn',s,'roi',r,'glm',g,'post',20);
                sss_imana('HRF:ROI_hrf_get','sn',s,'glm',g,'post',20,'bf',SPMf.xBF.bf);
                load(fullfile(baseDir,roiDir,sprintf('S%0.2d_glm%d_hrf.mat',s,glm))); % load T 
                T = getrow(T,T.region==r);
                pre = 10; post = 20;
                % Select a specific subset of things to plot 
                if glm==2
                    T.type(find(T.type==6))=5; %% BothRep
                    T.type(find(T.type==4))=3; %% CueOnlyRep
                    subset  = find(T.type==3 | T.type==5);
                elseif glm==0
                    subset = find(T.type==1);
                end
                hold on;
                traceplot([-pre:post],T.y_hat,'linestyle','--',...
                'split',[T.type],'subset',subset,...
                'linewidth',3,'linecolor','r'); % ,
                drawline([-8 6 16],'dir','vert','linestyle','--');
                drawline([0],'dir','horz','linestyle','--');
%                 hold off;
                xlabel('TR');
                ylabel('activation');
%                 title(sprintf('ROI: %s',regname{roi}));
                drawline(0);
            end
            
%             figure;
%             for r=1:8
%                 subplot(2,4,r);
%                 sss_imana('HRF:ROI_hrf_plot','sn',s,'roi',r,'glm',g,'post',20);
%             end
%             
            
            %FIT = addstruct(FIT,F);
            
            % save
%             save(fullfile(baseDir, roiDir, sprintf('%s_hrf_ROI_timeseries_glm%d.mat', sprintf('S%02d', s), glm)), '-struct', 'T');
            save(fullfile(baseDir, roiDir, sprintf('%s_hrf_fit_glm%d.mat', sprintf('S%02d', s), glm)), '-struct', 'F');            
        end

    case 'HRF:ROI_hrf_get'                   % Extract raw and estimated time series from ROIs
                
        sn = [];
        ROI = 'all';
        pre=10;
%        post=30;
        atlas = 'SSS';
        glm = 3; % change this glm=1 for S01-06
        bf = [];
%         vararginoptions(varargin,{'ROI','pre','post', 'glm', 'sn', 'atlas'});
        vararginoptions(varargin,{'sn','glm','post','bf'});
        glmDir = fullfile(baseDir,sprintf(glmDir,glm));
        T=[];
    
        subj_id = pinfo.subj_id{pinfo.sn==sn};
        fprintf('%s\n',subj_id);
    
        % load SPM.mat
        cd(fullfile(glmDir,subj_id));
        SPM = load('SPM.mat'); SPM=SPM.SPM;
        if ~isempty(bf)
            SPM.xBF.bf = bf;
            SPM = fMRI_design_changeBF(SPM);
        end
        
        % load ROI definition (R)
        % R = load(fullfile(baseDir, roiDir,[subj_id '_' atlas '_regions.mat']));
        R = load(fullfile(baseDir, roiDir,sprintf('%s_Task_regions.glm_3.mat',subj_id)));
        R=R.R;
        
        % extract time series data
        [y_raw, y_adj, y_hat, y_res,B] = region_getts(SPM,R);
        
        D = spmj_get_ons_struct(SPM);
        
        for r=1:size(y_raw,2)
            for i=1:size(D.block,1)
                D.y_adj(i,:)=cut(y_adj(:,r),pre,round(D.ons(i))-1,post,'padding','nan')';
                D.y_hat(i,:)=cut(y_hat(:,r),pre,round(D.ons(i))-1,post,'padding','nan')';
                D.y_res(i,:)=cut(y_res(:,r),pre,round(D.ons(i))-1,post,'padding','nan')';
                D.y_raw(i,:)=cut(y_raw(:,r),pre,round(D.ons(i))-1,post,'padding','nan')';
%                 D.y_adj(i,:)=cut(y_adj(:,r),pre,round(D.ons(i)),post,'padding','nan')';
%                 D.y_hat(i,:)=cut(y_hat(:,r),pre,round(D.ons(i)),post,'padding','nan')';
%                 D.y_res(i,:)=cut(y_res(:,r),pre,round(D.ons(i)),post,'padding','nan')';
%                 D.y_raw(i,:)=cut(y_raw(:,r),pre,round(D.ons(i)),post,'padding','nan')';
            end
            
            % Add the event and region information to tje structure. 
            len = size(D.event,1);                
            D.SN        = ones(len,1)*sn;
            D.region    = ones(len,1)*r;
            D.name      = repmat({R{r}.name},len,1);
%            D.hem       = repmat({R{r}.hem},len,1);
            D.type      = D.event; 
            T           = addstruct(T,D);
        end
        
        save(fullfile(baseDir,roiDir, sprintf('%s_glm%d_hrf.mat',subj_id,glm)),'T'); 
        varargout{1} = T;

   
    case 'HRF:ROI_hrf_plot'                 % Plot extracted time series
        % s = varargin{1};
        % roi = varargin{2}; 
        vararginoptions(varargin,{'sn','roi','glm','post','subset'});
        subj_id = pinfo.subj_id{pinfo.sn==sn};
        regname = {'SMA','PMv','PMd','M1','S1','SPLa','SPLp','DSVC','MT+','VSVC','EAC'};
        load(fullfile(baseDir,roiDir,sprintf('%s_glm%d_hrf.mat',subj_id,glm))); % load T 
        T = getrow(T,T.region==roi);
        pre = 10;
%        post = 30;
        % Select a specific subset of things to plot 
%         cond_name = {'MotorOnly-L','MotorOnly-S','CueOnly-L','CueOnly-S',...
%                             'BothRep-L','BothRep-S','NonRep-L','NonRep-S','Non-Interest'};
        if glm==2
%             T.type(find(T.type==6))=5; %% BothRep
%             T.type(find(T.type==4))=3; %% CueOnlyRep
            subset  = find(T.type==3 | T.type==5);
        elseif glm==0
            subset = find(T.type==1);
        end
       % figure;        
%         traceplot([-pre:post],T.y_adj,'errorfcn','stderr',...
%             'split',[T.type],'subset',subset,...
%             'leg',regname{roi},'leglocation','bestoutside'); % ,
        traceplot([-pre:post],T.y_adj,'errorfcn','stderr',...
            'split',[T.type],'subset',subset); % ,
        hold on;
        traceplot([-pre:post],T.y_hat,'linestyle','--',...
            'split',[T.type],'subset',subset,...
            'linewidth',3); % ,
        drawline([-8 6 16],'dir','vert','linestyle','--');
        drawline([0],'dir','horz','linestyle','--');
        hold off;
        xlabel('TR');
        ylabel('activation');
        title(sprintf('ROI: %s',regname{roi}));
        drawline(0);


    case 'ROI:findall' % use this... %rdm, glm3
        h=1;
        surf='32';
        
        % 
        ROI{1,1}={'L_SCEF_ROI','L_6ma_ROI','L_6mp_ROI'}; % SMA
        ROI{1,2}={'L_6v_ROI','L_PEF_ROI','L_55b_ROI'}; % PMv
        ROI{1,3}={'L_6a_ROI','L_6d_ROI','L_FEF_ROI'}; % PMd
        % I took M1 and S1 from the one that we have already
        ROI{1,6}={'L_AIP_ROI','L_7PC_ROI','L_LIPv_ROI','L_LIPd_ROI'}; % SPLa
        ROI{1,7}={'L_MIP_ROI','L_VIP_ROI','L_7PL_ROI'}; % SPLp
        ROI{1,8}={'L_IPS1_ROI','L_V3A_ROI','L_V3B_ROI','L_V7_ROI','L_V6A_ROI'}; % DSVC ,'L_V6_ROI'
        ROI{1,9}={'L_LO1_ROI','L_LO2_ROI','L_LO3_ROI','L_V3CD_ROI','L_V4t_ROI','L_FST_ROI','L_MT_ROI','L_MST_ROI','L_PH_ROI'}; % MT+
        ROI{1,10}={'L_FFC_ROI','L_VVC_ROI','L_V8_ROI','L_VMV1_ROI','L_VMV2_ROI','L_VMV3_ROI','L_PIT_ROI'}; % VSVC
        ROI{1,11}={'L_RI_ROI','L_MBelt_ROI','L_PBelt_ROI','L_A1_ROI','L_LBelt_ROI'}; % EAC
        
        ROI_name={'','SMA','PMv','PMd','M1','S1','SPLa','SPLp','DSVC','MT+','VSVC','EAC'};  % important! add the first ''
                
        pathtosurf=fullfile(atlasDir,sprintf('FS_LR_%s',surf));
        P_glasser=gifti(fullfile(pathtosurf,sprintf('Glasser_2016.%sk.%s.label.gii',surf,hem{h})));
        P_brodmann=gifti(fullfile(pathtosurf,sprintf('ROI.%sk.%s.label.gii',surf,hem{h})));
        P=zeros(size(P_glasser.cdata));
        
        for i=1:length(ROI_name)-1
            if (i==4||i==5) % M1 and S1
                if i==4
                    P(ismember(P_brodmann.cdata,2),1)=i;
                else
                    P(ismember(P_brodmann.cdata,1),1)=i;
                end
            else
                roi_num=[];
                for j=1:length(ROI{1,i})
                    roi_num=cat(2,roi_num,P_glasser.labels.key(strcmp(P_glasser.labels.name,ROI{1,i}{j})));
                end
                P(ismember(P_glasser.cdata,roi_num),1)=i;
            end
        end
        
        
        % create a label
        G=surf_makeLabelGifti(P,'labelNames',ROI_name);
        groupDir=fullfile(baseDir,wbDir,sprintf('group%sk',surf));
        file=fullfile(groupDir,sprintf('ROI.%s.%s.label.gii',hem{h},'SSS'));
        save(G,file);
        fprintf(1,'Done.\n');
        
        varargout={P,ROI_name};        
    case 'ROI:redefine' % use this... %rdm, glm3
        % fmask_data.mat file is necessary, imana_sss('WB:vol2surf_resample','sn',s)? ->
        hemis=1;
        % sn=subj_vec;
        glm=0;  % change this glm=1 for S01-S06
        surf='32';        
        [P,ROI_name]=sss_imana('ROI:findall'); 
        PP = load('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/surfaceWB/group32k/fmask_data.mat');
        P = (PP.mask'~=0).*P;
        vararginoptions(varargin,{'sn','glm'});
        n_roi=length(ROI_name)-1;
        
        % loop over subjects
        % for s=sn
        subj_id = char(pinfo.subj_id(pinfo.sn==sn));
        S_id = strrep(subj_id,'R','S');
        fprintf('%s...\n',subj_id);
        R=cell(1,1);
        surfDir=fullfile(baseDir, wbDir,S_id);
        
        idx=0;
        
        for h=hemis
            file=fullfile(baseDir,sprintf(glmDir,glm),subj_id,'mask.nii');
            pathtosurf=fullfile(atlasDir,sprintf('FS_LR_%s',surf));
            surface=fullfile(pathtosurf,sprintf('fs_LR.%sk.%s.flat.surf.gii',surf,hem{h}));
            for r=1:n_roi
                idx=idx+1;
                R{idx}.type     = 'surf_nodes';
                R{idx}.white    = fullfile(surfDir,sprintf('%s.%s.white.%sk.surf.gii',S_id,hem{h},surf));
                R{idx}.pial     = fullfile(surfDir,sprintf('%s.%s.pial.%sk.surf.gii',S_id,hem{h},surf));
                R{idx}.flat     = surface;
                R{idx}.linedef  = [5,0,1];                             % take 5 steps along node between white (0) and pial (1) surfaces
                R{idx}.image    = file;
                R{idx}.name     = [subj_id '_' ROI_name{r+1} '_' hem{h}];
                R{idx}.location = find(P==r);
            end
        end
        R=region_calcregions(R);
        savename=fullfile(baseDir,roiDir,sprintf('%s_Task_regions.glm_%d.mat',subj_id,glm));
%       savename=fullfile(baseDir,roiDir,sprintf('%s_SSS_regions.mat',sprintf('S%02d',sn)));

        save(savename,'R');
        
        fprintf('\nROIs have been defined for %s \n',subj_id);
        clear R
        % end  


    case 'ROI:getpsc' % use this... %rdm, glm3
        % sss_imana('GLM:psc','sn',s), sss_imana('') ->
        %glm=1;
%        sn=subj_vec;
        vararginoptions(varargin, {'sn','glm'});

        cwd=pwd;
        
        T=[];
        % harvest
        for s=sn % for each subj
            fprintf('\nSubject: s%02d\n', s) % output to user
            
            % change directory to subject glm
            cd(fullfile(baseDir,sprintf(glmDir,glm),sprintf('S%02d',s)))
%             temp = dir('psc_Transition*nii');
%             O = {};            for i=1:length(temp) O{i}=temp(i).name;end;
%             cond = [1:64]';
%             O = {'psc_Motor2-1.nii','psc_Cue2-1.nii','psc_BothLetter2-1.nii',...
%                 'psc_BothSpatial2-1.nii','psc_CueLetter2-1.nii','psc_CueSpatial2-1.nii',...
%                 'psc_NRepMotor2-1.nii','psc_NRepCue2-1.nii','psc_NRep2-1.nii',...
%                 'psc_Letter.nii','psc_Spatial.nii'};
%             O = {'psc_MotorR.nii','psc_MotorN.nii','psc_MotorR-MotorN.nii',...
%                 'psc_CueR.nii','psc_CueN.nii','psc_CueR-CueN.nii',...
%                 'psc_BothR.nii','psc_BothN.nii','psc_BothR-BothN.nii',...
% %                 };
%             O = {'psc_MotorR.nii','psc_MotorN.nii',...
%                 'psc_CueR.nii','psc_CueN.nii',...
%                 'psc_BothR.nii','psc_BothN.nii',...
%                 };
%             O = {'psc_NRep-L.nii','psc_MotorRep-L.nii','psc_CueRep-L.nii','psc_BothRep-L.nii',...
%                 'psc_NRep-S.nii','psc_MotorRep-S.nii','psc_CueRep-S.nii','psc_BothRep-S.nii'};
             O = {'psc_Trial-state 1.nii','psc_Trial-state 2.nii','psc_Trial-state 3.nii','psc_Trial-state 4.nii',...
                 'psc_Trial-state 5.nii','psc_Trial-state 6.nii','psc_Trial-state 7.nii','psc_Trial-state 8.nii',...
                 'psc_Trial-state 9.nii'};
             SeqType = [1 1 1 1 2 2 2 2]';
             RepType = [1 2 3 4 1 2 3 4]';
% load ROI
            load(fullfile(baseDir,roiDir,sprintf('%s_%s_regions.mat',sprintf('S%02d',s),'Task')));
%             load(fullfile(baseDir,roiDir,sprintf('%s_%s_regions.mat',sprintf('S%02d',s),'SSS')));

            %cond = [1 2 3 1 2 3]';
%             cond = [1 2 3 4 5 6 7 8 9 10 11]';
%             cond = [1 2 3 4 5 6]';
%             cond = [1 2]';
            V=spm_vol(char(O));            
            % get raw data for voxels in region
            for r=1:11
%             for r=1:length(R)
            %             for r=1:length(R) % for each region
                
                Y=region_getdata(V,R{r});  
                S.psc=nanmean(Y,2);  % use nanmean, SKim
                S.nanperc = (length(find(isnan(Y)==1))/prod(size(Y)))*ones(length(O),1);
                S.hemi=repmat(hemi,length(O),1);
                S.roi=repmat(r,length(O),1);
                S.SN=repmat(s,length(O),1);
                S.SeqType= SeqType;
                S.RepType = RepType;
                S.cond = [1:8]';

                T=addstruct(T,S);
                fprintf('%d.',r)
            end
        end
        % save T
%         savename=fullfile(baseDir,roiDir,sprintf('psc_glm%d_%s_N=%d.mat',glm,'SSS',numel(sn)));

        savename=fullfile(baseDir,roiDir,sprintf('psc_glm%d_%s_N=%d.mat',glm,'Task',numel(sn)));
        save(savename,'-struct','T');
        fprintf('\n')
    case 'ROI:plotSSS_avg' % use this...
        %glm=1;
%        sn=subj_vec;
        vararginoptions(varargin, {'sn','glm','ptype'});
        
        % load
        loadname=fullfile(baseDir,roiDir,sprintf('psc_glm%d_%s_N=%d.mat',glm,'Task',numel(sn)));
        T=load(loadname);
        
        T=normData(T,{'psc'},'sub');
        
        CAT.facecolor={
            [0.3 0.3 0.6],...
            [0.3 0.6 0.3],...
            [0.6 0.3 0.3],...
            [0.3 0.3 0.3]
            };
        CAT.edgecolor = 'none';
        
                  % 1     2     3    4    5     6      7      8     9
%         ROI_name={'SMA','PMv','PMd','M1','S1','SPLa','SPLp','DSVC','MT+','VSVC'};
        ROI_name={'SMA','PMv','PMd','M1','S1','SPLa','SPLp','DSVC','MT+','VSVC','EAC'};
%          ROI_name={'SMA','PMv','PMd','M1','S1','SPLa','SPLp'};

        % region interaction
      %  anovaMixed(T.psc,T.SN,'within',[T.cond T.roi],{'cond','roi'},'subset',T.cond~=2&ismember(T.roi,[3 6]));
        
        figure()
        set(gcf,'color','w');
%         x_coord = barplot([T.roi],T.normpsc,'split',[T.cond],'CAT',CAT,'gapwidth',[0.8 0. 0 0],'barwidth',1,...
%             'leg',{'Motor','Cue'},'subset',T.cond<3);
% 
%             % 'leg',{'Letter','Spatial'},'subset',T.cond~=2&T.roi~=8);
%         legend('FontSize',fsl*2,'FontName',fontname);        
%         set(gca,'YLim',[-0.25 0.25],'TickLength',[0.01 0.01],'YTick',[-0.2:0.1:0.2],...
%             'XTick',(x_coord(1:2:end)+x_coord(2:2:end))/2,'XTickLabel',ROI_name,...
%             'XLim',[x_coord(1)-1.2 x_coord(end)+1.2],'FontSize',fs*2,'LineWidth',1,'FontName', fontname);
%         ylabel('RS (exe2-exe1): percent signal change','FontSize',fsl*2,'FontName',fontname);
%         print(gcf,fullfile(pathToSave,'ROI_MotorvsCue'),sprintf('-d%s','svg'),'-r1000');
%         save(fullfile(pathToSave,'ROI_MotorvsCue.png'));
        if ptype==1
            legend_text = {'NRep','Seq-Rep','Cue-Rep','Both-Rep'};
            title_text = 'Repetition effects';
            x_coord = barplot([T.roi],T.psc,'split',[T.RepType],'CAT',CAT,'gapwidth',[0.8 0. 0 0],'barwidth',1,...
                'leg',legend_text,'subset',T.SeqType==2);
        elseif ptype==2
            legend_text = {'Letter','Spatial'};
            title_text = 'RS effects';
            x_coord = barplot([T.roi],T.normpsc,'split',[T.SeqType],'CAT',CAT,'gapwidth',[0.8 0. 0 0],'barwidth',1,...
                'leg',legend_text);            
        elseif ptype==3
            legend_text = {'Letter','Spatial'};
            title_text = 'Repetition effects (Cue)';
            x_coord = barplot([T.roi],T.normpsc,'split',[T.cond],'CAT',CAT,'gapwidth',[0.8 0. 0 0],'barwidth',1,...
                'leg',legend_text,'subset',T.cond>=5 & T.cond<=6);    
        elseif ptype==4
            legend_text = {'NRep-Motor','NRep-Cue','NRep-Both'};
            title_text = 'Non-repetition effects';
            x_coord = barplot([T.roi],T.normpsc,'split',[T.cond],'CAT',CAT,'gapwidth',[0.8 0. 0 0],'barwidth',1,...
                'leg',legend_text,'subset',T.cond>=7 & T.cond<=9);              
        elseif ptype==5
            legend_text = {'Letter','Spatial'};
            title_text = 'Non-repetition effects';
            x_coord = barplot([T.roi],T.normpsc,'split',[T.cond],'CAT',CAT,'gapwidth',[0.8 0. 0 0],'barwidth',1,...
                'leg',legend_text,'subset',T.cond>=10 & T.cond<=11);              
        elseif ptype==6
            legend_text = {'Repetition','Non-repetition'};
            title_text = 'Sequence repetition effects';
            x_coord = barplot([T.roi],T.normpsc,'split',[T.cond],'CAT',CAT,'gapwidth',[0.8 0. 0 0],'barwidth',1,...
                'leg',legend_text,'subset',T.cond>=1 & T.cond<=2);  
        elseif ptype==7
            legend_text = {'Repetition','Non-repetition'};
            title_text = 'Cue repetition effects';
            x_coord = barplot([T.roi],T.normpsc,'split',[T.cond],'CAT',CAT,'gapwidth',[0.8 0. 0 0],'barwidth',1,...
                'leg',legend_text,'subset',T.cond>=3 & T.cond<=4);    
        elseif ptype==8
            legend_text = {'Both_Rep','Non_Rep'};
            title_text = 'Sequence and Cue repetition effects';
            x_coord = barplot([T.roi],T.normpsc,'split',[T.cond],'CAT',CAT,'gapwidth',[0.8 0. 0 0],'barwidth',1,...
                'leg',legend_text,'subset',T.cond>=5 & T.cond<=6);    
        end
        L = legend; L.AutoUpdate = 'off';

        % legend(legend_text,'FontSize',fsl,'FontName',fontname);        
        set(gca,'YLim',[0 2],'TickLength',[0.01 0.01],'YTick',[0:0.5:2],...
            'XTick',(x_coord(1:2:end)+x_coord(2:2:end))/2,'XTickLabel',ROI_name,...
            'XLim',[x_coord(1)-1.2 x_coord(end)+1.2],'FontSize',fs,'LineWidth',1,'FontName', fontname);
%         ylabel('RS (exe2-exe1): percent signal change','FontSize',fsl,'FontName',fontname);
        ylabel('Percent signal change','FontSize',fsl,'FontName',fontname);

        title(title_text,'Fontsize',fsl*2);
        drawline(0,'dir','horz'); 
        print(gcf,fullfile(pathToSave,title_text),sprintf('-d%s','svg'),'-r1000');
        % save(fullfile(pathToSave,title_te));
        
        set(gcf,'Units','inches','PaperUnits','inches')
        set(gcf,'PaperPosition',[2 2 3.4 2.5]);% 2.2
        %     [A,~,C]=pivottable(T.SN,[T.roi T.cond],T.psc,'mean','subset',T.cond~=2);
%        roi=6;
%        ttest(A(:,C(:,2)==1&C(:,1)==roi),A(:,C(:,2)==3&C(:,1)==roi),2,'paired');
   case 'ROI_get_prewhitened_beta'         % ROI MVA 1: Extracts raw timeseries and get prewhitened beta %rdm, glm3
        % sss_imana('ROI:redefine','sn',s,'glm',glm) ->
        glm = 3;
        roi = [1:11];
        vararginoptions(varargin, {'sn','glm','roi'});
        subj_id = char(pinfo.subj_id(pinfo.sn==sn));
        S_id = strrep(subj_id,'R','S');
        load(fullfile(baseDir, sprintf(glmDir,glm),subj_id,'SPM.mat'));
%         load(fullfile(baseDir,roiDir,sprintf('%s_SSS_regions.mat',sprintf('S%02d',sn))));
        load(fullfile(baseDir,roiDir,sprintf('%s_Task_regions.glm_%d.mat',subj_id,glm)));

        data = region_getdata(SPM.xY.VY,R(roi));
        
        Data = [];
        
        for r=roi
            % Get pre-whitened beta
            %B = mva_prewhiten_beta(data,SPM); 
            normmode = 'runwise';
            [B,~,Sw] = rsa.spm.noiseNormalizeBeta(data{r},SPM,'normmode',normmode);
            regNInt = [9:9:72 73:80]; % beta 80개. [9:9:72 73-80] 로 수정. 뒤에 8개 constant
            regInt = [1:80];regInt(regNInt) = []; % [1:80] 로 수정
            %regNInt = [9:9:144 145:160];
            %regInt = [1:160];regInt(regNInt) = [];
            S.beta = {B};
%             S.beta = {B(regInt,:)};
%             S.beta_nointerest = {B(regNInt,:)};
            S.glm = glm;
            S.subj = sn;
            S.numVox = size(data{r},2);
            S.name = {R{r}.name};
            S.normmode = normmode;
            Data = addstruct(Data,S);
        end
        
%         save(fullfile(baseDir,'patterns',sprintf('%s_ROI_pwhBeta.glm%d.mat',sprintf('S%02d',sn), glm)),'Data');
        save(fullfile(baseDir,'patterns',sprintf('%s_fROI_pwhBeta.glm%d.mat',subj_id, glm)),'Data');

%        varargout = Data;
   case 'ROI_calc_rdm'         % ROI MVA 1: Extracts raw timeseries and get prewhitened beta %rdm, glm3
        % sss_imana('ROI_get_prewhitened_beta','sn',s,'glm',glm) ->
        glm = 3;
        roi = [1:11];
%         normmode = 'runwise';
%         normmethod = 'multivariate';
        vararginoptions(varargin, {'sn','glm','roi'});
        subj_id = char(pinfo.subj_id(pinfo.sn==sn));
        load(fullfile(baseDir, sprintf(glmDir,glm),subj_id,'SPM.mat'));
%         load(fullfile(baseDir,roiDir,sprintf('%s_SSS_regions.mat',sprintf('S%02d',sn))));
        load(fullfile(baseDir,roiDir,sprintf('%s_Task_regions.glm_%d.mat',subj_id,glm)));
        data = region_getdata(SPM.xY.VY,R(roi));
       
        Data = [];
%        condvec = [repmat([1:8 0],1,16) zeros(1,16)]'; % should be a column vector
        condvec = [repmat([1:8 0],1,8) zeros(1,8)]';
%        partvec = [reshape(repmat([1:16],9,1),1,144) zeros(1,16)]'; % should be a column vector
        partvec = [reshape(repmat([1:8],9,1),1,72) zeros(1,8)]'; % beta 80개 run 1,...,1:10개, 2,...,2:10개, ...
%         condvec = repmat([1:8],1,16);
%         partvec = reshape(repmat([1:16],8,1),1,128);

%         load(fullfile(baseDir,'patterns',sprintf('%s_ROI_pwhBeta.glm%d.mat',sprintf('S%02d',sn), glm)));
        load(fullfile(baseDir,'patterns',sprintf('%s_fROI_pwhBeta.glm%d.mat',subj_id, glm)));
        for r=roi
            RDM  = rsa.distanceLDC(Data.beta{r},partvec,condvec,SPM.xX.xKXs.X); 
%             RDM = rsa.spm.distanceLDCraw(data{r},SPM,condvec,'normmode','runwise','normmethod','multivariate');
            S.glm = glm;
            S.subj = sn;
            S.numVox = size(data{r},2);
            S.name = {R{r}.name};
%             S.normmode = normmode;
%             S.normmethod = normmethod;
            S.RDM = {RDM};
            Data = addstruct(Data,S);      
        end
%         save(fullfile(baseDir,'patterns',sprintf('%s_RDM.glm%d.mat',sprintf('S%02d',sn), glm)),'Data');
        save(fullfile(baseDir,'patterns',sprintf('%s_fRDM.glm%d.mat',subj_id, glm)),'Data');
       
end





%% Added 
%% ------------------------- Functions -------------------------------------

function calc_PSC_formula = pscFormula(num_run)
    if num_run < 2
        error('num_run must be 2 or greater.');
    end

    numerator = sprintf('i%d', num_run + 1);
    variables = cell(1, num_run);
    for i = 1:num_run
        variables{i} = sprintf('i%d', i);
    end

    denominator = sprintf('(%s)', strjoin(variables, '+'));

    calc_PSC_formula = sprintf('%s./(%s/%d)', numerator, denominator, num_run);


function [et1, et2, tert] = spmj_et1_et2_tert(dataDir, subj_name)
    
    % vararginoptions(varargin,{'sn'});
    
    fmapDir = fullfile(dataDir, "BIDS", ['sub-' subj_name], "fmap");
    funcDir = fullfile(dataDir, "BIDS", ['sub-' subj_name], "func");
    
    phasediff = jsondecode(fileread(fullfile(fmapDir, ['sub-' subj_name '_phasediff.json'])));
    func1st = jsondecode(fileread(fullfile(funcDir, ['sub-' subj_name '_task-task_run-01_bold.json'])));
    
    %% compute tert
    % total EPI readout time = = echo spacing (in ms) * base resolution 
    % (also knows as number of echos). If you use GRAPPA acceleration, 
    % you need to divide the total number of echos by two:
    base_resolution = func1st.BaseResolution;
    echo_spacing = func1st.EffectiveEchoSpacing * 1000;                        % compute echo spacing in milliseconds
    tert = echo_spacing * base_resolution;
    
    %% retrieve et1 et2 (in milliseconds)
    et1 = phasediff.EchoTime1 * 1000;
    et2 = phasediff.EchoTime2 * 1000;

