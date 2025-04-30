function varargout = sss_imana(what,varargin)
%%
% After mounting the diedrichsen datashare on a mac computer.
if ismac
    rootdir='/Volumes/Diedrichsen_data$/data/SeqSpatialSupp_fMRI';
% After mounting the diedrichsen datashare on the CBS server.
elseif isunix
    rootdir='/srv/diedrichsen/data/SeqSpatialSupp_fMRI';
elseif ispc
    rootDir='F:\SeqSpatialSupp_fMRI';
    atlasDir='D:/mobaxterm/sungbeenpark/github/diedrichsenlab/atlas';
    dir_git = '\\wsl.localhost/ubuntu-22.04/home/sungbeenpark/github';
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
regDir          = 'RegionOfInterest';
freesurferDir   = 'freesurf';
wbDir           = 'surfaceWB';
glmDir          = '/glm_%d';
roiDir          = 'ROI';

if exist(dir_git, 'dir') && ~contains(path, dir_git)
    addpath(genpath(dir_git));
    addpath(genpath(fullfile('cifti_write-matlab/private')));
end
atlasDir = fullfile(dir_git,'SeqSpatialSupp_fMRI/atlas');

%% input arguments
sn = [];
rtm = 0;
prefix = 'u';
vararginoptions(varargin,{'sn', 'rtm', 'prefix'})
if isempty(sn)
    error('BIDS:move_unzip_raw_func -> ''sn'' must be passed to this function.')
end

%% subject info
% Read info from participants .tsv file
pinfo = dload(fullfile(dir_git,'SeqSpatialSupp_fMRI/participants.tsv'));
subj_id = char(pinfo.subj_id(pinfo.sn==sn));
S_id = strrep(subj_id,'R','S');
list_run = cellfun(@(x) sscanf(x, '%d.'), pinfo.runlist(sn), 'UniformOutput', false);
list_run = list_run{1};

% ROI
hemi = {'lh','rh'}; % left & right hemi folder names/prefixes
hem = {'L', 'R'}; % hemisphere: 1=LH 2=RH
hname = {'CortexLeft', 'CortexRight'}; % 'CortexLeft', 'CortexRight', 'Cerebellum'

% Scanner
TR = 1000;
numDummys  = 8;  % dummy images at the start of each run (these are discarded)
numTRs = 410; %% Adjust this for experiment
endTR = 410; % S.Park: prevent the error in sss_imana('FUNC:realign_unwarp','sn',s,'rtm',0)

% nTRs = [410*ones(1,10) 401 406 410 404 410 410 385]; % For S11
% nTRs = [410*ones(1,16) 385]; % for S09

%% MAIN OPERATION 
switch(what)
    case 'PREP:step0'           
        % All preprocessing steps by just one go (AC coordinates (loc_AC) are prerequisite)
        % handling input args:

        % BIDS transfers:
        sss_imana('BIDS:move_unzip_raw_anat');
        sss_imana('BIDS:move_unzip_raw_func');
        sss_imana('BIDS:move_unzip_raw_fmap');
        % ANAT functions:
        sss_imana('ANAT:reslice_LPI');
        %  after this step, find the AC coordinate and enter it into
        %  participants.tsv

        fprintf('Manually align the mean bias corrected image and the anatomical\n')
    
    case 'PREP:step1'           
        % All preprocessing steps by just one go (AC coordinates (loc_AC) are prerequisite)
  
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
        % loop on subjects and preprocessing the data:
        for s = sn
            % FUNC:
            sss_imana('FUNC:coreg','sn',s,'prefix','u','rtm',0);
            sss_imana('FUNC:make_samealign','sn',s);
            sss_imana('FUNC:make_maskImage','sn',s);
        end


    case 'BIDS:move_unzip_raw_anat'
        % Moves, unzips and renames raw functional (BOLD) images from 
        % BIDS directory. After you run this function you will find
        % nRuns Nifti files named <subj_id>_run_XX.nii in the 
        % <project_id>/imaging_data_raw/<subj_id>/ directory.

        fname_raw = fullfile(baseDir, BIDSDir,['sub-' subj_id], 'anat', [char(pinfo.AnatRawName(pinfo.sn==sn)) '.nii.gz']);

        dir_output = fullfile(baseDir,anatomicalDir,subj_id);
        fname_output = fullfile(dir_output,[subj_id '_anatomical_raw.nii.gz']);

        if ~exist(dir_output,"dir")
            mkdir(dir_output);
        end

        % copy file to destination:
        [status,msg] = copyfile(fname_raw,fname_output);
        if ~status  
            error('ANAT:move_anatomical -> subj %d raw anatmoical was not moved from BIDS to the destenation:\n%s',sn,msg)
        end

        % unzip the .gz files to make usable for SPM:
        gunzip(fname_fname_output);

        % delete the compressed file:
        delete(fname_fname_output);

    case 'BIDS:move_unzip_raw_func'
        % Moves, unzips and renames raw functional (BOLD) images from BIDS
        % directory
        
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
