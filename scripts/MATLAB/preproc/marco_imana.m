function varargout = template_imana(what,varargin)
% Template function for preprocessing of the fMRI data.
% Rename this function to <experiment_name>_imana.m 
% Don't forget to add path the required tools!


% Directory specification

% Define the data base directory 

% automatic detection of datashare location:
% After mounting the diedrichsen datashare on a mac computer.
if isfolder('/Volumes/Diedrichsen_data$/data/<project_name>')
    workdir='/Volumes/Diedrichsen_data$/data/<project_name>';
% After mounting the diedrichsen datashare on the CBS server.
elseif isfolder('/cifs/diedrichsen/data/<project_name>')
    workdir='/cifs/diedrichsen/data/<project_name>';
else
    fprintf('Workdir not found. Mount or connect to server and try again.');
end

baseDir         = (sprintf('%s/',workdir));     % Base directory of the project
BIDS_dir        = 'BIDS';                       % Raw data post AutoBids conversion
behaviourDir    = 'behavioural_data';           % Timing data from the scanner
imagingRawDir   = 'imaging_data_raw';           % Temporary directory for raw functional data
imagingDir      = 'imaging_data';               % Preprocesses functional data
anatomicalDir   = 'anatomicals';                % Preprocessed anatomical data (LPI + center AC + segemnt)
fmapDir         = 'fieldmaps';                  % Fieldmap dir after moving from BIDS and SPM make fieldmap
suitDir         = 'suit';
regDir          = 'RegionOfInterest';

%% subject info

% Read info from participants .tsv file 
pinfo = dload(fullfile(baseDir,'participants.tsv'));

%% MAIN OPERATION 
switch(what)
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
            % BIDS transfers:
            template_imana('BIDS:move_unzip_raw_anat','sn',s);
            template_imana('BIDS:move_unzip_raw_func','sn',s);
            template_imana('BIDS:move_unzip_raw_fmap','sn',s);

            % ANAT functions:
            template_imana('ANAT:reslice_LPI','sn',s);
            template_imana('ANAT:centre_AC','sn',s);
            template_imana('ANAT:segmentation','sn',s);
            
            % FUNC functions:
            template_imana('FUNC:make_fmap','sn',s);
            template_imana('FUNC:realign_unwarp','sn',s,'rtm',0);
            template_imana('FUNC:move_realigned_images','sn',s,'prefix','u','rtm',0);
            template_imana('FUNC:meanimage_bias_correction','sn',s,'prefix','u','rtm',0);
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
            template_imana('FUNC:coreg','sn',s,'prefix','u','rtm',0);
            template_imana('FUNC:make_samealign','sn',s);
            template_imana('make_maskImage','sn',s);
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
        anat_raw_path = fullfile(baseDir,BIDS_dir,sprintf('sub-s%.02d',sn),'ses-01','anat',[char(pinfo.AnatRawName(pinfo.sn==sn)) '.nii.gz']);

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
        for sess = 1:pinfo.numSess(pinfo.sn==sn)
            
            % pull list of runs from the participant.tsv:
            run_list = eval(['pinfo.runsSess',num2str(sess),'(pinfo.sn==',num2str(sn),')']);
            run_list = str2double(split(run_list,'.'));

            % loop on runs of sess:
            for i = 1:length(run_list)
                
                % pull functional raw name from the participant.tsv:
                FuncRawName_tmp = eval(['pinfo.FuncRawNameSess',num2str(sess),'(pinfo.sn==',num2str(sn),')']);
                FuncRawName_tmp = [char(FuncRawName_tmp) '.nii.gz'];

                % add run number to the name of the file:
                FuncRawName_tmp = replace(FuncRawName_tmp,'XX',sprintf('%.02d',i));

                % path to the subj func data:
                func_raw_path = fullfile(baseDir,BIDS_dir,sprintf('sub-s%.02d',sn),sprintf('ses-%.02d',sess),'func',FuncRawName_tmp);
        
                % destination path:
                output_folder = fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn)),['sess' num2str(sess)]);
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
        end    

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
            fmapMagnitudeName_tmp = eval(['pinfo.fmapMagnitudeNameSess',num2str(sess),'(pinfo.sn==',num2str(sn),')']);
            magnitude = [char(fmapMagnitudeName_tmp) '.nii.gz'];
            
            fmapPhaseName_tmp = eval(['pinfo.fmapPhaseNameSess',num2str(sess),'(pinfo.sn==',num2str(sn),')']);
            phase = [char(fmapPhaseName_tmp) '.nii.gz'];

            % path to the subj fmap data:
            magnitude_path = fullfile(baseDir,BIDS_dir,sprintf('sub-s%.02d',sn),sprintf('ses-%.02d',sess),'fmap',magnitude);
            phase_path = fullfile(baseDir,BIDS_dir,sprintf('sub-s%.02d',sn),sprintf('ses-%.02d',sess),'fmap',phase);
    
            % destination path:
            output_folder = fullfile(baseDir,fmapDir,char(pinfo.subj_id(pinfo.sn==sn)),['sess' num2str(sess)]);
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
        R = V.mat(1:3,1:3);
        t = -1 * R * [pinfo.locACx(pinfo.sn==sn),pinfo.locACy(pinfo.sn==sn),pinfo.locACz(pinfo.sn==sn)];
        V.mat(1:3,4) = t;

        % writing the image with the changed header:
        spm_write_vol(V,dat);

    
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
            J.warp.write = [0 0];
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
        for sess = 1:pinfo.numSess(pinfo.sn==sn)
            % Prefix of the functional files:
            prefixepi  = '';

            % echo times of the gradient eho sequence:
            et1 = 4.92;
            et2 = 7.38;
    
            % total EPI readout time = = echo spacing (in ms) * base resolution 
            % (also knows as number of echos). If you use GRAPPA acceleration, 
            % you need to divide the total number of echos by two:
            tert = 90 * 0.7 / 2;

            % pull list of runs from the participant.tsv:
            run_list = eval(['pinfo.runsSess',num2str(sess),'(pinfo.sn==',num2str(sn),')']);
            run_list = split(run_list,'.');
            run_list = cellfun(@(x) sprintf('%.02d',str2double(x)), run_list, 'UniformOutput', false);

            subfolderFieldmap = sprintf('sess%d',sess);
            % function to create the makefieldmap job and passing it to the SPM
            % job manager:
            spmj_makefieldmap(baseDir,char(pinfo.subj_id(pinfo.sn==sn)),run_list, ...
                              'et1', et1, ...
                              'et2', et2, ...
                              'tert', tert, ...
                              'image', numDummys+1, ...
                              'prefix',prefixepi, ...
                              'rawdataDir',fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess)), ...
                              'subfolderFieldmap',subfolderFieldmap);
        end

    case 'FUNC:realign_unwarp'      
        % Do spm_realign_unwarp
        
        % handling input args:
        sn = [];
        rtm = 0;
        vararginoptions(varargin,{'sn','rtm'})
        if isempty(sn)
            error('FUNC:make_fmap -> ''sn'' must be passed to this function.')
        end

        % loop on sessions:
        for sess = 1:pinfo.numSess(pinfo.sn==sn)
            % Prefix of the functional files:
            prefixepi  = '';

            % pull list of runs from the participant.tsv:
            run_list = eval(['pinfo.runsSess',num2str(sess),'(pinfo.sn==',num2str(sn),')']);
            run_list = split(run_list,'.');
            run_list = cellfun(@(x) sprintf('%.02d',str2double(x)), run_list, 'UniformOutput', false);
            
            subfolderFieldmap = sprintf('sess%d',sess);
            spmj_realign_unwarp(baseDir,char(pinfo.subj_id(pinfo.sn==sn)),run_list, startTR, inf, ...
                                'prefix',prefixepi,...
                                'rawdataDir',fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess)),...
                                'subfolderFieldmap',subfolderFieldmap,...
                                'rtm',rtm);
        end

    case 'FUNC:realign'          
        % realign functional images
        % SPM realigns all volumes to the mean volume of first run
                
        for s = sn
            spm_jobman('initcfg')
            
            data = {};
                % initialize data cell array which will contain file names for runs/TR images
                func_ses_subj_dir = fullfile(imaging_dir ,subj_id);
                                
                for r = runs
                    % Obtain the number of TRs for the current run
                    for j = 1:numTRs - numDummys
                        data{r}{j,1} = fullfile(func_ses_subj_dir,sprintf('%s_run-%02d.nii,%d', subj_id, r,j));
                    end % j (TRs/images)
                end % r (runs)            
            spmj_realign(data);
            fprintf('- runs realigned for %s  ',subj_id);

        end % s (sn)
        
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
        for sess = 1:pinfo.numSess(pinfo.sn==sn)
            % pull list of runs from the participant.tsv:
            run_list = eval(['pinfo.runsSess',num2str(sess),'(pinfo.sn==',num2str(sn),')']);
            run_list = split(run_list,'.');
            run_list = cellfun(@(x) sprintf('%.02d',str2double(x)), run_list, 'UniformOutput', false);
            
            % loop on runs of the session:
            for r = 1:length(run_list)
                % realigned (and unwarped) images names:
                file_name = [prefix, char(pinfo.subj_id(pinfo.sn==sn)), '_run_', run_list{r}, '.nii'];
                source = fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),file_name);
                dest = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess));
                if ~exist(dest,'dir')
                    mkdir(dest)
                end
                dest = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),file_name);
                % move to destination:
                [status,msg] = movefile(source,dest);
                if ~status  
                    error('BIDS:move_realigned_images -> %s',msg)
                end

                % realign parameters names:
                source = fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),['rp_', char(pinfo.subj_id(pinfo.sn==sn)), '_run_', run_list{r}, '.txt']);
                dest = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),['rp_', char(pinfo.subj_id(pinfo.sn==sn)), '_run_', run_list{r}, '.txt']);
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
                source = fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),['mean', prefix, char(pinfo.subj_id(pinfo.sn==sn)), '_run_', run_list{1}, '.nii']);
                dest = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),['mean', prefix, char(pinfo.subj_id(pinfo.sn==sn)), '_run_', run_list{1}, '.nii']);
            else        % if registered to mean image of each run:
                source = fullfile(baseDir,imagingRawDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),[prefix, 'meanepi_', char(pinfo.subj_id(pinfo.sn==sn)), '.nii']);
                dest = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),[prefix, 'meanepi_', char(pinfo.subj_id(pinfo.sn==sn)), '.nii']);
            end
            % move to destination:
            [status,msg] = movefile(source,dest);
            if ~status  
                error('BIDS:move_realigned_images -> %s',msg)
            end
        end
    
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
        for sess = 1:pinfo.numSess(pinfo.sn==sn)
            % pull list of runs from the participant.tsv:
            run_list = eval(['pinfo.runsSess',num2str(sess),'(pinfo.sn==',num2str(sn),')']);
            run_list = split(run_list,'.');
            run_list = cellfun(@(x) sprintf('%.02d',str2double(x)), run_list, 'UniformOutput', false);

            if rtm==0   % if registered to first volume of each run:
                P{1} = fullfile(baseDir, imagingDir, char(pinfo.subj_id(pinfo.sn==sn)), sprintf('sess%d',sess), ['mean', prefix,  char(pinfo.subj_id(pinfo.sn==sn)), '_run_', run_list{1}, '.nii']);
            else        % if registered to mean image of each run:
                P{1} = fullfile(baseDir, imagingDir, char(pinfo.subj_id(pinfo.sn==sn)), sprintf('sess%d',sess), [prefix, 'meanepi_', char(pinfo.subj_id(pinfo.sn==sn)), '.nii']);
            end
            spmj_bias_correct(P);
        end

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
        for sess = 1:pinfo.numSess(pinfo.sn==sn)
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
            run_list = eval(['pinfo.runsSess',num2str(sess),'(pinfo.sn==',num2str(sn),')']);
            run_list = split(run_list,'.');
            run_list = cellfun(@(x) sprintf('%.02d',str2double(x)), run_list, 'UniformOutput', false);

            if rtm==0   % if registered to first volume
                mean_file_name = sprintf('rbmean%s%s_run_%s.nii', prefix, char(pinfo.subj_id(pinfo.sn==sn)), run_list{1});
            else    % if registered to the mean image
                mean_file_name = sprintf('rb%smeanepi_%s.nii', prefix, char(pinfo.subj_id(pinfo.sn==sn)));
            end
            J.source = {fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),mean_file_name)}; 
            J.ref = {fullfile(baseDir,anatomicalDir,char(pinfo.subj_id(pinfo.sn==sn)),[char(pinfo.subj_id(pinfo.sn==sn)), '_anatomical','.nii'])};
            J.other = {''};
            J.eoptions.cost_fun = 'nmi';
            J.eoptions.sep = [4 2];
            J.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            J.eoptions.fwhm = [7 7];
            matlabbatch{1}.spm.spatial.coreg.estimate=J;
            spm_jobman('run',matlabbatch);
            
            % (3) Check alignment manually by using fsleyes similar to step
            % one.
        end

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
        for sess = 1:pinfo.numSess(pinfo.sn==sn)
            % pull list of runs from the participant.tsv:
            run_list = eval(['pinfo.runsSess',num2str(sess),'(pinfo.sn==',num2str(sn),')']);
            run_list = split(run_list,'.');
            run_list = cellfun(@(x) sprintf('%.02d',str2double(x)), run_list, 'UniformOutput', false);
            
            % select the reference image:
            if rtm==0
                P{1} = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),['rb' 'mean' prefix char(pinfo.subj_id(pinfo.sn==sn)) '_run_' run_list{1} '.nii']);
            else
                P{1} = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),['rb' prefix 'meanepi_' char(pinfo.subj_id(pinfo.sn==sn)) '.nii']);
            end

            % select images to be realigned:
            Q = {};
            for r = 1:length(run_list)
                for i = 1:numTRs
                     Q{end+1} = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),[prefix char(pinfo.subj_id(pinfo.sn==sn)) '_run_' run_list{r} '.nii,' num2str(i)]);
                end
            end

            spmj_makesamealign_nifti(char(P),char(Q));
        end
    
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
        for sess = 1:pinfo.numSess(pinfo.sn==sn)
            run_list = eval(['pinfo.runsSess',num2str(sess),'(pinfo.sn==',num2str(sn),')']);
            run_list = split(run_list,'.');
            run_list = cellfun(@(x) sprintf('%.02d',str2double(x)), run_list, 'UniformOutput', false);
        
            % bias corrected mean epi image:
            if rtm==0
                nam{1} = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),['rb' 'mean' prefix char(pinfo.subj_id(pinfo.sn==sn)) '_run_' run_list{1} '.nii']);
            else
                nam{1} = fullfile(baseDir,imagingDir,char(pinfo.subj_id(pinfo.sn==sn)),sprintf('sess%d',sess),['rb' prefix 'meanepi_' char(pinfo.subj_id(pinfo.sn==sn)) '.nii']);
            end
            nam{2}  = fullfile(baseDir,anatomicalDir, subj_id{sn}, ['c1' subj_id{sn}, '_anatomical.nii']);
            nam{3}  = fullfile(baseDir,anatomicalDir, subj_id{sn}, ['c2' subj_id{sn}, '_anatomical.nii']);
            nam{4}  = fullfile(baseDir,anatomicalDir, subj_id{sn}, ['c3' subj_id{sn}, '_anatomical.nii']);
            spm_imcalc(nam, 'rmask_noskull.nii', 'i1>1 & (i2+i3+i4)>0.2')
            
            source = fullfile(baseDir,imagingDir,subj_id{sn}, 'rmask_noskull.nii');
            dest = fullfile(baseDir,anatomicalDir,subj_id{sn},'rmask_noskull.nii');
            movefile(source,dest);
            
            % gray matter mask for covariance estimation
            % ------------------------------------------
            nam={};
            nam{1}  = fullfile(imagingDir,subj_id{sn}, 'sess1', ['rb' prefix 'meanepi_' subj_id{sn} '.nii']);
            nam{2}  = fullfile(anatomicalDir, subj_id{sn}, ['c1' subj_id{sn}, '_anatomical.nii']);
            spm_imcalc(nam, 'rmask_gray.nii', 'i1>1 & i2>0.4')
            
            source = fullfile(imagingDir,subj_id{sn}, 'rmask_gray.nii');
            dest = fullfile(anatomicalDir,subj_id{sn},'rmask_gray.nii');
            movefile(source,dest);
        end


    case 'GLM:make_glm_1'    % design glm
        % make the design matrix for the glm
        % models each condition as a separate regressors
        % For conditions with multiple repetitions, one regressor
        % represents all the instances
        
        sn = [1:length(pinfo.participant_id)];
        hrf_cutoff = Inf;
        prefix = 'r'; % prefix of the preprocessed epi we want to use
        glm = 1;
        vararginoptions(varargin, {'sn', 'hrf_cutoff', 'ses'});
        

        % get the info file that specifies the the tasks and order?
        Dd = dload(fullfile(base_dir, 'task_description.tsv'));
        
        for s = sn
                func_subj_dir = fullfile(base_dir, func_dir,subj_str{s});
 
                % loop through runs within the current sessions
                itaskUni = 0;
                for ses = [1]
                 % create a directory to save the design
                  subj_est_dir = fullfile(base_dir, glm_first_dir,subj_str{s}, sprintf('ses-%02d',ses));
                  dircheck(subj_est_dir)
                  
                  T = []; % task/condition + session + run info
                  J = []; % structure with SPM fields to make the design
                 
                  J.dir            = {subj_est_dir};
                  J.timing.units   = 'secs';
                  J.timing.RT      = 1.3;
                  J.timing.fmri_t  = 16;
                  J.timing.fmri_t0 = 8;
                  
                    % get the list of runs for the current session
                    runs = run_list{ses};
                    for run = 1:2 %length(runs)
                       %V = spm_vol(fullfile(base_dir,func_dir, subj_str{s},sprintf('ses-%02d', ses),sprintf('r%s_run-%02d.nii', subj_str{s}, run)));
                       %numTRs = length(V);
             
                       % get the path to the tsv file
                       tsv_path = fullfile(base_dir, func_dir,subj_str{s});
                       % get the tsvfile for the current run
                       D = dload(fullfile(tsv_path,sprintf('ses-%02d',ses), sprintf('run%d.tsv', run)));
                       
                       % Get the onset and duration of the last sentence
                       lastSentenceOnset = D.onset(end);
                       lastSentenceDuration = D.duration(end);
                       
                       % Convert the end time of the last sentence to TRs
                       endTimeInTRs = ceil((lastSentenceOnset + lastSentenceDuration) / J.timing.RT);


                       % Define scans up to the last sentence's end time
                       N = cell(endTimeInTRs - numDummys, 1);
                       
                       for i = 1:(endTimeInTRs - numDummys)
                           N{i} = fullfile(func_subj_dir, sprintf('ses-%02d', ses), sprintf('%s%s_run-%02d.nii, %d', prefix, subj_str{s}, run, i+numDummys)); % to exclude dummy volumes
                       end % i (image numbers)
                       J.sess(run).scans = N; % scans in the current runs
                        
                       % loop over trials within the current run and build up
                       % the design matrix
                       for ic = 1:length(Dd.task_name)
                           itaskUni = itaskUni+1;
                           % get the indices corresponding to the current
                           % condition.
                           % this line is necessary as there are some
                           % conditions with more than 1 repetition
                           idx = strcmp(D.trial_type, Dd.task_name{ic});
                           fprintf('* %d instances found for condition %s in run %02d\n', sum(idx), Dd.task_name{ic}, run)
                            
                           %
                           % filling in "reginfo"
                           TT.sn        = s;
                           TT.sess      = ses;
                           TT.run       = run;
                           TT.task_name = Dd.task_name(ic);
                           TT.task      = ic;
                           TT.taskUni   = itaskUni;
                           TT.n_rep     = sum(idx);
                            
                           % filling in fields of J (SPM Job)
                           J.sess(run).cond(ic).name = Dd.task_name{ic};
                           J.sess(run).cond(ic).tmod = 0;
                           J.sess(run).cond(ic).orth = 0;
                           J.sess(run).cond(ic).pmod = struct('name', {}, 'param', {}, 'poly', {});
                            
                           % get onset and duration (should be in seconds)
                           onset    = D.onset(idx) - (J.timing.RT*numDummys);
                           fprintf("The onset is %f\n", onset)
                           if onset < 0
                               warning("negative onset found")
                           end
                           duration = D.duration(idx);
                           fprintf("The duration is %f\n", duration);
                            
                           J.sess(run).cond(ic).onset    = onset;
                           J.sess(run).cond(ic).duration = duration;
                            
                           % add the condition info to the reginfo structure
                           T = addstruct(T, TT);
                            
                            
                        end % ic (conditions)
                        
                        % Regressors of no interest 
                       J.sess(run).multi     = {''};
                       J.sess(run).regress   = struct('name', {}, 'val', {});
                       J.sess(run).multi_reg = {''};
                       J.sess(run).hpf       = hrf_cutoff; % set to 'inf' if using J.cvi = 'FAST'. SPM HPF not applied
                   end % run (runs of current session)
                
                
               J.fact             = struct('name', {}, 'levels', {});
               J.bases.hrf.derivs = [0 0];
               J.bases.hrf.params = [4.5 11];                                  % set to [] if running wls
               J.volt             = 1;
               J.global           = 'None';
               J.mask             = {fullfile(func_subj_dir,'ses-01','rmask_noskull.nii')};
               J.mthresh          = 0.05;
               J.cvi_mask         = {fullfile(func_subj_dir, 'ses-01', 'rmask_gray.nii')};
               J.cvi              =  'fast';
                
               spm_rwls_run_fmri_spec(J);
                
                
               dsave(fullfile(J.dir{1},sprintf('%s_reginfo.tsv', subj_str{s})), T);
               fprintf('- estimates for glm_%d session %d has been saved for %s \n', glm, ses, subj_str{s});
             end % ses (session)
            
            
        end % sn (subject)   
    
    case 'GLM:estimate'      % estimate beta values
        
        sn       = subj_id; % subject list
        sessions   = [1];       % session number
        
        vararginoptions(varargin, {'sn', 'sessions'})
        
        for s = sn
         
            for ses = sessions
                fprintf('- Doing glm estimation for session %02d %s\n', ses, subj_str{s});
                subj_est_dir = fullfile(base_dir, glm_first_dir,subj_str{s}, sprintf('ses-%02d', ses));         
            
                load(fullfile(subj_est_dir,'SPM.mat'));
                SPM.swd = subj_est_dir;
            
                spm_rwls_spm(SPM);
            end
        end % s (sn),  
         
    case 'GLM:T_contrast'    % make T contrasts for each condition
        %%% Calculating contrast images.
        
        sn             = subj_id;    % subjects list
        ses            = 1;              % task number
        glm            = 1;              % glm number
        baseline       = 'rest';         % contrast will be calculated against base (available options: 'rest')
        
        vararginoptions(varargin, {'sn', 'glm', 'ses', 'baseline'})
        
        for s = sn
            
            % get the subject id folder name
            fprintf('Contrasts for session %02d %s\n', ses, subj_str{s})
            glm_dir = fullfile(base_dir, glm_first_dir, subj_str{s}, ses_str{ses}); 
            
            cd(glm_dir);
            
            % load the SPM.mat file
            load(fullfile(glm_dir, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
            T    = dload(fullfile(glm_dir, sprintf('%s_reginfo.tsv', subj_str{s})));
            
            % t contrast for each condition type
            utask = unique(T.task)';
            idx = 1;
            for ic = utask
                switch baseline
                    case 'myBase' % contrast vs future baseline :)))
                        % put your new contrasts here!
                    case 'rest' % contrast against rest
                        con                          = zeros(1,size(SPM.xX.X,2));
                        con(:,logical((T.task == ic)& (T.n_rep>0))) = 1;
%                         n_rep = length(T.run(T.task == ic));
%                         n_rep_t = T.n_rep(T.task == ic);
%                         name = unique(T.task_name(T.task == ic));
%                         fprintf('- task is %s: \n', name{1});
%                         fprintf('number of reps in all runs = %d\n', n_rep);
%                         fprintf('numberof reps recorded in tsv = %d\n', n_rep_t);
                        con                          = con/abs(sum(con));            
                end % switch base

                % set the name of the contrast
                contrast_name = sprintf('%s-%s', char(unique(T.task_name(T.task == ic))), baseline);
                SPM.xCon(idx) = spm_FcUtil('Set', contrast_name, 'T', 'c', con', SPM.xX.xKXs);
                
                idx = idx + 1;
            end % ic (conditions)
            
            SPM = spm_contrasts(SPM,1:length(SPM.xCon));
            save('SPM.mat', 'SPM','-v7.3');
            SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            save(fullfile(glm_dir, 'SPM_light.mat'), 'SPM')

            % rename contrast images and spmT images
            conName = {'con','spmT'};
            for i = 1:length(SPM.xCon)
                for n = 1:numel(conName)
                    oldName = fullfile(glm_dir, sprintf('%s_%2.4d.nii',conName{n},i));
                    newName = fullfile(glm_dir, sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                    movefile(oldName, newName);
                end % conditions (n, conName: con and spmT)
            end % i (contrasts)
        end % sn
    
    case 'SURF:reconall' % Freesurfer reconall routine
        % Calls recon-all, which performs, all of the
        % FreeSurfer cortical reconstruction process
        
        sn   = subj_id; % subject list
        
        vararginoptions(varargin, {'sn'});
        
        % Parent dir of anatomical images    
        for s = sn
            fprintf('- recon-all %s\n', subj_str{s});
                        % Get the directory of subjects anatomical;
            freesurfer_reconall(fs_dir, subj_str{s}, ...
                      fullfile(anatomical_dir, subj_str{s}, 'anatomical.nii'));
        end % s (sn)
        
    case 'SURF:fs2wb'          % Resampling subject from freesurfer fsaverage to fs_LR
        
        sn   = subj_id; % subject list
        res  = 32;          % resolution of the atlas. options are: 32, 164
        hemi = [1, 2];      % list of hemispheres
       
        vararginoptions(varargin, {'sn', 'res', 'hemi'});

        for s = sn 
            % get the subject id folder name
            outDir   = fullfile(baseDir, 'surfaceWB', 'data'); dircheck(outDir);
            surf_resliceFS2WB(subj_str{s}, fs_dir, outDir, 'hemisphere', hemi, 'resolution', sprintf('%dk', res))
        end % s (sn)  

end




%%  =======================Project-specific Cases==================================

switch(what)
    case 'SUIT:isolate_segment'  
    % Segment cerebellum into grey and white matter
    
        sn = subj_id;
        
        vararginoptions(varargin, {'sn'});
        
        for s = sn
            fprintf('- Isolate and segment the cerebellum for %s\n', ...
                subj_str{s})
            spm_jobman('initcfg')
            
            % Get the file of subjects anatomical
            anat_subj_dir  = fullfile(anatomical_dir, subj_str{s});
            anat_name = 'anatomical.nii'
    
            % Define suit folder
            suit_dir = fullfile(baseDir, 'suit/anatomicals',subj_str{s});
            % Create suit folder if it does not exist
            if ~exist(suit_dir, 'dir')
                mkdir (suit_dir)
            end
            
            % Copy anatomical_raw file to suit folder
            source = fullfile(anat_subj_dir, anat_name);
            dest   = fullfile(suit_dir, anat_name);           
            copyfile(source, dest);
            
            % go to subject directory for suit and isolate segment
            suit_isolate_seg({dest}, 'keeptempfiles', 1);
        end % s (sn)
    
    case 'SUIT:normalise_dartel' % SUIT normalization using dartel
        % LAUNCH SPM FMRI BEFORE RUNNING!!!!!
        sn = subj_id; %subjNum
        vararginoptions(varargin, 'sn');
        
        for s = sn
            suit_subj_dir = fullfile(baseDir, 'suit/anatomicals', subj_str{s});
            job.subjND.gray       = {fullfile(suit_subj_dir,'c_anatomical_seg1.nii')};
            job.subjND.white      = {fullfile(suit_subj_dir,'c_anatomical_seg2.nii')};
            job.subjND.isolation  = {fullfile(suit_subj_dir,'c_anatomical_pcereb.nii')};
            suit_normalize_dartel(job);
    
        end % s (subjects)

    case 'SUIT:save_dartel_def'    
        sn = subj_id; %subjNum
        % Saves the dartel flow field as a deformation file. 
        for s = sn
            cd(fullfile(baseDir,'suit/anatomicals', subj_str{s}));
            anat_name = 'anatomical';
            suit_save_darteldef(anat_name);
        end

    case 'SUIT:reslice'            % Reslice stuff into suit space 
        % run the case with 'anatomical' to check the suit normalization
        % make sure that you reslice into 2mm^3 resolution
        
        sn   = subj_id;
        type = 'con';  % 'betas' or 'con' or 'ResMS' or 'cerebellarGrey' or 'anatomical'
        mask = 'c_anatomical_pcereb'; % 'cereb_prob_corr_grey' or 'cereb_prob_corr' or 'dentate_mask' or 'pcereb'
        glm  = 1;             % glm number. Used for reslicing betas and contrasts 
        
        vararginoptions(varargin, {'sn', 'type', 'mask', 'glm'})
        
        for s = sn
            suit_dir = fullfile(baseDir, 'suit/anatomical',subj_str{s});
            switch type
                case 'anatomical'
                    subj_dir = suit_dir;
                    % Get the name of the anatpmical image
                    files2map = sprintf('%s_T1w_lpi.nii', subj_str{s});
                    
                    job.subj.resample = {sprintf('%s,1', files2map)};
                 case 'betas'
                    glmSubjDir = fullfile(glm_first_dir,sprintf('glm_%d',glm),subj_str{s});
                    images='resbeta_0';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                case 'con'
                    glmSubjDir = fullfile(glm_first_dir,sprintf('glm_%d',glm),subj_str{s});
                    images='con_';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                case 'spmT'
                    glmSubjDir = fullfile(glm_first_dir,sprintf('glm_%d',glm),subj_str{s});
                    images='spmT_';
                    source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
            end
            job.subj.affineTr = {fullfile(baseDir,'suit','anatomicals',subj_str{s},'Affine_c_anatomical_seg1.mat')};
            job.subj.flowfield= {fullfile(baseDir,'suit','anatomicals',subj_str{s},'u_a_c_anatomical_seg1.nii')};
            job.subj.resample = {source.name};
            job.subj.mask     = {fullfile(baseDir,'suit','anatomicals',subj_str{s},sprintf('%s.nii',mask))};
            job.vox           = [1 1 1];
            % Replace Nans with zeros to avoid big holes in the the data 
            for i=1:length(source)
                V=spm_vol(source(i).name); 
                X=spm_read_vols(V); 
                X(isnan(X))=0; 
                spm_write_vol(V,X); 
            end; 
            suit_reslice_dartel(job);
            
            source=fullfile(glmSubjDir,'*wd*');
            destination=fullfile(baseDir,'suit',sprintf('glm_%d',glm),subj_str{s});
            movefile(source,destination);
    
            fprintf('%s have been resliced into suit space \n',type)
        end





 





