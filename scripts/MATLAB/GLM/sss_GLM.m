function varargout = sss_GLM(what,varargin)
%%
% After mounting the diedrichsen datashare on a mac computer.

if ismac
    workdir = '/Volumes/Diedrichsen_data$/data/SeqSpatialSupp_fMRI';
    atlasDir = '/Volumes/Diedrichsen_data$/data/Atlas_templates';
    dir_git = '/Users/sungbeenpark/github';
    % dir_tool = '/Users/sungbeenpark/Documents/MATLAB';
% After mounting the diedrichsen datashare on the CBS server.
elseif isunix
    workdir = '/srv/diedrichsen/data/SeqSpatialSupp_fMRI';
    atlasDir = '/srv/diedrichsen/data/Atlas_templates';
    standardmeshDir = '/home/ROBARTS/skim2764/imaging_tools/surfAnalysis/standard_mesh';
elseif ispc
    workdir = 'F:/SeqSpatialSupp_fMRI';
    atlasDir = 'D:/mobaxterm/sungbeenpark/github/diedrichsenlab/atlas';
    standardmeshDir = 'D:/mobaxterm/sungbeenpark/github/surfAnalysis/standard_mesh';
    dir_git = 'D:/mobaxterm/sungbeenpark/github';
    % dir_tool = dir_git;
    
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
glmDir = 'glm_%d';
roiDir = 'ROI';

addpath(genpath(dir_git));
% addpath(genpath(dir_tool));

%% subject info

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

fsl=16;
fs=13;
fontname='Arial';

%% MAIN OPERATION 
switch(what)
    case 'GLM:design'
        %% dependency:
            % https://github.com/spm/spm.git
            % https://github.com/jdiedrichsen/rwls.git
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

        fname = fullfile(dir_git,'SeqSpatialSupp_fMRI/participants.tsv');
        [subj_id, S_id] = get_id(fname, sn);
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
        fname = fullfile(dir_git,'SeqSpatialSupp_fMRI/participants.tsv');
        [subj_id, S_id] = get_id(fname, sn);
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
        fname = fullfile(dir_git,'SeqSpatialSupp_fMRI/participants.tsv');
        [subj_id, S_id] = get_id(fname, sn);

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
        fname = fullfile(dir_git,'SeqSpatialSupp_fMRI/participants.tsv');
        [subj_id, S_id] = get_id(fname, sn);
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

    case 'WB:vol2surf_indiv' % map indiv vol contrasts (.nii) onto surface (.gifti)
            % sss_imana('GLM:psc','sn',s) ->
            % glm = 3;
            surf = '32k';
            hemis = [1 2];
            hem = {'L','R'};
            % map = 'psc'; % 't'; 'con'; 'dist'; 'psc';
            % vararginoptions(varargin,{'sn', 'glm', 'hemis', 'map', 'surf'});
            vararginoptions(varargin,{'sn', 'glm' 'map'});
            fname = fullfile(dir_git,'diedrichsenlab/SeqSpatialSupp_fMRI/participants.tsv');
            [subj_id, S_id] = get_id(fname, sn);
            for h = hemis
                surfDir = fullfile(baseDir, wbDir);
                white = fullfile(surfDir,S_id,sprintf('%s.%s.white.%s.surf.gii', S_id, hem{h}, surf));
                pial = fullfile(surfDir,S_id,sprintf('%s.%s.pial.%s.surf.gii', S_id, hem{h}, surf));
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
                    case 'beta' % beta maps (univariate GLM)
                        fnames      = cell(1,numel(SPM.Vbeta));
                        con_name    = cell(1,numel(SPM.Vbeta));
                        for j=1:numel(fnames)
                            fnames{j}   = fullfile(subjGLM, SPM.Vbeta(j).fname);
                            con_name{j} = SPM.Vbeta(j).fname;
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
                outfile = fullfile(surfDir,sprintf('glm%d',glm),sprintf('%s.%s.glm%d.%s.func.gii',subj_id,hem{h},glm,map));
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
   case 'ROI_get_prewhitened_beta'
       % Mahalanobis distance를 구하기 위해 (spatial) prewhitening 한 beta를 계산.
        % sss_imana('ROI:redefine','sn',s,'glm',glm) ->
        glm = 3;
        roi = [1:11];
        vararginoptions(varargin, {'sn','glm','roi'});
        fname = fullfile(dir_git,'diedrichsenlab/SeqSpatialSupp_fMRI/participants.tsv');
        [subj_id, S_id] = get_id(fname, sn);
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
            regInt = [1:80]; regInt(regNInt) = []; % [1:80] 로 수정
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
        save(fullfile(baseDir,'patterns',sprintf('glm%d/%s_fROI_pwhBeta.glm%d.mat',glm, subj_id, glm)),'Data');

%        varargout = Data;
   case 'ROI_calc_rdm'         % ROI MVA 1: Extracts raw timeseries and get prewhitened beta %rdm, glm3
        % sss_imana('ROI_get_prewhitened_beta','sn',s,'glm',glm) ->
        glm = 3;
        roi = [1:11];
%         normmode = 'runwise';
%         normmethod = 'multivariate';
        vararginoptions(varargin, {'sn','glm','roi'});
        fname = fullfile(dir_git,'diedrichsenlab/SeqSpatialSupp_fMRI/participants.tsv');
        [subj_id, S_id] = get_id(fname, sn);
        load(fullfile(baseDir, sprintf(glmDir,glm),subj_id,'SPM.mat'));
%         load(fullfile(baseDir,roiDir,sprintf('%s_SSS_regions.mat',sprintf('S%02d',sn))));
        load(fullfile(baseDir,roiDir,sprintf('%s_Task_regions.glm_%d.mat',subj_id,glm)));
        data = region_getdata(SPM.xY.VY,R(roi));
       
        Data = [];
%        condvec = [repmat([1:8 0],1,16) zeros(1,16)]'; % should be a column vector
        condvec = [repmat([1:8 0],1,8) zeros(1,8)]';
%        partvec = [reshape(repmat([1:16],9,1),1,144) zeros(1,16)]'; % should be a column vector
        % beta 80개 = run 1,...,1:9개 + 2,...,2:9개 + ... + 0, ..., 0:8개
        partvec = [reshape(repmat([1:8],9,1),1,72) zeros(1,8)]';
%         condvec = repmat([1:8],1,16);
%         partvec = reshape(repmat([1:16],8,1),1,128);

%         load(fullfile(baseDir,'patterns',sprintf('%s_ROI_pwhBeta.glm%d.mat',sprintf('S%02d',sn), glm)));
        load(fullfile(baseDir,'patterns',sprintf('glm%d/%s_fROI_pwhBeta.glm%d.mat', glm, subj_id, glm)));
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
        save(fullfile(baseDir,'patterns',sprintf('glm%d/%s_fRDM.glm%d.mat',glm, subj_id, glm)),'Data');
       
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


function [subj_id, S_id] = get_id(fname, sn)
    pinfo = dload(fname);
    subj_id = char(pinfo.subj_id(pinfo.sn==sn));
    S_id = strrep(subj_id,'R','S');