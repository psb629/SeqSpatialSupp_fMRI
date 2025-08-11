function varargout = sss_hrf(what,varargin)

if ispc
    cd '\\wsl.localhost/ubuntu-24.04/home/sungbeenpark/github/SeqSpatialSupp_fMRI/scripts/MATLAB'
elseif ismac
    cd '/Users/sungbeenpark/github/SeqSpatialSupp_fMRI/scripts/MATLAB'
end

sss_init;

%% HRF parameters
hrf_params_default = [6 16 1 1 6 0 32];
hrf_params = [];

%% argument inputs
sn = [];
glm = [];
LR = 'L';
pre = 10;
post = 20;
run = 1;
roi = 'M1';
vararginoptions(varargin,{'sn','glm','hrf_params','LR','run','pre','post','roi'});
hrf_params = [hrf_params hrf_params_default(length(hrf_params)+1:end)];
if isempty(sn)
    error('GLM:design -> ''sn'' must be passed to this function.')
end
[subj_id, S_id] = get_id(fullfile(dir_git,'SeqSpatialSupp_fMRI/participants.tsv'), sn);

if isempty(glm)
    error('GLM:design -> ''glm'' must be passed to this function.')
end
glmDir = sprintf('glm_%d',glm);

%% ROI
hemi = {'lh','rh'}; % left & right hemi folder names/prefixes
hem = {'L', 'R'}; % hemisphere: 1=LH 2=RH
hname = {'CortexLeft', 'CortexRight'}; % 'CortexLeft', 'CortexRight', 'Cerebellum'

%% MAIN OPERATION 
switch(what)
    case 'HRF:all'
        

    case 'HRF:get_mean_ts'
        % ROI 내의 평균 time series를 GLM 모델과 비교하여 출력.
        DoSave = false;
        vararginoptions(varargin,{'sn','glm','LR','DoSave'});
        DoSave = logical(DoSave);

        [subj_id, S_id] = get_id(sn);
        sprintf('HRF:get_mean_ts: %s...\n',subj_id);
        %% load R.mat (ROI information)
        workDir = fullfile(baseDir,roiDir,sprintf('glm%d',glm),subj_id);
        fname = fullfile(workDir,sprintf('%s.Task_regions.glm%d.mat',subj_id,glm));
        R = load(fname); % load the variable R
        R = R.R;
        % [R, V] = sss_hrf('ROI:deform','sn',sn,'glm',glm,'LR',LR);

        %% load SPM.mat (GLM information)
        SPM = load(fullfile(baseDir,sprintf('glm_%d',glm),subj_id,'SPM.mat'));
        SPM = SPM.SPM;
        % Find onsets for all events
        [D, start_sess] = spmj_get_ons_struct(SPM);
        
        %% extract an average time series for each ROI
        fprintf('Extration Y_raw for each ROI...');
        [y_raw, y_adj, y_hat, y_res, B] = region_getts(SPM,R); fprintf('done!\n');
        
        %% arrange variables as a Q
        Q = {};
        Q.R = R;
        Q.Y = {};
        Q.Y.D = D;
        Q.Y.D.start_sess = start_sess;
        Q.Y.B = B;
        Q.Y.y_raw = y_raw;
        Q.Y.y_hat = y_hat;
        Q.Y.y_res = y_res;
        Q.Y.y_adj = y_adj;

        varargout = {Q};
        if DoSave
            fname = fullfile(dir_work,sprintf('%s.glm%d.%drois.mat',subj_id,glm,length(R)));
            save(fname,'Q','-v7.3');
        end

    % case 'HRF:fit_params'
    %     LR = 'L';
    %     vararginoptions(varargin,{'sn','glm','roi','LR'});
    % 
    %     [subj_id, S_id] = get_id(sn);
    %     dir_work = fullfile(baseDir,roiDir,sprintf('glm%d',glm),subj_id);
    % 
    %     %% load y_raw
    %     fname = fullfile(dir_work,sprintf('cifti.L.glm%d.%s.%s.y_raw.nii',glm,subj_id,roi));
    %     cii = cifti_read(fname);
    %     Yraw = double(cii.cdata(:,:)');
    %     clear cii
    % 
    %     %% load SPM.mat (GLM information)
    %     SPM = load(fullfile(baseDir,sprintf('glm_%d',glm),subj_id,'SPM.mat'));
    %     SPM = SPM.SPM;
    % 
    %     %% optimization
    %     % fit_method = 'gamma';
    %     fit_method = 'library';
    %     [SPM,Yhat,Yres,p_opt] = spmj_fit_hrfparams(SPM,Yraw,fit_method);
    % 
    %     varargout = {SPM,Yhat,Yres,p_opt};

    case 'HRF:example'
        dir_work = fullfile(baseDir,roiDir,glmDir,subj_id);

        %% load y_raw
        fname = fullfile(dir_work,sprintf('cifti.L.glm%d.%s.%s.y_raw.nii',glm,subj_id,roi));
        cii = cifti_read(fname);
        Yraw = double(cii.cdata(:,:)');
        clear cii

        %% load SPM.mat (GLM information)
        SPM = load(fullfile(baseDir,sprintf('glm_%d',glm),subj_id,'SPM.mat'));
        SPM = SPM.SPM;

        % Get the hemodynamic response in micro-time resolution
        SPM.xBF.UNITS    = 'secs'; % units of the hrf
        SPM.xBF.T        = 16; % microtime resolution: number of time samples per scan
        SPM.xBF.dt       = SPM.xY.RT/SPM.xBF.T;  % Delta-t per sample 
        SPM.xBF.T0       = 1; % first time bin
        SPM.xBF.name     = 'fitted_hrf';
        SPM.xBF.order    = 1;
        SPM.xBF.Volterra = 1;  % volterra expansion order?
        SPM.xBF.bf = spm_hrf(SPM.xBF.dt,hrf_params);
        % p(1) - delay of response (relative to onset)          6
        % p(2) - delay of undershoot (relative to onset)       16
        % p(3) - dispersion of response                         1
        % p(4) - dispersion of undershoot                       1
        % p(5) - ratio of response to undershoot                6
        % p(6) - onset {seconds}                                0
        % p(7) - length of kernel {seconds}                    32
        SPM.xBF.length = size(SPM.xBF.bf,1)*SPM.xBF.dt; % support in seconds 
        
        % Reconvolve the design matrix with new HRF
        SPM = spmj_glm_convolve(SPM);
        
        % Restimate the betas and get predicted and residual response
        [beta, Yhat, Yres] = spmj_glm_fit(SPM,Yraw);
        
        figure;
        sgtitle(sprintf('%s (%s)',subj_id,roi));
        % Diagnostic plots 
        subplot(2,2,1);
        t=[SPM.xBF.dt*SPM.xBF.T0:SPM.xBF.dt:SPM.xBF.length]; 
        plot(t,SPM.xBF.bf);
        drawline([0],'dir','horz','linestyle','-');
        title('Basis Function(s)'); 
        
        % i-th run, convolved design matrix, sum of regressors of interest
        subplot(2,2,2);
        X=SPM.xX.X(SPM.Sess(run).row,SPM.Sess(run).col);
        plot(sum(X,2));
        title(sprintf('Run %d - overall response',run));
        
        subplot(2,2,3);
        % 각 ROI별 복셀 평균을 하여 Y들을 구함
        Yhat_run = mean(Yhat(SPM.Sess(run).row,:),2);
        Yres_run = mean(Yres(SPM.Sess(run).row,:),2);
        Yadj_run = Yres_run+Yres_run; 
        t= SPM.Sess(run).row; 
        plot(t,Yhat_run,'b',t,Yadj_run,'b:'); 
        title(sprintf('Run %d - overall response',run));
        
        % Get onset structure, cut-out the trials of choice, and plot evoked response
        subplot(2,2,4);
        D = spmj_get_ons_struct(SPM);
        Yadj = Yhat+Yres;
        for i=1:size(D.block,1)
            % 각 ROI별 복셀 평균을 하여 Y들을 onset time (TR) 기준 [pre,post] 범위로 절단.
            % 그후 각 Block 별(열) 이어 붙임.
            D.y_adj(i,:)=cut(mean(Yadj,2),pre,round(D.ons(i))-1,post,'padding','nan')';
            D.y_hat(i,:)=cut(mean(Yhat,2),pre,round(D.ons(i))-1,post,'padding','nan')';
            D.y_res(i,:)=cut(mean(Yres,2),pre,round(D.ons(i))-1,post,'padding','nan')';
        end
        T = getrow(D,mod(D.num,2)==1); % Get the first onset for each double (홀수인 행 추출) 
        traceplot([-pre:post],T.y_adj,'errorfcn','stderr'); % ,
        hold on;
        traceplot([-pre:post],T.y_hat,'linestyle','--',...
                'linewidth',3); % ,
        drawline([-8 8 16],'dir','vert','linestyle',':');
        drawline([0],'dir','vert','linestyle','--');
        drawline([0],'dir','horz','linestyle','-');
        hold off;
        xlabel('TR');
        ylabel('activation');
        A = mean(D.y_adj,2);
        B = mean(D.y_hat,2);
        coef = corrcoef(A(~isnan(A)),B(~isnan(B)));
        D.y_res(isnan(D.y_res)) = 0;
        epsq = trace(D.y_res'*D.y_res);
        title(sprintf('%s\nr=%.3f, |e|^{2}=%g',mat2str(hrf_params),coef(1,2),epsq));
        
%     case 'HRF:ROI_hrf_get'  % Extract raw and estimated time series from ROIs
%         sn = [];
%         ROI = 'all';
%         pre=10;
% %        post=30;
%         atlas = 'SSS';
%         glm = 3; % change this glm=1 for S01-06
%         bf = [];
% %         vararginoptions(varargin,{'ROI','pre','post', 'glm', 'sn', 'atlas'});
%         vararginoptions(varargin,{'sn','glm','post','bf'});
%         glmDir = fullfile(baseDir,sprintf(glmDir,glm));
%         T=[];
% 
%         [subj_id, ~] = get_id(sn);
%         fprintf('%s\n',subj_id);
% 
%         % load SPM.mat
%         cd(fullfile(glmDir,subj_id));
%         SPM = load('SPM.mat'); SPM=SPM.SPM;
%         if ~isempty(bf)
%             SPM.xBF.bf = bf;
%             SPM = fMRI_design_changeBF(SPM);
%         end
% 
%         % load ROI definition (R)
%         % R = load(fullfile(baseDir, roiDir,[subj_id '_' atlas '_regions.mat']));
%         R = load(fullfile(baseDir,roiDir,sprintf('%s_Task_regions.glm_3.mat',subj_id)));
%         R=R.R;
% 
%         % extract time series data from GifTi
%         [y_raw, y_adj, y_hat, y_res, B] = region_getts(SPM,R);
% 
%         D = spmj_get_ons_struct(SPM);
% 
%         for r=1:size(y_raw,2)
%             for i=1:size(D.block,1)
%                 D.y_adj(i,:)=cut(y_adj(:,r),pre,round(D.ons(i))-1,post,'padding','nan')';
%                 D.y_hat(i,:)=cut(y_hat(:,r),pre,round(D.ons(i))-1,post,'padding','nan')';
%                 D.y_res(i,:)=cut(y_res(:,r),pre,round(D.ons(i))-1,post,'padding','nan')';
%                 D.y_raw(i,:)=cut(y_raw(:,r),pre,round(D.ons(i))-1,post,'padding','nan')';
% %                 D.y_adj(i,:)=cut(y_adj(:,r),pre,round(D.ons(i)),post,'padding','nan')';
% %                 D.y_hat(i,:)=cut(y_hat(:,r),pre,round(D.ons(i)),post,'padding','nan')';
% %                 D.y_res(i,:)=cut(y_res(:,r),pre,round(D.ons(i)),post,'padding','nan')';
% %                 D.y_raw(i,:)=cut(y_raw(:,r),pre,round(D.ons(i)),post,'padding','nan')';
%             end
% 
%             % Add the event and region information to tje structure. 
%             len = size(D.event,1);                
%             D.SN        = ones(len,1)*sn;
%             D.region    = ones(len,1)*r;
%             D.name      = repmat({R{r}.name},len,1);
% %            D.hem       = repmat({R{r}.hem},len,1);
%             D.type      = D.event; 
%             T           = addstruct(T,D);
%         end
% 
%         save(fullfile(baseDir,roiDir, sprintf('%s_glm%d_hrf_post%d.mat',subj_id,glm,post)),'T'); 
%         varargout{1} = T;
%     case 'HRF:fit' % finding optimal parameters for hrf
%         % we have a instruction and a movement trial type and we want to
%         % find a set of parameter that works for both. The important thing
%         % is that the duration of underlying neural event is the only thing
%         % that is different between these two trial types.
%         glm = 3;
%         regN = [1:11]; % SMA, PMv, PMd, M1, S1, aSPL, pSPL, DSVC, MT+, VSVC, EAC
%         duration = 1;
%         onsetshift = 0;
%         pre = 5;
%         post = 25;
%         fig = 1;
%         fitMethod = 'library';
% %         vararginoptions(varargin, {'sn', 'glm', 'regN', 'duration', 'onsetshift', 'pre', 'post', 'fig'});
%         vararginoptions(varargin,{'sn','glm','regN','post','subset','fitMethod'});
% %         cwd = pwd;
%         [subj_id, ~] = get_id(sn);
% 
%         % loop over the subject - TODO: group fitting
%         for s = sn
%             % initialize
%             %T = [];
%             %FIT = [];
% 
%             fprintf('fitting a hrf for subject %s\n', subj_id);
% 
%             cd(fullfile(baseDir,sprintf(glmDir,glm),subj_id)); % cd to subject's GLM dir
%             load SPM;
% 
%             % % updating the filenames - since I ran the code on server
%             % if ~strcmp(fullfile(baseDir,sprintf(glmDir, glm), sprintf('S%02d', s)), SPM.swd) % need to rename SPM
%             %     SPM = spmj_move_rawdata(SPM, fullfile(imagingDir, sprintf('S%02d', s)));
%             % endS
% 
%             load(fullfile(baseDir,roiDir,sprintf('%s_Task_regions.glm_%d.mat',subj_id,glm)));
%             %load(fullfile(baseDir, roiDir, sprintf('%s_SSS_regions.mat', sprintf('S%02d', s)))); % Load R
%             %load(fullfile(roiDir, sprintf('%s_Wang_regions.mat', sprintf('s%02d', s))));
%             R = R(regN); % only use the specified regions
% 
%             Data = region_getdata(SPM.xY.VY, R);
% 
%             reg=[];
%             data=[];
%             for i = 1:length(Data)
%                 reg = [reg ones(1,size(Data{i},2))*regN(i)];
%                 data = [data Data{i}];
%             end
%             clear Data
% 
%             Y = spm_filter(SPM.xX.K, SPM.xX.W*data); % filter out low-frequence trends in Y
%             Yres = spm_sp('r', SPM.xX.xKXs, Y);
%             err_before = sum(sum(Yres.^2))/numel(Yres);
% 
% %             for r = 1:length(SPM.nscan)
% %                 for u=1:length(SPM.Sess(r).U)
% %                     SPM.Sess(r).U(u).dur = ones(size(SPM.Sess(r).U(u).dur))*duration;
% %                     SPM.Sess(r).U(u).ons = SPM.Sess(r).U(u).ons+onsetshift;
% %                 end
% %                 SPM.Sess(r).U=spm_get_ons(SPM,r);
% %             end
% 
%             % Fit a common hrf for specified regions
% %             [SPMf, Yhat, Yres] = fit_hrf(SPM, data);  % 'fit',[1,2]'  
%             [SPMf, Yhat, Yres, p_opt] = spmj_fit_hrfparams(SPM,data,fitMethod);  % 'fit',[1,2]'  
%             % Check Error after
%             err_after = sum(sum(Yres.^2))/numel(Yres);
%             F.region = regN;
%             F.bf_before = SPM.xBF.bf;
%             F.bf_after = SPMf.xBF.bf;
%             F.params_before = SPM.xBF.params;
%             F.params_after = p_opt;
% 
%             F.err_before = err_before;
%             F.err_after = err_after;
%             [err_before err_after p_opt]
%             %g = 0;
% %             figure;
% %             set(gcf,'color','w');
% %             for r=1:8
% %                 subplot(2,4,r);
% %                 sss_hrf('HRF:ROI_hrf_get','sn',s,'glm',glm,'post',post);
% %                 sss_hrf('HRF:ROI_hrf_plot','sn',s,'roi',r,'glm',glm,'post',post,'subset',[]);
% %                 sss_hrf('HRF:ROI_hrf_get','sn',s,'glm',glm,'post',post,'bf',SPMf.xBF.bf);
% %                 load(fullfile(baseDir,roiDir,sprintf('%s_glm%d_hrf_post%d.mat',subj_id,glm,post))); % load T 
% %                 T = getrow(T,T.region==r);
% %                 pre = 10; post = 20;
% %                 % Select a specific subset of things to plot 
% %                 % if glm==2
% %                 %     T.type(find(T.type==6))=5; %% BothRep
% %                 %     T.type(find(T.type==4))=3; %% CueOnlyRep
% %                 %     subset  = find(T.type==3 | T.type==5);
% %                 % elseif glm==0
% %                 %     subset = find(T.type==1);
% %                 % end
% %                 hold on;
% %                 traceplot([-pre:post],T.y_hat,'linestyle','--',...
% %                 'split',[T.type],'subset',subset,...
% %                 'linewidth',3,'linecolor','r'); % ,
% %                 drawline([-8 6 16],'dir','vert','linestyle','--');
% %                 drawline([0],'dir','horz','linestyle','--');
% % %                 hold off;
% %                 xlabel('TR');
% %                 ylabel('activation');
% % %                 title(sprintf('ROI: %s',regname{roi}));
% %                 drawline(0);
% %             end
% 
% %             figure;
% %             for r=1:8
% %                 subplot(2,4,r);
% %                 sss_imana('HRF:ROI_hrf_plot','sn',s,'roi',r,'glm',g,'post',20);
% %             end
% %             
% 
%             %FIT = addstruct(FIT,F);
% 
%             % save
% %             save(fullfile(baseDir, roiDir, sprintf('%s_hrf_ROI_timeseries_glm%d.mat', sprintf('S%02d', s), glm)), '-struct', 'T');
%             save(fullfile(baseDir,roiDir,sprintf('%s_hrf_fit_glm%d_post%d',subj_id,glm,post)),'-struct','F');            
%         end
% 
%     case 'HRF:ROI_hrf_plot'                 % Plot extracted time series
%         % s = varargin{1};
%         % roi = varargin{2};
%         subset = [1:8];
%         roi = [1:11];
%         vararginoptions(varargin,{'sn','roi','glm','post','subset'});
%         [subj_id, ~] = get_id(sn);
%         regname = {'SMA','PMv','PMd','M1','S1','SPLa','SPLp','DSVC','MT+','VSVC','EAC'};
%         load(fullfile(baseDir,roiDir,sprintf('%s_glm%d_hrf_post%d.mat',subj_id,glm,post))); % load T
% 
%         pre = 10;
% %        post = 30;
%         % Select a specific subset of things to plot 
% %         cond_name = {'MotorOnly-L','MotorOnly-S','CueOnly-L','CueOnly-S',...
% %                             'BothRep-L','BothRep-S','NonRep-L','NonRep-S','Non-Interest'};
%         if glm==2
% %             T.type(find(T.type==6))=5; %% BothRep
% %             T.type(find(T.type==4))=3; %% CueOnlyRep
%             subset  = find(T.type==3 | T.type==5);
%         elseif glm==0
%             subset = find(T.type==1);
%         end
%        % figure;        
% %         traceplot([-pre:post],T.y_adj,'errorfcn','stderr',...
% %             'split',[T.type],'subset',subset,...
% %             'leg',regname{roi},'leglocation','bestoutside'); % ,
%         for r = roi
%             figure(r);
%             T_ = getrow(T,T.region==r);
%             traceplot([-pre:post],T_.y_adj,'errorfcn','stderr',...
%                 'split',[T_.type],'subset',subset);
%             hold on;
%             traceplot([-pre:post],T_.y_hat,'linestyle','--',...
%                 'split',[T_.type],'subset',subset,'linewidth',3);
%             drawline([-8 6 16],'dir','vert','linestyle','--');
%             drawline([0],'dir','horz','linestyle','--');
%             hold off;
%             xlabel('TR');
%             ylabel('activation');
%             title(sprintf('ROI: %s',regname{r}));
%             drawline(0);
%         end

end

%%  =======================Project-specific Cases==================================

% switch(what)
%     case 'SUIT:isolate_segment'  
%     % Segment cerebellum into grey and white matter
% 
%         sn = subj_id;
% 
%         vararginoptions(varargin, {'sn'});
% 
%         for s = sn
%             fprintf('- Isolate and segment the cerebellum for %s\n', ...
%                 subj_str{s})
%             spm_jobman('initcfg')
% 
%             % Get the file of subjects anatomical
%             anat_subj_dir  = fullfile(anatomical_dir, subj_str{s});
%             anat_name = 'anatomical.nii'
% 
%             % Define suit folder
%             suit_dir = fullfile(baseDir, 'suit/anatomicals',subj_str{s});
%             % Create suit folder if it does not exist
%             if ~exist(suit_dir, 'dir')
%                 mkdir (suit_dir)
%             end
% 
%             % Copy anatomical_raw file to suit folder
%             source = fullfile(anat_subj_dir, anat_name);
%             dest   = fullfile(suit_dir, anat_name);           
%             copyfile(source, dest);
% 
%             % go to subject directory for suit and isolate segment
%             suit_isolate_seg({dest}, 'keeptempfiles', 1);
%         end % s (sn)
% 
%     case 'SUIT:normalise_dartel' % SUIT normalization using dartel
%         % LAUNCH SPM FMRI BEFORE RUNNING!!!!!
%         sn = subj_id; %subjNum
%         vararginoptions(varargin, 'sn');
% 
%         for s = sn
%             suit_subj_dir = fullfile(baseDir, 'suit/anatomicals', subj_str{s});
%             job.subjND.gray       = {fullfile(suit_subj_dir,'c_anatomical_seg1.nii')};
%             job.subjND.white      = {fullfile(suit_subj_dir,'c_anatomical_seg2.nii')};
%             job.subjND.isolation  = {fullfile(suit_subj_dir,'c_anatomical_pcereb.nii')};
%             suit_normalize_dartel(job);
% 
%         end % s (subjects)
% 
%     case 'SUIT:save_dartel_def'    
%         sn = subj_id; %subjNum
%         % Saves the dartel flow field as a deformation file. 
%         for s = sn
%             cd(fullfile(baseDir,'suit/anatomicals', subj_str{s}));
%             anat_name = 'anatomical';
%             suit_save_darteldef(anat_name);
%         end
% 
%     case 'SUIT:reslice'            % Reslice stuff into suit space 
%         % run the case with 'anatomical' to check the suit normalization
%         % make sure that you reslice into 2mm^3 resolution
% 
%         sn   = subj_id;
%         type = 'con';  % 'betas' or 'con' or 'ResMS' or 'cerebellarGrey' or 'anatomical'
%         mask = 'c_anatomical_pcereb'; % 'cereb_prob_corr_grey' or 'cereb_prob_corr' or 'dentate_mask' or 'pcereb'
%         glm  = 1;             % glm number. Used for reslicing betas and contrasts 
% 
%         vararginoptions(varargin, {'sn', 'type', 'mask', 'glm'})
% 
%         for s = sn
%             suit_dir = fullfile(baseDir, 'suit/anatomical',subj_str{s});
%             switch type
%                 case 'anatomical'
%                     subj_dir = suit_dir;
%                     % Get the name of the anatpmical image
%                     files2map = sprintf('%s_T1w_lpi.nii', subj_str{s});
% 
%                     job.subj.resample = {sprintf('%s,1', files2map)};
%                  case 'betas'
%                     glmSubjDir = fullfile(glm_first_dir,sprintf('glm_%d',glm),subj_str{s});
%                     images='resbeta_0';
%                     source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
%                     cd(glmSubjDir);
%                 case 'con'
%                     glmSubjDir = fullfile(glm_first_dir,sprintf('glm_%d',glm),subj_str{s});
%                     images='con_';
%                     source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
%                     cd(glmSubjDir);
%                 case 'spmT'
%                     glmSubjDir = fullfile(glm_first_dir,sprintf('glm_%d',glm),subj_str{s});
%                     images='spmT_';
%                     source=dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
%                     cd(glmSubjDir);
%             end
%             job.subj.affineTr = {fullfile(baseDir,'suit','anatomicals',subj_str{s},'Affine_c_anatomical_seg1.mat')};
%             job.subj.flowfield= {fullfile(baseDir,'suit','anatomicals',subj_str{s},'u_a_c_anatomical_seg1.nii')};
%             job.subj.resample = {source.name};
%             job.subj.mask     = {fullfile(baseDir,'suit','anatomicals',subj_str{s},sprintf('%s.nii',mask))};
%             job.vox           = [1 1 1];
%             % Replace Nans with zeros to avoid big holes in the the data 
%             for i=1:length(source)
%                 V=spm_vol(source(i).name); 
%                 X=spm_read_vols(V); 
%                 X(isnan(X))=0; 
%                 spm_write_vol(V,X); 
%             end; 
%             suit_reslice_dartel(job);
% 
%             source=fullfile(glmSubjDir,'*wd*');
%             destination=fullfile(baseDir,'suit',sprintf('glm_%d',glm),subj_str{s});
%             movefile(source,destination);
% 
%             fprintf('%s have been resliced into suit space \n',type)
%         end





 





