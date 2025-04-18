function varargout = sss_hrf(what,varargin)

if ispc
    workdir='F:\SeqSpatialSupp_fMRI';
elseif isfolder('/Volumes/Diedrichsen_data$/data/<project_name>')
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
roiDir          = 'ROI';
glmDir          = 'glm_%d';

dir_git = 'D:/mobaxterm/sungbeenpark/github';
if exist(dir_git, 'dir') && ~contains(path, dir_git)
    addpath(genpath(dir_git));
end

% Read info from participants .tsv file 
hem = {'L', 'R'};

%% MAIN OPERATION 
switch(what)
    case 'ROI:findall' % use this... %rdm, glm3
        % ROI 11개를 모아놓은 ROI.L.SSS.label.gii 파일 만들기
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
        % ROI{1,12}={}; % V1
        
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
        % R=region_calcregions(R) 함수 사용
        % fmask_data.mat file is necessary, imana_sss('WB:vol2surf_resample','sn',s)? ->
        hemis=1;
        % sn=subj_vec;
        glm=0;  % change this glm=1 for S01-S06
        surf='32';        
        [P,ROI_name]=sss_hrf('ROI:findall');
        % 왜 그룹 마스크? 개인 마스크가 아니고?
        PP = load('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/surfaceWB/group32k/fmask_data.mat');
        P = (PP.mask'~=0).*P;
        vararginoptions(varargin,{'sn','glm'});
        n_roi=length(ROI_name)-1;
        
        % loop over subjects
        [subj_id, S_id] = get_id(sn);

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
% region_savasimage
        save(savename,'R');
        
        fprintf('\nROIs have been defined for %s \n',subj_id);
        clear R
        % end  

    case 'ROI:make_cifti'
        TR = 1;
        % dtype = 'scalars';
        dtype = 'series';
        dnames = 'y_raw';
        vararginoptions(varargin,{'fname_load','fname_save','fname_vol','data'});
        R = load(fname_load);
        V = spm_vol(fname_vol);
        cii = region_make_cifti(R,V,'data',data,'dtype',dtype,'dnames',dnames,'TR',TR);

    case 'HRF:ROI_hrf_get'  % Extract raw and estimated time series from ROIs
        % [y_raw, y_adj, y_hat, y_res, B] = region_getts(SPM,R) 함수 사용
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
    
        [subj_id, ~] = get_id(sn);
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
        R = load(fullfile(baseDir,roiDir,sprintf('%s_Task_regions.glm_3.mat',subj_id)));
        R=R.R;
        
        % extract time series data from GifTi
        [y_raw, y_adj, y_hat, y_res, B] = region_getts(SPM,R);
        
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
        
        save(fullfile(baseDir,roiDir, sprintf('%s_glm%d_hrf_post%d.mat',subj_id,glm,post)),'T'); 
        varargout{1} = T;
    case 'HRF:fit' % finding optimal parameters for hrf
        % we have a instruction and a movement trial type and we want to
        % find a set of parameter that works for both. The important thing
        % is that the duration of underlying neural event is the only thing
        % that is different between these two trial types.
        glm = 3;
        regN = [1:11]; % SMA, PMv, PMd, M1, S1, aSPL, pSPL, DSVC, MT+, VSVC, EAC
        duration = 1;
        onsetshift = 0;
        pre = 5;
        post = 25;
        fig = 1;
        fitMethod = 'library';
%         vararginoptions(varargin, {'sn', 'glm', 'regN', 'duration', 'onsetshift', 'pre', 'post', 'fig'});
        vararginoptions(varargin,{'sn','glm','regN','post','subset','fitMethod'});
%         cwd = pwd;
        [subj_id, ~] = get_id(sn);
        
        % loop over the subject - TODO: group fitting
        for s = sn
            % initialize
            %T = [];
            %FIT = [];
            
            fprintf('fitting a hrf for subject %s\n', subj_id);
            
            cd(fullfile(baseDir,sprintf(glmDir,glm),subj_id)); % cd to subject's GLM dir
            load SPM;
            
            % % updating the filenames - since I ran the code on server
            % if ~strcmp(fullfile(baseDir,sprintf(glmDir, glm), sprintf('S%02d', s)), SPM.swd) % need to rename SPM
            %     SPM = spmj_move_rawdata(SPM, fullfile(imagingDir, sprintf('S%02d', s)));
            % endS
            
            load(fullfile(baseDir,roiDir,sprintf('%s_Task_regions.glm_%d.mat',subj_id,glm)));
            %load(fullfile(baseDir, roiDir, sprintf('%s_SSS_regions.mat', sprintf('S%02d', s)))); % Load R
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
            [SPMf, Yhat, Yres, p_opt] = spmj_fit_hrfparams(SPM,data,fitMethod);  % 'fit',[1,2]'  
            % Check Error after
            err_after = sum(sum(Yres.^2))/numel(Yres);
            F.region = regN;
            F.bf_before = SPM.xBF.bf;
            F.bf_after = SPMf.xBF.bf;
            F.params_before = SPM.xBF.params;
            F.params_after = p_opt;
           
            F.err_before = err_before;
            F.err_after = err_after;
            [err_before err_after p_opt]
            %g = 0;
%             figure;
%             set(gcf,'color','w');
%             for r=1:8
%                 subplot(2,4,r);
%                 sss_hrf('HRF:ROI_hrf_get','sn',s,'glm',glm,'post',post);
%                 sss_hrf('HRF:ROI_hrf_plot','sn',s,'roi',r,'glm',glm,'post',post,'subset',[]);
%                 sss_hrf('HRF:ROI_hrf_get','sn',s,'glm',glm,'post',post,'bf',SPMf.xBF.bf);
%                 load(fullfile(baseDir,roiDir,sprintf('%s_glm%d_hrf_post%d.mat',subj_id,glm,post))); % load T 
%                 T = getrow(T,T.region==r);
%                 pre = 10; post = 20;
%                 % Select a specific subset of things to plot 
%                 % if glm==2
%                 %     T.type(find(T.type==6))=5; %% BothRep
%                 %     T.type(find(T.type==4))=3; %% CueOnlyRep
%                 %     subset  = find(T.type==3 | T.type==5);
%                 % elseif glm==0
%                 %     subset = find(T.type==1);
%                 % end
%                 hold on;
%                 traceplot([-pre:post],T.y_hat,'linestyle','--',...
%                 'split',[T.type],'subset',subset,...
%                 'linewidth',3,'linecolor','r'); % ,
%                 drawline([-8 6 16],'dir','vert','linestyle','--');
%                 drawline([0],'dir','horz','linestyle','--');
% %                 hold off;
%                 xlabel('TR');
%                 ylabel('activation');
% %                 title(sprintf('ROI: %s',regname{roi}));
%                 drawline(0);
%             end
            
%             figure;
%             for r=1:8
%                 subplot(2,4,r);
%                 sss_imana('HRF:ROI_hrf_plot','sn',s,'roi',r,'glm',g,'post',20);
%             end
%             
            
            %FIT = addstruct(FIT,F);
            
            % save
%             save(fullfile(baseDir, roiDir, sprintf('%s_hrf_ROI_timeseries_glm%d.mat', sprintf('S%02d', s), glm)), '-struct', 'T');
            save(fullfile(baseDir,roiDir,sprintf('%s_hrf_fit_glm%d_post%d',subj_id,glm,post)),'-struct','F');            
        end
   
    case 'HRF:ROI_hrf_plot'                 % Plot extracted time series
        % s = varargin{1};
        % roi = varargin{2};
        subset = [1:8];
        roi = [1:11];
        vararginoptions(varargin,{'sn','roi','glm','post','subset'});
        [subj_id, ~] = get_id(sn);
        regname = {'SMA','PMv','PMd','M1','S1','SPLa','SPLp','DSVC','MT+','VSVC','EAC'};
        load(fullfile(baseDir,roiDir,sprintf('%s_glm%d_hrf_post%d.mat',subj_id,glm,post))); % load T
        
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
        for r = roi
            figure(r);
            T_ = getrow(T,T.region==r);
            traceplot([-pre:post],T_.y_adj,'errorfcn','stderr',...
                'split',[T_.type],'subset',subset);
            hold on;
            traceplot([-pre:post],T_.y_hat,'linestyle','--',...
                'split',[T_.type],'subset',subset,'linewidth',3);
            drawline([-8 6 16],'dir','vert','linestyle','--');
            drawline([0],'dir','horz','linestyle','--');
            hold off;
            xlabel('TR');
            ylabel('activation');
            title(sprintf('ROI: %s',regname{r}));
            drawline(0);
        end

end


function [subj_id, S_id] = get_id(sn)
    pinfo = dload('D:/mobaxterm/sungbeenpark/github/diedrichsenlab/SeqSpatialSupp_fMRI/participants.tsv');
    subj_id = char(pinfo.subj_id(pinfo.sn==sn));
    S_id = strrep(subj_id,'R','S');

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





 





