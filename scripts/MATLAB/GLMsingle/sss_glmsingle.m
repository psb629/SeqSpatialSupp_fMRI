function varargout = sss_glmsingle(what,varargin)

if ispc
    cd '\\wsl.localhost/ubuntu-24.04/home/sungbeenpark/github/SeqSpatialSupp_fMRI/scripts/MATLAB'
elseif ismac
    cd '/Users/sungbeenpark/github/SeqSpatialSupp_fMRI/scripts/MATLAB'
end
sss_init;

% Add GLMsingle to the MATLAB path (in case the user has not already done so).
GLMsingle_dir = fullfile(dir_git,'GLMsingle');

addpath(fullfile(GLMsingle_dir,'matlab'));
% addpath(fullfile(GLMsingle_dir,'matlab', 'utilities'));

% if the submodules were installed we try to add their code to the path
addpath(fullfile(GLMsingle_dir, 'matlab', 'fracridge', 'matlab'));
clear GLMsingle_dir;

% check that the dependencies are in the path
tmp = which('fracridge.m');
if isempty(tmp)
  error('fracridge is missing. Please install from: https://github.com/nrdg/fracridge.git')
end

%% ROI
hemi = {'lh','rh'}; % left & right hemi folder names/prefixes
hem = {'L', 'R'}; % hemisphere: 1=LH 2=RH
hname = {'CortexLeft', 'CortexRight'}; % 'CortexLeft', 'CortexRight', 'Cerebellum'

%% argument inputs
sn = [];
glm = [];
% stimdur = 0.1;
stimdur = 2;
vararginoptions(varargin,{'sn','glm','stimdur','map'});

if isempty(sn)
    error('GLM:design -> ''sn'' must be passed to this function.')
end
[subj_id, S_id] = get_id(fullfile(baseDir,'participants.tsv'), sn);
if isempty(glm)
    error('GLM:design -> ''glm'' must be passed to this function.')
end

glmDir = sprintf('glm_%d',glm);
SPM_folder  = fullfile(baseDir,glmDir,subj_id);

dir_glmsingle = fullfile(baseDir,'GLMsingle');

switch(what)
    case 'GLM:estimate'
        % Name of directory to which outputs will be saved
        outputdir = fullfile(dir_glmsingle, glmDir);
        
        % load SPM file
        SPM = load(fullfile(SPM_folder,'SPM.mat'));
        SPM = SPM.SPM;
        TR = SPM.xY.RT;

        mask = niftiread(fullfile(baseDir,'glm_1',subj_id,'mask.nii'));

        % load design matrix
        design = cell(1,length(SPM.Sess));
        for zz=1:length(SPM.Sess)  % for each run

          ncond = length(SPM.Sess(zz).U);    % number of conditions
          nvol = length(SPM.Sess(zz).row);   % number of volumes

          design{zz} = zeros(nvol,ncond);

          for yy=1:length(SPM.Sess(zz).U)    % for each condition
            design{zz}(round(SPM.Sess(zz).U(yy).ons/TR)+1,yy) = 1;  % set all of the onsets
          end
        end

        % load fMRI data
        data = cell(1,length(SPM.Sess));
        fname = unique(struct2table(SPM.xY.VY).fname);
        for zz=1:length(fname)
            tmp = niftiread(fname{zz});
            mask4d = repmat(mask, [1,1,1,size(tmp,4)]);
            tmp = tmp .* cast(mask4d, 'like', tmp);
            data{zz} = tmp;
        end

        % 
        opt = struct('wantmemoryoutputs',[0 0 0 1]);
        [results] = GLMestimatesingletrial(design,data,stimdur,TR,fullfile(outputdir,subj_id),opt);

        % =========== Ali's stuff after glm single fit
        % copy subject glm mask to glmsingle direcotry:
        copyfile(fullfile(baseDir,'glm_1',subj_id,'mask.nii'),fullfile(outputdir,subj_id,'mask.nii'));
        
        % Save betas as nifti:
        % load glmsingle model
        modelname = 'TYPED_FITHRF_GLMDENOISE_RR.mat';
        m = load(fullfile(outputdir, subj_id, modelname));
        
        % get event onsets and sort chronological:
        D = spmj_get_ons_struct(SPM);
        D.ons = D.ons-1;
        % sort based on onsets:
        blocks = unique(D.block)';
        for b = blocks
            % sorting:
            rows = D.block==b;

            ons = D.ons(rows);
            event = D.event(rows);
            eventname = D.eventname(rows);
            num = D.num(rows);

            % sorting based on onset:
            [~, ix] = sort(ons);
            ons = ons(ix);
            event = event(ix);
            eventname = eventname(ix);
            num = num(ix);
            iti = diff(ons);

            % adding to dataframe:
            D.ons(rows) = ons;
            D.event(rows) = event;
            D.eventname(rows) = eventname;
            D.num(rows) = num;
            idx = find(rows);
            D.iti(idx(2:end),1) = iti;
        end
        
        info_base = niftiinfo(fullfile(baseDir, glmDir, subj_id, 'beta_0001.nii'));
        sz = size(m.modelmd);
        
        niftidir = fullfile(outputdir, subj_id);
        for i = 1:sz(4)
            % make nifti:
            info = info_base;
            info.Filename = [];
            info.Filemoddate = [];
            info.Filesize = [];
            descrip = sprintf('glmsingle:beta (%.4d) - Sn(%d) %s', i, D.block(i), D.eventname{i});
            info.Description = descrip;
            info.raw.descrip = descrip;
            nii = m.modelmd(:,:,:,i);

            % save nifti:
            niftiwrite(nii,fullfile(niftidir, sprintf('beta_%.4d.nii',i)), info);
        end

        % make reginfo.tsv:
        reginfo = {};
        reginfo.sn = sn * ones(size(D.block));
        reginfo.run = D.block;
        reginfo.name = D.eventname;
        reginfo.ons = D.ons;
        dsave(fullfile(outputdir, subj_id, 'reginfo.tsv'),reginfo);

        save R2:
        R2 = m.R2;
        info = info_base;
        info.Filename = [];
        info.Filemoddate = [];
        info.Filesize = [];
        descrip = 'glmsingle:R2 percent';
        info.Description = descrip;
        info.raw.descrip = descrip;
        niftiwrite(R2, fullfile(niftidir,'R2.nii'), info);

        % save HRF:
        HRFindex = m.HRFindex;
        info = info_base;
        info.Filename = [];
        info.Filemoddate = [];
        info.Filesize = [];
        info.Datatype = 'double';
        descrip = 'glmsingle:HRF index';
        info.Description = descrip;
        info.raw.descrip = descrip;
        niftiwrite(HRFindex, fullfile(niftidir,'HRFindex.nii'), info);

    case 'ROI:make_cifti.y_series~'
        %% load SPM file
        SPM = load(fullfile(SPM_folder,'SPM.mat'));
        SPM = SPM.SPM;
        
        %% load reginfo
        reginfo = dload(fullfile(dir_glmsingle, glmDir, subj_id, 'reginfo.tsv'));

        %% load glmsingle model
        typeA = load(fullfile(dir_glmsingle,glmDir,subj_id,'TYPEA_ONOFF.mat'));
        typeD = load(fullfile(dir_glmsingle,glmDir,subj_id,'TYPED_FITHRF_GLMDENOISE_RR.mat'));
        designinfo = load(fullfile(dir_glmsingle,glmDir,subj_id,'DESIGNINFO.mat'));
        designSINGLE = designinfo.designSINGLE;
        stimdur = designinfo.stimdur;
        tr = designinfo.tr;

        %% load mask
        mask = niftiread(fullfile(dir_glmsingle, glmDir, subj_id, 'mask.nii'));
        
        %% 피험자의 surface 공간에서 각 ROI의 node (2-D) 정보 
        fname = fullfile(baseDir,roiDir,subj_id,sprintf('%s.Task_regions.mat',subj_id));
        R = load(fname); R = R.R;

        %% get region voxels, "assuming all images have the SAME AFFINE":
        X=[];Y=[];Z=[];
        for r=1:length(R)
            if (~isempty(R{r}))
                % idx coord' = Affine' * real coord
                [x,y,z]=spmj_affine_transform(R{r}.data(:,1),R{r}.data(:,2),R{r}.data(:,3),inv(SPM.xY.VY(1).mat));
                X=[X;x];Y=[Y;y];Z=[Z;z];
                from(r)=size(X,1)+1;
                to(r)=size(X,1);
            else 
                from(r)=size(X,1)+1;
                to(r)=size(X,1);
            end
        end
        % length(to) = the # of ROIs
        % length(X) = the total summation of # of voxels in all ROIs

        %% load y_raw
        SPM = rename_source(SPM);
        data = cell(1,length(SPM.Sess));
        fname = unique(struct2table(SPM.xY.VY).fname);
        for zz=1:length(fname)
            tmp = niftiread(fname{zz});
            mask4d = repmat(mask, [1,1,1,size(tmp,4)]);
            tmp = tmp .* int16(mask4d);
            data{zz} = tmp;
        end

        %% make y_adj
        y_adj = zeros(length(X), 7400);
        for ii=1:length(designSINGLE)
            data_run = data{ii};
            nTR = SPM.nscan(ii);
        
            % --- A. Get the GLMdenoise Principal Components for this run ---
            % 'pcregressors' is a cell array, one cell per run
            noise_pcs = typeD.pcregressors{ii}; 
            
            % The number of PCs used in the final model is stored in 'pcnum'
            num_pcs_to_use = typeD.pcnum;
            
            % Keep only the PCs that were actually used in the model
            if num_pcs_to_use > 0
                selected_pcs = noise_pcs(:, 1:num_pcs_to_use);
            else
                selected_pcs = []; % No PCs were selected
            end
            
            % --- B. Generate the Polynomial Regressors ---
            num_timepoints = nTR;
            run_duration_min = (num_timepoints * tr) / 60;
            poly_degree = round(run_duration_min / 2);
            
            % Create the polynomial regressors (including the constant term)
            polynomial_regressors = zeros(num_timepoints, poly_degree + 1);
            for p = 0:poly_degree
                polynomial_regressors(:, p+1) = (linspace(-1, 1, num_timepoints)').^p;
            end
            
            % --- C. Combine into a single nuisance regressor matrix ---
            % Also add any 'extraregressors' you might have used
            extra_regs = []; % Populate this if you used opt.extraregressors
            nuisance_regressors = [selected_pcs, polynomial_regressors, extra_regs];
            
            volume_size = [size(data_run, 1), size(data_run, 2), size(data_run, 3)];
            linear_indices = sub2ind(volume_size, round(X), round(Y), round(Z));
            
            num_voxels_total = size(data_run, 1) * size(data_run, 2) * size(data_run, 3); % 90*90*32 = 259200
            data_reshaped = reshape(data_run, num_voxels_total, size(data_run, 4));
            
            selected_data = double(data_reshaped(linear_indices, :)); % voxels x time data
            
            betas_nuisance = nuisance_regressors \ selected_data'; 
            y_noise = nuisance_regressors * betas_nuisance;
            % y_adj = [y_adj (selected_data' - y_noise)'];
            y_adj(:,(ii-1)*740+1:ii*740) = (selected_data' - y_noise)';
        end

    case 'GLM:t_contrast~'
        dir_work = fullfile(dir_glmsingle, glmDir, subj_id);
        % Make t-maps:
        mask = niftiread(fullfile(baseDir,'glm_1',subj_id,'mask.nii'));
        % load design:
        designinfo = load(fullfile(dir_work, 'DESIGNINFO.mat'));
        % load reginfo:
        T = dload(fullfile(dir_work, 'reginfo.tsv'));
        
        %% load betas:
        betafiles = dir(fullfile(dir_work, 'beta*.nii'));
        beta = {};
        info = {};
        for i = 1:length(betafiles)
            beta{i,1} = niftiread(fullfile(betafiles(i).folder, betafiles(i).name));
            beta{i,1} = beta{i,1} .* single(mask);
            info{i,1} = niftiinfo(fullfile(betafiles(i).folder, betafiles(i).name));
        end
        
        %% estimate t-maps:
        conditions = unique(T.name);
        xCon = {};
        for c = 1:length(conditions)
            cond = conditions{c};
            
            % contrast
            xcon = zeros(size(designinfo.stimorder,2),1);
            xcon(strcmp(T.name, cond)) = 1;
            xcon = xcon / sum(xcon);
            
            % design matrix Xs
            Xs=0; for i=1:length(designinfo.designSINGLE); Xs=Xs+designinfo.designSINGLE{i}; end;
            
            if isempty(xCon) % empty
                xCon = spm_FcUtil('Set', cond, 'T', 'c', xcon, Xs);
            elseif sum(strcmp(cond, {xCon.name})) > 0 % replace
                idx = find(strcmp(cond, {xCon.name}));
                xCon(idx) = spm_FcUtil('Set', cond, 'T', 'c', xcon, Xs);
            else % new
                xCon(end+1) = spm_FcUtil('Set', cond, 'T', 'c', xcon, Xs);
            end
            % idx = find(strcmp(T.name,cond));
            % select_betas = beta(idx);
            % cat4d = cat(4, select_betas{:});
            % tstats = nanmean(cat4d,4) ./ (std(cat4d,[],4)./sqrt(length(idx)));
            % tstats(isnan(tstats)) = 0;
            % % save tstats:
            % infotmp = info{1};
            % infotmp.Filename = [];
            % infotmp.Filemoddate = [];
            % infotmp.Filesize = [];
            % descrip = sprintf('t-stats:%s',conditions{i});
            % infotmp.Description = descrip;
            % infotmp.raw.descrip = descrip;
            % niftiwrite(tstats,fullfile(niftidir,sprintf('tmap_%s.nii',replace(conditions{i},":", "-"))), infotmp);
        end
        % save(fullfile(dir_work, 'xCon.mat'), 'xCon', '-v7.3');
        
        %% run contrasts
        SPM = spm_contrasts(SPM, 1:length(xCon));
    
    case 'WB:vol2surf'
        % https://github.com/DiedrichsenLab/surfAnalysis.git
        % https://github.com/nno/surfing.git
        dir_work = fullfile(dir_glmsingle,glmDir,wbDir,subj_id);
        if (~exist(dir_work,'dir'))
            mkdir(dir_work);
        end
        dir_glm = fullfile(dir_glmsingle,glmDir,subj_id);

        V = {};
        cols = {};
        switch map
            case 'beta' % beta maps (univariate GLM)
                fnames = dir(fullfile(dir_glm,'beta_*.nii'));
                for f = 1:length(fnames)
                    V{f} = fullfile(fnames(f).folder, fnames(f).name);
                    cols{f} = fnames(f).name;
                end
            case 't' % t-values maps (univariate GLM)
                fnames = dir(fullfile(dir_glm,'tmap_*.nii'));
                for f = 1:length(fnames)
                    V{f} = fullfile(fnames(f).folder, fnames(f).name);
                    cols{f} = fnames(f).name;
                end
        end
        writecell(cols, fullfile(dir_work,sprintf('%s.%s_orders.csv',subj_id,map)), 'Delimiter', 'tab');

        for h = [1 2]
            white = fullfile(baseDir,wbDir,S_id,sprintf('%s.%s.white.32k.surf.gii',S_id,hem{h}));
            pial = fullfile(baseDir,wbDir,S_id,sprintf('%s.%s.pial.32k.surf.gii',S_id,hem{h}));
            C1 = gifti(white);
            C2 = gifti(pial);
          
            output = fullfile(dir_work,sprintf('%s.%s.glm_%d.%s.func.gii',subj_id,hem{h},glm,map));
            G = surf_vol2surf(C1.vertices, C2.vertices, V, 'anatomicalStruct', hname{h}, 'exclude_thres', 0.9, 'faces', C2.faces, 'ignore_zeros', 0);
            G = surf_makeFuncGifti(G.cdata,'anatomicalStruct', hname{h}, 'columnNames', cols);
            save(G, output);

            fprintf('mapped %s %s glm_%d \n', subj_id, hem{h}, glm);
        end

    case 'ROI:make_cifti'
        % 피험자의 volume 공간으로 매핑한 ROI 마스크를 이용하여 피험자의
        % ResMS 데이터를 voxel 단위로 얻고, 그 결과를 cifti로 저장

        % 피험자의 surface 공간에서 각 ROI의 node (2-D) 정보 
        fname = fullfile(baseDir,roiDir,S_id,sprintf('%s.Task_regions.mat',S_id));
        % [R, V] = sss_hrf('ROI:deform','sn',sn,'glm',glm,'LR',LR);
        R = load(fname); R = R.R;

        % 피험자 EPI 의 3-D 정보
        % VolFile = fullfile(baseDir,glmDir,subj_id,'mask.nii');
        VolFile = R{1,1}.image;
        V = spm_vol(VolFile);

        %% load glmsingle model
        % typeA = load(fullfile(dir_glmsingle,glmDir,subj_id,'TYPEA_ONOFF.mat'));
        % typeD = load(fullfile(dir_glmsingle, glmDir, subj_id, 'TYPED_FITHRF_GLMDENOISE_RR.mat'));

        %% data
        fprintf('extrating %s from each ROI for subject %s...\n', map, subj_id);
        switch map
            case {'R2','r2'}
                fnames = dir(fullfile(dir_glmsingle, glmDir, subj_id, 'R2.nii'));
            case {'HRF','hrf'}
                fnames = dir(fullfile(dir_glmsingle, glmDir, subj_id, 'HRFindex.nii'));
            case {'BETA','beta','Beta'}
                fnames = dir(fullfile(dir_glmsingle, glmDir, subj_id, 'beta_*.nii'));
        end
        
        dir_output = fullfile(dir_glmsingle, glmDir, roiDir);
        if ~exist(dir_output,'dir')
            mkdir(dir_output);
        end

        for ii = 1:length(fnames)
            folder = fnames(ii).folder;
            name = fnames(ii).name;
            fname = fullfile(folder, name);

            % volume info
            dataset(ii) = spm_vol(fname);
        end
        % R의 정보를 토대로 추출한 1D data
        D = region_getdata(dataset, R);

        % length(D) = the number of ROIs (Left/Right)
        area = 'OTHER';

        for jj = 1:length(D)
            roi = R{jj}.name;
            hem = R{jj}.hem;

            if strcmpi(map,'beta')
                sss = sprintf('cifti.%s.%s.%s.%s.dtseries.nii', hem, subj_id, roi, upper(map));
                dtype = 'series';
                TR = 5;
            else
                sss = sprintf('cifti.%s.%s.%s.%s.dscalar.nii', hem, subj_id, roi, upper(map));
                dtype = 'scalars';
                TR = 1;
            end

            fname = fullfile(dir_output,sss);
            if ~exist(fname,'file')
                cii = region_make_cifti(R{jj},V,'data',D{jj}','dtype',dtype,'struct',area,'TR',TR);
            end
            cifti_write(cii, fname);

            clear cii
        end
end
