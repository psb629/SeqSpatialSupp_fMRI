function varargout = sss_GLM(what,varargin)

if ispc
    cd '\\wsl.localhost/ubuntu-22.04/home/sungbeenpark/github/SeqSpatialSupp_fMRI/scripts/MATLAB'
elseif ismac
    cd '/Users/sungbeenpark/github/SeqSpatialSupp_fMRI/scripts/MATLAB'
end
sss_init;
%% Scanner
TR = 1; % unit: second
numDummys  = 3;  % dummy images at the start of each run (these are discarded)

%% Adjust this for experiment
nTRs = 410;
% nTRs = [410*ones(1,8)];
% nTRs = [410*ones(1,10) 401 406 410 404 410 410 385]; % For S11
% nTRs = [410*ones(1,16) 385]; % for S09

%% HRF parameters
hrf_params_default = [5 15 1 1 6 0 32];
hrf_params = [];
hrf_cutoff = 128;

map = 'beta';
%% argument inputs
sn = [];
glm = [];
vararginoptions(varargin,{'sn','glm','nTRs','hrf_params','map'});
hrf_params = [hrf_params hrf_params_default(length(hrf_params)+1:end)];
if length(hrf_params)~=length(hrf_params_default)
    error('Wrong hrf_param [%s].',num2str(hrf_params))
end

if isempty(sn)
    error('GLM:design -> ''sn'' must be passed to this function.')
end
[subj_id, S_id] = get_id(fullfile(baseDir,'participants.tsv'), sn);

if isempty(glm)
    error('GLM:design -> ''glm'' must be passed to this function.')
end
glmDir = sprintf('glm_%d',glm);

%% ROI
hemi = {'lh','rh'}; % left & right hemi folder names/prefixes
hem = {'L', 'R'}; % hemisphere: 1=LH 2=RH
hname = {'CortexLeft', 'CortexRight'}; % 'CortexLeft', 'CortexRight', 'Cerebellum'

%% prefix of output files
prefix = 'u';

%% conditions
cues = ["L", "S"]; % Letter / Spatial cue
sequences = [1, 2, 3, 4]; % 1: 32451, 2:35124, 3:13254, 4:14523

%% MAIN OPERATION 
switch(what)
    case 'SPM:resave'
        fprintf('Resaving SPM.mat of %s\n',subj_id);
        dir_work = fullfile(baseDir, glmDir, subj_id);
        tmp = load(fullfile(dir_work,'SPM.mat'));
        delete(fullfile(dir_work,'SPM.mat'));
        SPM = tmp.SPM;
        save(fullfile(dir_work,'SPM.mat'),'SPM','-v7.3');

    case 'GLM:all'
        spm_get_defaults('cmdline', true);  % Suppress GUI prompts, no request for overwirte

        %% Check for and delete existing SPM.mat file
        % spm_file = fullfile(baseDir,glmDir,subj_id,'SPM.mat');
        % if exist(spm_file, 'file')
        %     delete(spm_file);
        % end

        %% Run
        sss_GLM('GLM:design','sn',sn,'glm',glm,'hrf_params',hrf_params);
        sss_GLM('GLM:estimate','sn',sn,'glm',glm);
        % sss_GLM('GLM:t_contrast','sn',sn,'glm',glm);
        % sss_GLM('WB:vol2surf','sn',sn,'glm',glm,'map','beta'); % https://github.com/nno/surfing.git, spm nanmean
        % sss_GLM('WB:vol2surf','sn',sn,'glm',glm,'map','ResMS');
        % sss_GLM('WB:vol2surf','sn',sn,'glm',glm,'map','con');
        % sss_GLM('WB:vol2surf','sn',sn,'glm',glm,'map','t');

        %% HRF tunning
        % onset = 0;
        % list_param = [4 14;5 15;6 16];
        % for i=1:length(list_param)
        %     param = list_param(i,:);
        %     for j=[-1,0,1]
        %         tmp = [param(1), param(2)+j];
        %         for k=[6,3]
        %             hrf = [tmp,1,1,k,onset];
        %             disp(hrf);
        %             sss_GLM('GLM:HRF_tuner','sn',sn,'glm',glm,'hrf_params',hrf);
        %         end
        %     end
        % end

    case 'GLM:get_event'
        w = sprintf('GLM:make_glm_%d',glm);
        events = sss_GLM(w,'sn',sn,'glm',glm);

        % dir_work = fullfile(baseDir,behavDir,sprintf('sub-%s',subj_id));
        % if (~exist(dir_work,'dir'))
        %     mkdir(dir_work);
        % end

        % writetable(events, fullfile(dir_work, sprintf('glm_%d.tsv', glm)), 'FileType', 'text', 'Delimiter', '\t')
        varargout{1} = events;

    case 'GLM:init'
        D = dload(fullfile(baseDir,behavDir,sprintf('sub-%s/behav_info.tsv',subj_id)));
        % runs = 1:8;
        % D = getrow(D, ismember(D.BN, runs));

        events.BN = [];
        events.TN = [];
        events.onset = [];
        events.duration = [];
        events.seq = [];
        events.cue = [];
        events.eventname = [];

        varargout{1} = D;
        varargout{2} = events;

    case 'GLM:make_glm_1'

        [D, events] = sss_GLM('GLM:init','sn',sn,'glm',glm);

        %% GLM1 : Trial State
        % 1: (1,L), 2: (1,S), 3: (2,L), 4: (2,S), 5: (3,L), 6: (3,S), 7: (4,L), 8: (4,S)

        onset_shift = -(numDummys*TR) * 1000;

        for idx_s = 1:length(sequences)
            for idx_c = 1:length(cues)
                % row_idxs = ismember(D.cue,cues(idx_c)) & ismember(D.sequence,sequences(idx_s)) & D.TN<=64;
                row_idxs = ismember(D.cue,cues(idx_c)) & ismember(D.sequence,sequences(idx_s));

                events.BN = [events.BN; D.BN(row_idxs)];
                events.TN = [events.TN; D.TN(row_idxs)];
                events.onset = [events.onset; D.onset(row_idxs)+D.prepTime(row_idxs)+onset_shift];
                events.duration = [events.duration; repmat(2000, [sum(row_idxs), 1])];
                events.seq = [events.seq; D.sequence(row_idxs)];
                events.cue = [events.cue; D.cue(row_idxs)];
                events.eventname = [events.eventname; repmat(sprintf("(%d,%s)",sequences(idx_s),cues(idx_c)), [sum(row_idxs), 1])];
            end
        end
        
        events = struct2table(events);
        events.onset = events.onset .* 0.001;
        events.duration = events.duration .* 0.001;

        varargout{1} = events;
    
    case 'GLM:make_glm_2'

        [D, events] = sss_GLM('GLM:init','sn',sn,'glm',glm);
        
        %% GLM2 : Repetition Suppression
        % |(s,c)|(1,L)|(1,S)|(2,L)|(2,S)|(3,L)|(3,S)|(4,L)|(4,S)| 
        % |-----|-----|-----|-----|-----|-----|-----|-----|-----|
        % |(1,L)|(B00 | S01 | C02 | N03)| C04 | N05 | C06 | N07 |
        % |(1,S)|(    | B09 | N10 | C11)| N12 | C13 | N14 | C15 |
        % |(2,L)|(    |     | B18 | S19)| C20 | N21 | C22 | N23 |
        % |(2,S)|(    |     |     | B27)| N28 | C29 | N30 | C31 |
        % |(3,L)|     |     |     |     |(B36 | S37 | C38 | N39)|
        % |(3,S)|     |     |     |     |(    | B45 | N46 | C47)|
        % |(4,L)|     |     |     |     |(    |     | B54 | S55)|
        % |(4,S)|     |     |     |     |(    |     |     | B63)|
        % B: Both-Rep, S: Seq-Rep, C: Cue-Rep, N: No-Rep, (First Finger)

        onset_shift = -(numDummys*TR) * 1000;

        for t = 1:length(D.TN)
            trial = D.TN(t);
            
            events.BN = [events.BN; D.BN(t)];
            events.TN = [events.TN; trial];
            events.onset = [events.onset; D.onset(t)+D.prepTime(t)+onset_shift];
            events.duration = [events.duration; 2000];
            %% current trial
            seq_f = D.sequence(t);
            cue_f = string(D.cue(t));
            events.seq = [events.seq; seq_f];
            events.cue = [events.cue; cue_f];
            %% Drop the first trial of each block 
            if mod(trial,17) == 1
                events.eventname = [events.eventname; "Rest"];
                continue
            end
            [TS_f, ~] = get_TS(seq_f, cue_f);
            %% Previous trial
            seq_i = D.sequence(t-1);
            cue_i = string(D.cue(t-1));
            [TS_i, ~] = get_TS(seq_i, cue_i);
            %% Transition
            [~, Trans] = get_Trans(TS_i, TS_f);
            events.eventname = [events.eventname; Trans];
        end
        
        events = struct2table(events);
        events.onset = events.onset .* 0.001;
        events.duration = events.duration .* 0.001;

        varargout{1} = events;

    case 'GLM:make_glm_2-ii'

        % varargout{1} = sss_GLM('GLM:make_glm_2','sn',sn,'glm',2);

        [D, events] = sss_GLM('GLM:init','sn',sn,'glm',glm);
        
        %% GLM2-ii : Repetition Suppression 2
        onset_shift = -(numDummys*TR) * 1000;

        for t = 1:length(D.TN)-1
            trial = D.TN(t);

            events.BN = [events.BN; D.BN(t)];
            events.TN = [events.TN; trial];
            events.onset = [events.onset; D.onset(t)+D.prepTime(t)+onset_shift];
            events.duration = [events.duration; 2000];
            %% current trial
            seq_i = D.sequence(t);
            cue_i = string(D.cue(t));
            [TS_i, ~] = get_TS(seq_i, cue_i);
            %% Next trial
            seq_f = D.sequence(t+1);
            cue_f = string(D.cue(t+1));
            events.seq = [events.seq; seq_f];
            events.cue = [events.cue; cue_f];
            %% Drop the last trial of each block 
            if mod(trial,17) == 0
                events.eventname = [events.eventname; "Rest"];
                continue
            end
            [TS_f, ~] = get_TS(seq_f, cue_f);
            %% Transition
            [~, Trans] = get_Trans(TS_i, TS_f);
            events.eventname = [events.eventname; Trans];
        end
        
        events = struct2table(events);
        events.onset = events.onset .* 0.001;
        events.duration = events.duration .* 0.001;

        varargout{1} = events;

    case 'GLM:make_glm_3'

        [D, events2] = sss_GLM('GLM:init','sn',sn,'glm',glm);

        %% GLM3 : GLM1 + Preparation
        % 1: (1,L), 2: (1,S), 3: (2,L), 4: (2,S), 5: (3,L), 6: (3,S), 7: (4,L), 8: (4,S)
        events1 = sss_GLM('GLM:make_glm_1','sn',sn,'glm',1);

        % 1+8: p(1,L), 2+8: p(1,S), 3+8: p(2,L), 4+8: p(2,S), 5+8: p(3,L), 6+8: p(3,S), 7+8: p(4,L), 8+8: p(4,S)

        onset_shift = -(numDummys*TR) * 1000;

        for idx_s = 1:length(sequences)
            % preparation
            for idx_c = 1:length(cues)
                row_idxs = ismember(D.cue,cues(idx_c)) & ismember(D.sequence,sequences(idx_s));

                events2.BN = [events2.BN; D.BN(row_idxs)];
                events2.TN = [events2.TN; D.TN(row_idxs)];
                events2.onset = [events2.onset; D.onset(row_idxs)+onset_shift];
                events2.duration = [events2.duration; repmat(1000, [sum(row_idxs), 1])];
                events2.seq = [events2.seq; D.sequence(row_idxs)];
                events2.cue = [events2.cue; D.cue(row_idxs)];
                events2.eventname = [events2.eventname; repmat(sprintf("p(%d,%s)",sequences(idx_s),cues(idx_c)), [sum(row_idxs), 1])];
            end
        end
        
        events2 = struct2table(events2);
        events2.onset = events2.onset .* 0.001;
        events2.duration = events2.duration .* 0.001;

        varargout{1} = [events1; events2];

    case 'GLM:design'
        %% dependency:
        % https://github.com/spm/spm.git
        % https://github.com/jdiedrichsen/rwls.git

        %% import globals from spm_defaults
        % spm_get_defaults;
        global defaults; 
        if (isempty(defaults))
            spm_defaults;
        end

        %% load event file
        % events_file = fullfile(baseDir,behavDir,sprintf('sub-%s/glm_%d.tsv',subj_id,glm));
        Dd = sss_GLM('GLM:get_event','sn',sn,'glm',glm);

        regressors = unique(Dd.eventname);
        nRegr = length(regressors);

        %% cvi type
        % cvi_type = 'wls';
        cvi_type = 'fast';        
        
        %% initiate J
        J = [];
        J.dir = {fullfile(baseDir,glmDir,subj_id)};
        if ~exist(J.dir{1},"dir")
            mkdir(J.dir{1});
        end

        J.timing.units = 'secs'; % timing unit that all timing in model will be
        J.timing.RT = 1; % TR (in seconds, as per 'J.timing.units')
        % number of temporal bins in which the TR is divided,
        % defines the discrtization of the HRF inside each TR
        J.timing.fmri_t = 16;
        % slice number that corresponds to that acquired halfway in
        % each TR
        J.timing.fmri_t0 = 1;

        T = [];

        epi_files = dir(fullfile(baseDir,imagingDir,subj_id,sprintf('%s_run_*.nii',subj_id)));

        runs = 1:8;
        itaskUni = 0;
        for run=runs
            % Setup scans for current session
            N = {};
            cnt = 1;
            for ii = numDummys+1:nTRs
                N{cnt} = fullfile(epi_files(run).folder, sprintf('%s,%d',epi_files(run).name,ii));
                cnt = cnt + 1;
            end
            % N = {fullfile(epi_files(run).folder, epi_files(run).name)};
            J.sess(run).scans = N;

            % Preallocate memory for conditions
            J.sess(run).cond = repmat(struct('name', '', 'onset', [], 'duration', []), nRegr, 1);

            for regr = 1:nRegr
                itaskUni = itaskUni + 1;
                rows = find(Dd.BN == run & strcmp(Dd.eventname, regressors(regr)));

                % Regressor name
                J.sess(run).cond(regr).name = regressors{regr};

                % Define onset
                J.sess(run).cond(regr).onset  = Dd.onset(rows);
                
                % Define durationDuration(regr));
                J.sess(run).cond(regr).duration = Dd.duration(rows); % needs to be in seconds

                % Define time modulator
                % Add a regressor that account for modulation of betas over time
                J.sess(run).cond(regr).tmod=0;

                % Orthogonalize parametric modulator
                % Make the parametric modulator orthogonal to the main regressor
                J.sess(run).cond(regr).orth=0;
                
                % Define parametric modulators
                % Add a parametric modulators, like force or reaction time. 
                J.sess(run).cond(regr).pmod=struct('name',{},'param',{},'poly',{});

                % add the condition info to the reginfo structure

                % filling in "reginfo"
                TT.sn       = sn;
                TT.run      = run;
                TT.cond     = regr;
                TT.regIdx   = itaskUni;
                % TT.seq      = unique(Dd.seq(rows));
                % TT.cue      = unique(Dd.cue(rows));
                TT.reg      = cellstr(regressors(regr));
                % TT.n_rep     = sum(idx);
        
                T = addstruct(T, TT);
            end

            %% J.sess(run).multi
            % Purpose: Specifies multiple conditions for a session. Usage: It is used
            % to point to a file (.mat or .txt) that contains multiple conditions,
            % their onsets, durations, and names in a structured format. If you have a
            % complex design where specifying conditions manually within the script is
            % cumbersome, you can prepare this information in advance and just
            % reference the file here. Example Setting: J.sess(run).multi =
            % {'path/to/multiple_conditions_file.mat'}; If set to {' '}, it indicates
            % that you are not using an external file to specify multiple conditions,
            % and you will define conditions directly in the script (as seen with
            % J.sess(run).cond).
            J.sess(run).multi={''};

            %% J.sess(run).regress
            % Purpose: Allows you to specify additional regressors that are not
            % explicitly modeled as part of the experimental design but may account for
            % observed variations in the BOLD signal. Usage: This could include
            % physiological measurements (like heart rate or respiration) or other
            % variables of interest. Each regressor has a name and a vector of values
            % corresponding to each scan/time point.
            J.sess(run).regress=struct('name',{},'val',{});

            %% J.sess(run).multi_reg
            % Purpose: Specifies a file containing multiple
            % regressors that will be included in the model as covariates. Usage: This
            % is often used for motion correction, where the motion parameters
            % estimated during preprocessing are included as regressors to account for
            % motion-related artifacts in the BOLD signal. Example Setting:
            % J.sess(run).multi_reg = {'path/to/motion_parameters.txt'}; The file
            % should contain a matrix with as many columns as there are regressors and
            % as many rows as there are scans/time points. Each column represents a
            % different regressor (e.g., the six motion parameters from realignment),
            % and each row corresponds to the value of those regressors at each scan.
            J.sess(run).multi_reg={''};
            
            % Define high pass filter cutoff (in seconds): see glm cases.
            J.sess(run).hpf=hrf_cutoff;
        end

        % Specify factorial design
        J.fact = struct('name', {}, 'levels', {});

        % Specify hrf parameters for convolution with
        % regressors
        J.bases.hrf.derivs = [0 0];

        J.bases.hrf.params = hrf_params; % positive and negative peak of HRF - set to [] if running wls (?)
        defaults.stats.fmri.hrf=J.bases.hrf.params; 
        
        % Specify the order of the Volterra series expansion 
        % for modeling nonlinear interactions in the BOLD response
        % *Example Usage*: Most analyses use 1, assuming a linear
        % relationship between neural activity and the BOLD
        % signal.
        J.volt = 1;

        % Specifies the method for global normalization, which
        % is a step to account for global differences in signal
        % intensity across the entire brain or between scans.
        J.global = 'None';

        % remove voxels involving non-neural tissue (e.g., skull)
        J.mask = {fullfile(baseDir,anatomicalDir, S_id, 'rmask_noskull.nii,1')};
        
        % Set threshold for brightness threshold for masking 
        % If supplying explicit mask, set to 0  (default is 0.8)
        J.mthresh = 0.05;
        
        % Create map where non-sphericity correction must be applied
        J.cvi_mask = {fullfile(baseDir, anatomicalDir, S_id, 'rmask_gray.nii')};
        
        % Method for non sphericity correction
        J.cvi = cvi_type;
        
        % remove empty rows (e.g., when skipping runs)
        J.sess = J.sess(~arrayfun(@(x) all(structfun(@isempty, x)), J.sess));

        % save the variable T as reginfo.tsv
        dsave(fullfile(J.dir{1},'reginfo.tsv'), T);

        % run
        spm_rwls_run_fmri_spec(J);

        % Save the GLM file for this subject.
        tmp = load(fullfile(J.dir{1},'SPM.mat'));
        delete(fullfile(J.dir{1},'SPM.mat'));
        SPM = tmp.SPM;
        save(fullfile(J.dir{1},'SPM.mat'),'SPM','-v7.3');

        fprintf('- estimates for glm_%d subject %s has been saved for %s \n', glm, subj_id);
    
    case 'GLM:estimate' % estimate beta coefficient
        % Estimate the GLM from the appropriate SPM.mat file.
        % Make GLM files with case 'GLM_make'.
        dir_work = fullfile(baseDir, glmDir, subj_id);
        load(fullfile(dir_work, 'SPM.mat'));
        SPM.swd = dir_work;

        %% Non-interest index
        iB = SPM.xX.iB;
        save(fullfile(dir_work, "iB.mat"), "iB", '-v7.3');

        %% Run the GLM.
        spm_rwls_spm(SPM);
        % if fig==1
        %     dm = SPM.xX.xKXs.X(SPM.Sess(1).row,SPM.Sess(1).col); % design matrix for one run
        %     figure; imagesc(dm); axis square; colorbar;
        %     title('Design matrix'); xlabel('Regressors'); ylabel('Volumes'); set(gca, 'fontsize', 13)
        %     cv = cov(dm);                   % covariance matrix for one run
        %     varE = nanmean(diag(inv(cv)));  % mean variance of the regression estimates
        %     varX = nanmean(1./diag(cv));    % mean variance of the regressors in the design matrix
        %     vif = varE ./ varX;             % variance inflation factor
        %     figure;
        %     subplot(2,2,1)
        %     imagesc(cv); axis square; title('Covariance matrix'); colorbar;
        %     subplot(2,2,2)
        %     imagesc(inv(cv)); axis square; title('Inverse of the covariance'); colorbar;
        %     subplot(2,2,3)
        %     imagesc(1./diag(cv)); axis square; title('Mean variance per regressor'); colorbar;
        %     subplot(2,2,4)
        %     imagesc(inv(cv)./(1./cv)); axis square; title('Covariance inflation factor'); colorbar;
        %     fprintf(1, '\nVariance estimates: %2.3f\nVariance regressors: %2.3f\nVariance inflation factor: %2.3f (the closer to 1, the better)\n', varE, varX, vif);
        % end

        %% resave SPM
        tmp = load(fullfile(dir_work,'SPM.mat'));
        delete(fullfile(dir_work,'SPM.mat'));
        SPM = tmp.SPM;
        save(fullfile(dir_work,'SPM.mat'),'SPM','-v7.3');

        %% Save the prewhitened design matrix
        nKX = SPM.xX.nKX;
        save(fullfile(dir_work,'nKX_data.mat'),'nKX');

    case 'GLM:t_contrast'
        % get the subject id folder name
        fprintf('Contrasts for participant %s\n', subj_id);
        dir_work = fullfile(baseDir, glmDir, subj_id);

        SPM = load(fullfile(dir_work,'SPM.mat')); SPM=SPM.SPM;

        T = dload(fullfile(dir_work,'reginfo.tsv'));
        T.reg = cellstr(string(T.reg));
        contrasts = unique(T.reg);

        for c = 1:length(contrasts)
            contrast_name = contrasts{c};
            xcon = zeros(size(SPM.xX.X,2), 1);
            xcon(strcmp(T.reg, contrast_name)) = 1;
            xcon = xcon / sum(xcon);
            if ~isfield(SPM, 'xCon') | isempty(SPM.xCon)
                SPM.xCon = spm_FcUtil('Set', contrast_name, 'T', 'c', xcon, SPM.xX.xKXs);
                cname_idx = 1;
            elseif sum(strcmp(contrast_name, {SPM.xCon.name})) > 0
                idx = find(strcmp(contrast_name, {SPM.xCon.name}));
                SPM.xCon(idx) = spm_FcUtil('Set', contrast_name, 'T', 'c', xcon, SPM.xX.xKXs);
                cname_idx = idx;
            else
                SPM.xCon(end+1) = spm_FcUtil('Set', contrast_name, 'T', 'c', xcon, SPM.xX.xKXs);
                cname_idx = length(SPM.xCon);
            end
        end
        save(fullfile(dir_work,'SPM.mat'),'SPM','-v7.3');

        % extra contrasts
        SPM = sss_GLM('GLM:custom_contrast','sn',sn,'glm',glm);

        % run contrasts
        SPM = spm_contrasts(SPM,1:length(SPM.xCon));
            
        % rename contrast images and spmT images
        for idx = 1:length(SPM.xCon)
            % con
            oldName = fullfile(dir_work, sprintf('%s', SPM.xCon(idx).Vcon.fname));
            newName = fullfile(dir_work, sprintf('con_%s.nii', SPM.xCon(idx).name));
            movefile(oldName, newName);
            % t
            oldName = fullfile(dir_work, sprintf('%s', SPM.xCon(idx).Vspm.fname));
            newName = fullfile(dir_work, sprintf('spmT_%s.nii', SPM.xCon(idx).name));
            movefile(oldName, newName);
        end

        % SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
        % save(fullfile(dir_work,'SPM_light.mat'),'SPM')
        save(fullfile(dir_work,'SPM.mat'),'SPM','-v7.3');

    case 'GLM:custom_contrast'
        fprintf('(Additional) Contrasts for participant %s\n', subj_id);
        dir_work = fullfile(baseDir, glmDir, subj_id);

        SPM = load(fullfile(dir_work,'SPM.mat')); SPM=SPM.SPM;

        T = dload(fullfile(dir_work,'reginfo.tsv'));
        T.reg = cellstr(string(T.reg));

        nRun = length(SPM.Sess);
        nCond = length(SPM.Sess(1).Fc);
        switch glm
        case 1
            %% (1,L), (1,S), (2,L), (2,S), (3,L), (3,S), (4,L), (4,S)
            contrasts = {'Letter','Spatial','Letter-Spatial'};
            xcons = [1 0 1 0 1 0 1 0 ; 0 1 0 1 0 1 0 1 ; 1 -1 1 -1 1 -1 1 -1];
        case {2,"2-ii"}
            %% B_L, B_S, C_L, C_S, N_L, N_S, Rest, S_L, S_S
            contrasts = {
                'Letter','Spatial','Letter-Spatial', ...
                'wRS_L','wRS_S','wRS_L-S',...
                'acRS_L','acRS_S','acRS_L-S',...
            };
            xcons = [
                1 0 1 0 1 0 0 1 0 ; 0 1 0 1 0 1 0 0 1 ; 1 -1 1 -1 1 -1 0 1 -1 ;
                1 0 -1 0 0 0 0 0 0 ; 0 1 0 -1 0 0 0 0 0 ; 1 -1 -1 1 0 0 0 0 0 ;
                0 0 0 0 -1 0 0 1 0 ; 0 0 0 0 0 -1 0 0 1 ; 0 0 0 0 -1 1 0 1 -1 ;
            ];
        case 3
            %% (1,L), (1,S), (2,L), (2,S), (3,L), (3,S), (4,L), (4,S)
            contrasts = {'Letter','Spatial','Letter-Spatial'};
            xcons = [
                1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0; 
                0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0;
                1 -1 1 -1 1 -1 1 -1 0 0 0 0 0 0 0 0
            ];
        end

        for c = 1:length(contrasts)
            contrast_name = contrasts{c};
            xcon = xcons(c,:);
            cnt = sum(xcon > 0);
            xcon = [repmat(xcon,1,nRun) zeros(1,nRun)]'/(cnt*nRun);
            
            idx = nCond+c;
            SPM.xCon(idx) = spm_FcUtil('Set', contrast_name, 'T', 'c', xcon, SPM.xX.xKXs);
            % cname_idx = idx;
            % SPM = spm_contrasts(SPM,cname_idx);
        end

        varargout{1} = SPM;

    case 'WB:vol2surf' % map indiv vol contrasts (.nii) onto surface (.gii)
        dir_work = fullfile(baseDir,wbDir,glmDir,S_id);
        if (~exist(dir_work,'dir'))
            mkdir(dir_work);
        end
        dir_glm = fullfile(baseDir,glmDir,subj_id);

        V = {};
        cols = {};
        switch map
            case 'beta' % beta maps (univariate GLM)
                load(fullfile(dir_glm,'SPM.mat'));
                fnames = dir(fullfile(dir_glm,'beta_*.nii'));
                fnames = fnames(SPM.xX.iC);
                for f = 1:length(fnames)
                    V{f} = fullfile(fnames(f).folder, fnames(f).name);
                    cols{f} = fnames(f).name;
                end
            case 'ResMS' % residual
                V{1} = fullfile(dir_glm,'ResMS.nii');
                cols{1} = 'ResMS.nii';
            case 't' % t-values maps (univariate GLM)
                fnames = dir(fullfile(dir_glm,'spmT_*.nii'));
                for f = 1:length(fnames)
                    V{f} = fullfile(fnames(f).folder, fnames(f).name);
                    cols{f} = fnames(f).name;
                end
            case 'con' % contrast beta maps (univariate GLM)
                fnames = dir(fullfile(dir_glm,'con_*.nii'));
                for f = 1:length(fnames)
                    V{f} = fullfile(fnames(f).folder, fnames(f).name);
                    cols{f} = fnames(f).name;
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
        end
        writecell(cols, fullfile(dir_work,sprintf('%s.%s_orders.csv',subj_id,map)), 'Delimiter', 'tab');

        for h = [1 2]
            white = fullfile(baseDir,wbDir,S_id,sprintf('%s.%s.white.32k.surf.gii',S_id,hem{h}));
            pial = fullfile(baseDir,wbDir,S_id,sprintf('%s.%s.pial.32k.surf.gii',S_id,hem{h}));
            C1 = gifti(white);
            C2 = gifti(pial);
          
            output = fullfile(dir_work,sprintf('%s.%s.glm_%d.%s.func.gii',subj_id,hem{h},glm,map));
            % https://github.com/nno/surfing.git
            G = surf_vol2surf(C1.vertices, C2.vertices, V, 'anatomicalStruct', hname{h}, 'exclude_thres', 0.9, 'faces', C2.faces, 'ignore_zeros', 0);
            G = surf_makeFuncGifti(G.cdata,'anatomicalStruct', hname{h}, 'columnNames', cols);
            save(G, output);

            fprintf('mapped %s %s glm_%d \n', subj_id, hem{h}, glm);
        end

    case 'GLM:HRF_tuner'
        %% HRF tunning
        % params = mat2str(hrf_params(1:6));
        fprintf('GLM:HRF_tuner - %s: %s\n',subj_id,mat2str(hrf_params));
        dir_output = fullfile(baseDir,glmDir,subj_id,'hrf_tune');
        if (~exist(dir_output,'dir'))
            mkdir(dir_output);
        end
        % library = dload(fullfile(dir_SSS,'getcanonicalhrflibrary.tsv'));
        % library = readtable(fullfile(dir_SSS,'getcanonicalhrflibrary.tsv'), 'FileType', 'text', 'Delimiter', '\t');

        %% load SPM.mat (GLM information)
        SPM = load(fullfile(baseDir,glmDir,subj_id,'SPM.mat'));
        SPM = SPM.SPM;

        %% Get the hemodynamic response in micro-time resolution
        SPM.xBF.UNITS    = 'secs'; % units of the hrf
        SPM.xBF.T        = 16; % microtime resolution: number of time samples per scan
        SPM.xBF.dt       = SPM.xY.RT/SPM.xBF.T;  % Delta-t per sample 
        SPM.xBF.T0       = 1; % first time bin
        SPM.xBF.name     = 'fitted_hrf';
        SPM.xBF.order    = 1;
        SPM.xBF.Volterra = 1;  % volterra expansion order?
        SPM.xBF.params = hrf_params;
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

        %% 피험자의 surface 공간에서 각 ROI의 node (2-D) 정보 
        fname = fullfile(baseDir,roiDir,S_id,sprintf('%s.Task_regions.mat',S_id));
        R = load(fname); R = R.R;

        %% 피험자 EPI 의 3-D 정보
        % VolFile = R{1,1}.image;
        VolFile = fullfile(baseDir,glmDir,subj_id,'mask.nii');
        V = spm_vol(VolFile);

        for i=1:length(R)
            roi = R{i}.name;
            area = 'OTHER';

            hemisphere = R{i}.hem;
            % if hemisphere=='L'
            %     area = [roi '_LEFT'];
            % elseif hemisphere=='R'
            %     area = [roi '_RIGHT'];
            % end

            %% load y_raw
            dir_yraw = fullfile(baseDir,roiDir,subj_id);
            fname = fullfile(dir_yraw, sprintf('cifti.%s.%s.%s.y_raw.dtseries.nii',hemisphere,subj_id,roi));
            cii = cifti_read(fname);
            Yraw = double(cii.cdata(:,:)');
            tmp = numDummys+1:nTRs;
            idx = [];
            for j=0:7
                idx = [idx tmp+nTRs*j];
            end
            Yraw = Yraw(idx,:);
            clear cii
            
            % Restimate the betas and get predicted and residual response
            [beta, Yhat, Yres] = spmj_glm_fit(SPM,Yraw);

            params = sprintf('[%d',hrf_params(1));
            for ii=2:length(hrf_params)
                params = [params sprintf(',%d',hrf_params(ii))];
            end
            params = [params ']'];

            %% y_hat
            cii = region_make_cifti(R{i},V,'data',Yhat','dtype','series','struct',area,'TR',TR);
            fname = fullfile(dir_output, sprintf('cifti.%s.%s.%s.%s.%s.y_hat.dtseries.nii',hemisphere,glmDir,params,subj_id,roi));
            cifti_write(cii, fname);
            clear cii

            %% y_res
            cii = region_make_cifti(R{i},V,'data',Yres','dtype','series','struct',area,'TR',TR);
            fname = fullfile(dir_output, sprintf('cifti.%s.%s.%s.%s.%s.y_res.dtseries.nii',hemisphere,glmDir,params,subj_id,roi));
            cifti_write(cii, fname);
            clear cii

            %% beta
            cii = region_make_cifti(R{i},V,'data',beta(SPM.xX.iC,:)','dtype','scalars','struct',area,'TR',TR);
            fname = fullfile(dir_output, sprintf('cifti.%s.%s.%s.%s.%s.beta.dscalar.nii',hemisphere,glmDir,params,subj_id,roi));
            cifti_write(cii, fname);
            clear cii
        end
        xBF = SPM.xBF;
        save(fullfile(dir_output,sprintf('xBF_%s.mat',params)),'xBF','-v7.3');

        %% save vector information
        % % column head
        % tmp = {};
        % for run=1:length(SPM.Sess)
        %     tmp = cat(2,tmp,SPM.Sess(run).U(:).name);
        % end
        % vectors = tmp;
        % 
        % for run=1:length(SPM.Sess)
        %     cols = SPM.Sess(run).col;
        %     % partition vector
        %     vectors(2,cols) = {run};
        % 
        %     for col=cols
        %         idx = mod(col,length(cols));
        %         if idx==0
        %             idx = idx+length(cols);
        %         end
        %         % condition vector
        %         tmp = cellfun(@(x) sscanf(x,'Trial-State %d'), SPM.Sess(run).U(idx).name, 'UniformOutput', false);
        %         if ~isempty(tmp{1})
        %             vectors(3,col) = tmp;
        %         else
        %             vectors(3,col) = {0};
        %         end
        %     end
        % end
        % save(fullfile(dir_output,sprintf('vec_%s.mat',params)),'vectors','-v7.3');

    case 'GLM:psc' % calculate percent signal change for selected contrasts

        % Go to subject's directory and load SPM info
        cd(fullfile(baseDir,sprintf('glm_%d', glm), subj_id));
        if glm~=5
            load SPM;
        else
            load(fullfile(baseDir,sprintf('glm_%d', 1), subj_id,'SPM.mat'));
        end
%        load('SPM_info.mat');
        %load('SPM_info.mat');

        X = (SPM.xX.X(:,SPM.xX.iC));      % Design matrix - raw
        h = median(max(X));               % Height of response;
        numB = length(SPM.xX.iB);         % Partitions - runs
        P = cell(1,numB+1);
        run = 0;
        for p=SPM.xX.iB
            run=run+1;
            P{run}=sprintf('beta_%4.4d.nii',p);       % get the intercepts and use them to calculate the baseline (mean images)
        end

        for con=1:length(SPM.xCon)    % all contrasts
            P{numB+1}=sprintf('con_%s.nii',SPM.xCon(con).name);           % replace with T.con_name{con} in case
            outname=sprintf('psc_%s.nii',SPM.xCon(con).name);
            formula = sprintf(['100.*%f.*' pscFormula(numB)],h);
            spm_imcalc_ui(P,outname,formula,{0,[],spm_type(16),[]});       % Calculate percent signal change
            fprintf('Contrast: %s\n',SPM.xCon(con).name);
        end

        fprintf('Subject %s - Done\n', subj_id);

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
            cd(fullfile(baseDir,sprintf('glm_%d',glm),sprintf('S%02d',s)))
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
            for run=1:11
%             for r=1:length(R)
            %             for r=1:length(R) % for each region
                
                Y=region_getdata(V,R{run});  
                S.psc=nanmean(Y,2);  % use nanmean, SKim
                S.nanperc = (length(find(isnan(Y)==1))/prod(size(Y)))*ones(length(O),1);
                S.hemi=repmat(hemi,length(O),1);
                S.roi=repmat(run,length(O),1);
                S.SN=repmat(s,length(O),1);
                S.SeqType= SeqType;
                S.RepType = RepType;
                S.cond = [1:8]';

                T=addstruct(T,S);
                fprintf('%d.',run)
            end
        end
        % save T
%         savename=fullfile(baseDir,roiDir,sprintf('psc_glm%d_%s_N=%d.mat',glm,'SSS',numel(sn)));

        savename=fullfile(baseDir,roiDir,sprintf('psc_glm%d_%s_N=%d.mat',glm,'Task',numel(sn)));
        save(savename,'-struct','T');
        fprintf('\n')

   case 'ROI_calc_rdm'         % ROI MVA 1: Extracts raw timeseries and get prewhitened beta %rdm, glm3
        % sss_imana('ROI_get_prewhitened_beta','sn',s,'glm',glm) ->
        glm = 3;
        roi = [1:11];
%         normmode = 'runwise';
%         normmethod = 'multivariate';
        vararginoptions(varargin, {'sn','glm','roi'});
        fname = fullfile(dir_git,'diedrichsenlab/SeqSpatialSupp_fMRI/participants.tsv');
        [subj_id, S_id] = get_id(fname, sn);
        load(fullfile(baseDir, sprintf('glm_%d',glm),subj_id,'SPM.mat'));
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
        for run=roi
            RDM  = rsa.distanceLDC(Data.beta{run},partvec,condvec,SPM.xX.xKXs.X); 
%             RDM = rsa.spm.distanceLDCraw(data{r},SPM,condvec,'normmode','runwise','normmethod','multivariate');
            S.glm = glm;
            S.subj = sn;
            S.numVox = size(data{run},2);
            S.name = {R{run}.name};
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
