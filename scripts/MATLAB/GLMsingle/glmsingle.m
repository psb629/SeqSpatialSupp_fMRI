%% Start fresh
clear
clc
close all

if ispc
    cd '\\wsl.localhost/ubuntu-22.04/home/sungbeenpark/github/SeqSpatialSupp_fMRI/scripts/MATLAB'
elseif ismac
    cd '/Users/sungbeenpark/github/SeqSpatialSupp_fMRI/scripts/MATLAB'
end
sss_init;

% Add GLMsingle to the MATLAB path (in case the user has not already done so).
GLMsingle_dir = fullfile(dir_git,'GLMsingle');

addpath(fullfile(GLMsingle_dir, 'matlab'));
addpath(fullfile(GLMsingle_dir, 'matlab', 'utilities'));

% if the submodules were installed we try to add their code to the path
addpath(fullfile(GLMsingle_dir, 'matlab', 'fracridge', 'matlab'));

% check that the dependencies are in the path
tmp = which('fracridge.m');
if isempty(tmp)
  error('fracridge is missing. Please install from: https://github.com/nrdg/fracridge.git')
end

clear GLMsingle_dir;

sn = 1;
[subj_id, S_id] = get_id(fullfile(baseDir,'participants.tsv'), sn);

glm = 1;
glmDir = sprintf('glm_%d',glm);

SPM_folder  = fullfile(baseDir,glmDir,subj_id);
% stimdur     = 0.1;
stimdur     = 2;

% Name of directory to which outputs will be saved
outputdir   = fullfile(baseDir,'GLMsingle',glmDir);

%
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
copyfile(fullfile(baseDir,'glm1',subj_id,'mask.nii'),fullfile(outputdir,subj_id,'mask.nii'));


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

info_base = niftiinfo(fullfile(baseDir,'glm1',subj_id,'beta_0001.nii'));
sz = size(m.modelmd);

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
    niftidir = fullfile(outputdir,subj_id);
    niftiwrite(nii,fullfile(niftidir,sprintf('beta_%.4d.nii',i)), info);
end

% make reginfo.tsv:
reginfo = {};
reginfo.sn = sn * ones(size(D.block));
reginfo.run = D.block;
reginfo.name = D.eventname;
reginfo.ons = D.ons;
dsave(fullfile(outputdir, subj_id, 'reginfo.tsv'),reginfo);

% save R2:
R2 = m.R2;
info = info_base;
info.Filename = [];
info.Filemoddate = [];
info.Filesize = [];
descrip = 'glmsingle:R2 percent';
info.Description = descrip;
info.raw.descrip = descrip;
niftidir = fullfile(outputdir,subj_id);
niftiwrite(R2,fullfile(niftidir,'R2.nii'), info);


% Make t-maps:
% load reginfo:
reginfo = dload(fullfile(outputdir, subj_id, 'reginfo.tsv'));

% load betas:
betafiles = dir(fullfile(outputdir, subj_id, 'beta*.nii'));
beta = {};
info = {};
for i = 1:length(betafiles)
    beta{i,1} = niftiread(fullfile(betafiles(i).folder, betafiles(i).name));
    beta{i,1} = beta{i,1} .* single(mask);
    info{i,1} = niftiinfo(fullfile(betafiles(i).folder, betafiles(i).name));
end

% estimate t-maps:
conditions = unique(reginfo.name);
for i = 1:length(conditions)
    c = conditions{i};
    idx = find(strcmp(reginfo.name,c));
    select_betas = beta(idx);
    cat4d = cat(4, select_betas{:});
    tstats = nanmean(cat4d,4) ./ (std(cat4d,[],4)./sqrt(length(idx)));
    tstats(isnan(tstats)) = 0;
    % save tstats:
    infotmp = info{1};
    infotmp.Filename = [];
    infotmp.Filemoddate = [];
    infotmp.Filesize = [];
    descrip = sprintf('t-stats:%s',conditions{i});
    infotmp.Description = descrip;
    infotmp.raw.descrip = descrip;
    niftidir = fullfile(outputdir,subj_id);
    niftiwrite(tstats,fullfile(niftidir,sprintf('tmap_%s.nii',replace(conditions{i},":", "-"))), infotmp);
end

%% Define tmaps
