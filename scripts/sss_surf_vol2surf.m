function varargout = sss_surf_vol2surf(s, glm, map)
%%
addpath(genpath('D:/milli/diedrichsenlab/imaging_tools'));
addpath(genpath('D:/milli/diedrichsenlab/matlab'));
%%
dir_git = 'D:/mobaxterm/sungbeenpark/github';
pinfo = dload(fullfile(dir_git,'diedrichsenlab/SeqSpatialSupp_fMRI/participants.tsv'));
subj_id = char(pinfo.subj_id(pinfo.sn==s));
S_id = strrep(subj_id,'R','S');
%%
dir_work = 'D:/milli/diedrichsenlab/SeqSpatialSupp_fMRI';
dir_surf = fullfile(dir_work,'surfaceWB',S_id);
dir_glm = fullfile(dir_work,sprintf('glm_%d',glm));
%%
surf = '32k';
hemis = [1 2];
hem = {'L','R'};
for h = hemis
    white = fullfile(dir_surf, sprintf('%s.%s.white.%s.surf.gii', S_id, hem{h}, surf));
    pial = fullfile(dir_surf, sprintf('%s.%s.pial.%s.surf.gii', S_id, hem{h}, surf));
    C1 = gifti(white);
    C2 = gifti(pial);
    
	load(fullfile(subjGLM, 'SPM.mat')); 
    switch map
        case 't' % t-values maps (univariate GLM)
            fnames = cell(1,numel(SPM.xCon));
            con_name = cell(1,numel(SPM.xCon));
            for j=1:numel(fnames)
            	fnames{j} = fullfile(subjGLM, sprintf('spmT_%s.nii', SPM.xCon(j).name));
                con_name{j} = SPM.xCon(j).name;
            end
        case 'con' % contrast beta maps (univariate GLM)
            fnames = cell(1,numel(SPM.xCon));
            con_name = cell(1,numel(SPM.xCon));
            for j=1:numel(fnames)
            	fnames{j} = fullfile(subjGLM, sprintf('con_%s.nii', SPM.xCon(j).name));
            	con_name{j} = SPM.xCon(j).name;
            end    
        case 'psc' % percent signal change maps (univariate GLM)
            fnames = cell(1,numel(SPM.xCon));
            con_name = cell(1,numel(SPM.xCon));
        	for j=1:numel(fnames)
				fnames{j} = fullfile(subjGLM, sprintf('psc_%s.nii', SPM.xCon(j).name));
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