%%
if ispc
    cd '\\wsl.localhost/ubuntu-24.04/home/sungbeenpark/github/SeqSpatialSupp_fMRI/scripts/MATLAB'
elseif ismac
    cd '/Users/sungbeenpark/github/SeqSpatialSupp_fMRI/scripts/MATLAB'
end
sss_init;

dir_PALM = fullfile(dir_git,'PALM');
addpath(fullfile(dir_PALM));

dir_FS = fullfile(dir_git,'fs_LR_32');
addpath(fullfile(dir_FS));
midthickness = fullfile(dir_FS,'fs_LR.32k.L.midthickness.surf.gii');

dir_glmsingle = fullfile(baseDir,'GLMsingle');

%%
dir_work = fullfile(dir_glmsingle,'surfaceWB/glm_1');
cd(dir_work);

input = fullfile(dir_work,'cifti.L.glm_1.searchlight.mean_dist.first_finger.dscalar.nii');
design = fullfile(dir_SSS,'scripts/PALM/glm_1.searchlight/design.csv');
contrast = fullfile(dir_SSS,'scripts/PALM/glm_1.searchlight/contrast.csv');

%%
cmd = sprintf('wsl palm -i %s -s %s -d %s -t %s -n 5000 -C 3.1 -Cstat extent -o palm_RDM_L_cluster', input, midthickness, design, contrast);
system(cmd);