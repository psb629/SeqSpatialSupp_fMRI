function spmj_normalization_write(deformation_file, images,varargin)
% function spmj_normalization_write(deformation_file, images)
% INPUTS:
%   deformation_file:     the *_sn.mat file for the deformation
%   images:               a cell array of the images you want to resmape
% VARARGIN:
%   'outimages',cell:  names (and path) under which you want to store
%                         the resampled images. Otherwise they will be
%                         stored under the same directory with prefix 'w'

%12.09.2012 TW: change bounding box from [-78 -112 -50 78   76  85] to [-78 -112 -50 78   76  100]
if ischar(images)
    images=cellstr(images);
end;
outimages=[];

vararginoptions(varargin,{'outimages'});

J.subj.matname = {deformation_file};
J.subj.resample = images;
J.roptions.preserve = 0;
J.roptions.bb = [-78 -112 -65
    78   76  85];
J.roptions.vox = [2 2 2];
J.roptions.interp = 1;
J.roptions.wrap = [0 0 0];
J.roptions.prefix = 'w';

matlabbatch{1}.spm.spatial.normalise.write=J;
spm_jobman('run',matlabbatch);

if (~isempty(outimages))
    for j=1:length(images)
        % Move the image. If output is a one-part nii, automatically joint
        % the two parts 
        [dir,name,ext,num]=spm_fileparts(images{j});
        tempname=fullfile(dir,[J.roptions.prefix name ext num]); 
        V=spm_vol(tempname); 
        X=spm_read_vols(V); 
        V.fname=outimages{j}; 
        spm_write_vol(V,X); 
        delete(tempname); 
    end;
end;
