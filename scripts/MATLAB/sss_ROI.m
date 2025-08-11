function varargout = sss_ROI(what,varargin)

if ispc
    cd '\\wsl.localhost/ubuntu-24.04/home/sungbeenpark/github/SeqSpatialSupp_fMRI/scripts/MATLAB'
    dir_atlas = 'F:/Atlas/fs_LR_32';
elseif ismac
    cd '/Users/sungbeenpark/github/SeqSpatialSupp_fMRI/scripts/MATLAB'
    dir_atlas = '/Volumes/Diedrichsen_data$/data/Atlas_templates/fs_LR_32';
end

sss_init;

%% argument inputs
TR = 1;
sn = [];
glm = 1;
vararginoptions(varargin,{'sn','glm'});
if isempty(sn)
    error('''sn'' must be passed to this function.')
end
[subj_id, S_id] = get_id(fullfile(baseDir,'participants.tsv'), sn);

if isempty(glm)
    error('''glm'' must be passed to this function.')
end
glmDir = sprintf('glm_%d',glm);

%% ROI
hemi = {'lh','rh'}; % left & right hemi folder names/prefixes
hem = {'L', 'R'}; % hemisphere: 1=LH 2=RH
hname = {'CortexLeft', 'CortexRight'}; % 'CortexLeft', 'CortexRight', 'Cerebellum'

%% MAIN OPERATION 
switch(what)
    case 'ROI:init'
        if subj_id(1)=='S'
            sss_ROI('ROI:calc_region','sn',sn); % https://github.com/DiedrichsenLab/region.git
        end
        % sss_hrf('ROI:deform','sn',sn,'glm',glm);
        sss_ROI('ROI:make_cifti.y_raw','sn',sn); % https://github.com/Washington-University/cifti-matlab.git

    case 'ROI:glm'
        sss_ROI('ROI:make_cifti.ResMS','sn',sn,'glm',glm); % https://github.com/Washington-University/cifti-matlab.git

    % case 'ROI:findall' % use this... %rdm, glm3
    %     % ROI 11개를 모아놓은 ROI.L.SSS.label.gii 파일 만들기
    %     h=1;
    %     surf='32';
    % 
    %     % 
    %     ROI{1,1}={'L_SCEF_ROI','L_6ma_ROI','L_6mp_ROI'}; % SMA
    %     ROI{1,2}={'L_6v_ROI','L_PEF_ROI','L_55b_ROI'}; % PMv
    %     ROI{1,3}={'L_6a_ROI','L_6d_ROI','L_FEF_ROI'}; % PMd
    %     % I took M1 and S1 from the one that we have already
    %     ROI{1,6}={'L_AIP_ROI','L_7PC_ROI','L_LIPv_ROI','L_LIPd_ROI'}; % SPLa
    %     ROI{1,7}={'L_MIP_ROI','L_VIP_ROI','L_7PL_ROI'}; % SPLp
    %     ROI{1,8}={'L_IPS1_ROI','L_V3A_ROI','L_V3B_ROI','L_V7_ROI','L_V6A_ROI'}; % DSVC ,'L_V6_ROI'
    %     ROI{1,9}={'L_LO1_ROI','L_LO2_ROI','L_LO3_ROI','L_V3CD_ROI','L_V4t_ROI','L_FST_ROI','L_MT_ROI','L_MST_ROI','L_PH_ROI'}; % MT+
    %     ROI{1,10}={'L_FFC_ROI','L_VVC_ROI','L_V8_ROI','L_VMV1_ROI','L_VMV2_ROI','L_VMV3_ROI','L_PIT_ROI'}; % VSVC
    %     ROI{1,11}={'L_RI_ROI','L_MBelt_ROI','L_PBelt_ROI','L_A1_ROI','L_LBelt_ROI'}; % EAC
    %     % ROI{1,12}={}; % V1
    % 
    %     ROI_name={'','SMA','PMv','PMd','M1','S1','SPLa','SPLp','DSVC','MT+','VSVC','EAC'};  % important! add the first ''
    % 
    %     pathtosurf=fullfile(atlasDir,sprintf('FS_LR_%sk',surf));
    %     P_glasser=gifti(fullfile(pathtosurf,sprintf('Glasser_2016.%sk.%s.label.gii',surf,hem{h})));
    %     P_brodmann=gifti(fullfile(pathtosurf,sprintf('ROI.%sk.%s.label.gii',surf,hem{h})));
    %     P=zeros(size(P_glasser.cdata));
    % 
    %     for i=1:length(ROI_name)-1
    %         if (i==4||i==5) % M1 and S1
    %             if i==4
    %                 P(ismember(P_brodmann.cdata,2),1)=i;
    %             else
    %                 P(ismember(P_brodmann.cdata,1),1)=i;
    %             end
    %         else
    %             roi_num=[];
    %             for j=1:length(ROI{1,i})
    %                 roi_num=cat(2,roi_num,P_glasser.labels.key(strcmp(P_glasser.labels.name,ROI{1,i}{j})));
    %             end
    %             P(ismember(P_glasser.cdata,roi_num),1)=i;
    %         end
    %     end
    % 
    %     % create a label
    %     G=surf_makeLabelGifti(P,'labelNames',ROI_name);
    %     fname=fullfile(atlasDir,sprintf('fs_LR_32k/ROI.%s.%s.label.gii',hem{h},'SSS'));
    %     if ~exist(fname,'file')
    %         save(G,fname,'-v7.3');
    %     end
    %     fprintf(1,'Done.\n');
    % 
    %     varargout={P,ROI_name};

    case 'ROI:calc_region' % use this... %rdm, glm3
        % ROI template를 피험자의 공간(2D surface)에 매핑

        % [P,ROI_name]=sss_hrf('ROI:findall');
        % n_roi=length(ROI_name)-1;
        % 왜 그룹 마스크? 개인 마스크가 아니고? -> # of voxels 수를 맞추기 위해?
        %PP = load('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/surfaceWB/group32k/fmask_data.mat');
        %P = (PP.mask'~=0).*P;

        %% load Atlas
        atlasH = {'ROI.32k.L.label.gii', 'ROI.32k.R.label.gii'};
        atlas_gii = {gifti(fullfile(dir_atlas, atlasH{1})), gifti(fullfile(dir_atlas, atlasH{2}))};

        fprintf('%s...\n',subj_id);
        
        %% load mask.nii from GLM result
        Vol = fullfile(baseDir,glmDir,subj_id,'mask.nii');

        dir_surf = fullfile(baseDir,wbDir,S_id);
        Hem = {'L','R'};
        R = {};
        idx = 1;
        for h=[1 2]
            for reg = 2:length(atlas_gii{h}.labels.name)
                % pathtosurf=fullfile(atlasDir,sprintf('fs_LR_%sk',surf));
                % surface=fullfile(pathtosurf,sprintf('fs_LR.%sk.%s.flat.surf.gii',surf,hem{h}));
                % surface = fullfile('F:/Atlas/fs_LR_32/fs_LR.32k.L.flat.surf.gii');
                R{idx}.type     = 'surf_nodes_wb';
                R{idx}.white    = fullfile(dir_surf,sprintf('%s.%s.white.32k.surf.gii',S_id,hem{h})); % T1 space
                R{idx}.pial     = fullfile(dir_surf,sprintf('%s.%s.pial.32k.surf.gii',S_id,hem{h}));  % T1 space
                % R{idx}.flat     = surface;
                R{idx}.linedef  = [5,0,1];  % take 5 steps along node between white (0) and pial (1) surfaces
                R{idx}.image    = Vol;     % functional space
                R{idx}.name     = atlas_gii{h}.labels.name{reg};
                R{idx}.hem      = Hem{h};
                key             = atlas_gii{h}.labels.key(reg);
                R{idx}.location = find(atlas_gii{h}.cdata==key);
                idx = idx+1;
            end
        overlap = [1 2; 1 3; 1 4; 3 4; 7 8; 1 7];
        end

        % calculates the locations for regions of interest
        R = region_calcregions(R, 'exclude', [overlap; overlap+8], 'exclude_thres', .8);

        dir_work = fullfile(baseDir,roiDir,subj_id);
        if (~exist(dir_work,'dir'))
            mkdir(dir_work);
        end
        fname=fullfile(dir_work,sprintf('%s.Task_regions.mat',subj_id));
        %fname=fullfile(baseDir,roiDir,sprintf('%s_SSS_regions.mat',sprintf('S%02d',sn)));

        %% save R
        save(fname,'R','-v7.3');
        
        %% save R as an image file (Is it same as the deformation?)
        for r = 1:length(R)
            img = region_saveasimg(R{r}, Vol, 'name', ...
                fullfile(dir_work, sprintf('ROI.%s.%s.%s.nii',R{r}.hem,subj_id,R{r}.name)));
        end

        fprintf('\nROIs have been defined for %s \n',subj_id);
        % end  

    case 'ROI:deform'
        % 피험자의 공간에 매핑된 ROI.gii(surface)들을 다시 피험자의 Volume에 매핑
        fprintf('Deformation for %s...\n',subj_id);

        % workDir = fullfile(baseDir,roiDir);
        % fname = fullfile(workDir,sprintf('%s_Task_regions.glm_%d.mat',subj_id,glm));
        dir_work = fullfile(baseDir,roiDir,glmDir,subj_id);
        fname = fullfile(dir_work,sprintf('%s.Task_regions.glm%d.mat',subj_id,glm));
        R = load(fname); % load the variable R
        R = R.R;

        dir_anat = fullfile(baseDir,anatomicalDir,subj_id);
        deffile = fullfile(dir_anat,sprintf('iy_%s_anatomical.nii',subj_id));
        R1 = region_deformation(R, deffile, 'mask', R{1}.image); fprintf('done!\n');
        % fname = fullfile(workDir,sprintf('%s.Task_regions.glm%d.mat',subj_id,glm));
        % save(fname,'R1','-v7.3');

        % VolFile = fullfile(anatDir,sprintf('%s_anatomical.nii',subj_id));
        VolFile = fullfile(baseDir,glmDir,subj_id,'mask.nii');
        Vol = spm_vol(VolFile);
        varargout={R1, Vol};

        % save
        cd(dir_work);
        for i = 1:length(R1)
            R1{i}.name = sprintf('deform.%s.%s', R1{i}.hem, R1{i}.name);
            region_saveasimg(R1{i}, VolFile);
        end

    case 'ROI:make_cifti.y_raw'
        % 피험자의 volume 공간으로 매핑한 ROI 마스크를 이용하여 피험자의
        % y_raw 데이터를 voxel 단위로 얻고, 그 결과를 cifti로 저장

        % 피험자의 surface 공간에서 각 ROI의 node (2-D) 정보 
        fname = fullfile(baseDir,roiDir,subj_id,sprintf('%s.Task_regions.mat',subj_id));
        % [R, V] = sss_hrf('ROI:deform','sn',sn,'glm',glm,'LR',LR);
        R = load(fname); R = R.R;

        % 피험자 EPI 의 3-D 정보
        % VolFile = fullfile(baseDir,glmDir,subj_id,'mask.nii');
        VolFile = R{1,1}.image;
        V = spm_vol(VolFile);

        SPM = load(fullfile(baseDir,glmDir,subj_id,'SPM.mat'));
        SPM = SPM.SPM;

        fprintf('extrating Y_raw from each ROI for subject %s...\n',subj_id);
        % R의 정보를 토대로 추출한 2-D y data
        D = region_getdata(SPM.xY.VY,R);
        
        % dir_work = fullfile(baseDir, roiDir, glmDir);
        % if ~exist(dir_work,'dir')
        %     mkdir(dir_work);
        % end

        % length(D) = the number of ROIs (Left/Right)
        for i = 1:length(D)
            name = R{1,i}.name;
            hemisphere = R{1,i}.hem;

            %% y_raw
            fname = fullfile(baseDir,roiDir,subj_id,sprintf('cifti.%s.%s.%s.y_raw.dtseries.nii',hemisphere,subj_id,name));
            if ~exist(fname,'file')
                cii = region_make_cifti(R{i},V,'data',D{i}','dtype','series','TR',TR);
                cifti_write(cii, fname);
            end
            
            % %% others
            % [beta, Yhat, Yres] = spmj_glm_fit(SPM,D{i});
            % 
            % % y_hat
            % cii = region_make_cifti(R{i},V,'data',Yhat','dtype','series','TR',TR);
            % fname = fullfile(dir_work, sprintf('cifti.%s.glm%d.%s.%s.y_hat.dtseries.nii',hemisphere,glm,subj_id,name));
            % 
            % cifti_write(cii, fname);
            % 
            % % y_res
            % cii = region_make_cifti(R{i},V,'data',Yres','dtype','series','TR',TR);
            % fname = fullfile(dir_work, sprintf('cifti.%s.glm%d.%s.%s.y_res.dtseries.nii',hemisphere,glm,subj_id,name));
            % 
            % cifti_write(cii, fname);
        end

    case 'ROI:make_cifti.ResMS'
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

        dir_work = fullfile(baseDir,glmDir,subj_id);

        SPM = load(fullfile(dir_work,'SPM.mat'));
        SPM = SPM.SPM;

        fprintf('extrating Y_raw from each ROI for subject %s...\n',subj_id);
        SPM.VResMS.fname = fullfile(dir_work,'ResMS.nii');
        % R의 정보를 토대로 추출한 2-D y data
        D = region_getdata(SPM.VResMS,R);

        % length(D) = the number of ROIs (Left/Right)
        dir_output = fullfile(baseDir,roiDir,glmDir);
        if ~exist(dir_output,'dir')
            mkdir(dir_output);
        end
        area = 'OTHER';
        for i = 1:length(D)
            name = R{1,i}.name;
            hemisphere = R{1,i}.hem;
            
            %% beta
            fname = fullfile(dir_output,sprintf('cifti.%s.%s.%s.ResMS.dscalar.nii',hemisphere,subj_id,name));
            if ~exist(fname,'file')
                cii = region_make_cifti(R{i},V,'data',D{i}','dtype','scalars','struct',area,'TR',TR);
                cifti_write(cii, fname);
            end
            clear cii
        end

end




 





