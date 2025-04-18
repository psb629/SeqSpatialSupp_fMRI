% for s = [1:3 5 6 8:14]
%     sss_imana('WB:surf_resample','sn',s);
% end
% nTRs = [410*ones(1,10) 401 406 410 404 410 410 385]; % For S11
% nTRs = [410*ones(1,16) 385]; % for S09
nTRs = repmat([410*ones(1,16) 385],12,1);
nTRs(9,:) = [410*ones(1,10) 401 406 410 404 410 410 385];
subjlist = [1 2 3 5 6 8:14];
subjlist = [8:14];

k=1;
for s=subjlist
%     sss_imana('GLM:design','sn',s,'glm',0,'nTR',nTRs(k,end));
%     sss_imana('GLM:estimate','sn',s,'glm',0,'fig',0);
%     sss_imana('HRF:fit','sn',s,'regN',[1:7]);
%     load(fullfile('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/ROI',sprintf('S%02d_hrf_fit_glm0.mat',s)));
    params_after = [6 16];
    if s==11
        nTRs = [410*ones(1,10) 401 406 410 404 410 410 385];
    else
        nTRs = [410*ones(1,16) 385];
    end
    sss_imana('GLM:design','sn',s,'glm',2,'nTR',nTRs,'hrf_params',params_after);
    sss_imana('GLM:estimate','sn',s,'glm',2,'fig',0);

    sss_imana('GLM:tcontrast','sn',s,'glm',2,'opt',0);
    sss_imana('GLM:psc','sn',s,'glm',2);
    sss_imana('WB:vol2surf_indiv','sn',s,'glm',2,'map','psc');

%         sss_imana('WB:surf_resample','sn',s);
%           sss_imana('WB:vol2surf_indiv','sn',s,'glm',1,'map','psc');
%     sss_imana('GLM:tcontrast','sn',s,'glm',4);
%     sss_imana('WB:vol2surf_indiv','sn',s,'glm',2,'map','con');
%    sss_imana('ROI:redefine','sn',s);
    k = k + 1;
end

sss_imana('WB:vol2surf_group','sn',[1 2 3 5 6 8:14],'glm',2,'map','psc');
sss_imana('WB:vol2surf_stats','glm',2,'map','psc');
% sss_imana('WB:vol2surf_group','sn',[1 2 3 5 6 8:14],'glm',0,'map','con');
% sss_imana('WB:vol2surf_stats''glm',0,'map','con');
