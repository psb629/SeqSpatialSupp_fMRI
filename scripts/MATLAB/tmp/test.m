for s=[1:3 5 6 8:14]
% s = 14;
nTRs = [410*ones(1,10) 401 406 410 404 410 410 385]; % For S11
nTRs = [410*ones(1,16) 385]; % for S09

%     sss_imana('BIDS:move_unzip_raw_func','sn',s);
%     sss_imana('FUNC:make_fmap','sn',s);  
%     sss_imana('FUNC:realign_unwarp','sn',s,'rtm',0,'endTR',[410*ones(1,10) 401 406 410 404 410 410 385]);
%     sss_imana('FUNC:move_realigned_images','sn',s,'prefix','u','rtm',0);
%     sss_imana('FUNC:meanimage_bias_correction','sn',s,'prefix','u','rtm',0);
% end

% end
%     %for s=12

         sss_imana('GLM:design','sn',s,'glm',0,'nTR',nTRs(end));
         sss_imana('GLM:estimate','sn',s,'glm',0,'fig',0);
        sss_imana('HRF:fit','sn',s,'regN',[1:7]);
        load(fullfile('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/ROI',sprintf('S%02d_hrf_fit_glm0.mat',s)));
         sss_imana('GLM:design','sn',s,'glm',2,'nTR',nTRs,'hrf_params',params_after);
         sss_imana('GLM:estimate','sn',s,'glm',2,'fig',0);
        sss_imana('GLM:design','sn',s,'glm',3,'nTR',nTRs,'hrf_params',params_after);
        sss_imana('GLM:estimate','sn',s,'glm',3,'fig',0);
        sss_imana('ROI_get_prewhitened_beta','sn',s);
        sss_imana('ROI_calc_rdm','sn',s);
% %         if s==10
% %             load(fullfile('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/ROI',sprintf('S%02d_hrf_fit_glm0.mat',s)));
% %              sss_imana('GLM:design','sn',s,'glm',2,'hrf_params',params_after);
% %              sss_imana('GLM:estimate','sn',s,'glm',2,'fig',0);
% %         end
% %         
% %        
        sss_imana('GLM:tcontrast','sn',s,'glm',2);
        sss_imana('GLM:psc','sn',s,'glm',2);
%         sss_imana('WB:surf_resample','sn',s);
          sss_imana('WB:vol2surf_indiv','sn',s,'glm',2,'map','psc');
          sss_imana('WB:vol2surf_stats','sn',s,'glm',2);
% end        
% sss_imana('GLM:design','sn',s,'glm',0);
% sss_imana('GLM:estimate','sn',s,'glm',0,'fig',0);
% sss_imana('ROI:redefine','sn',s);

% sss_imana('WB:vol2surf_group','sn',[1 2 3 5 6 7 8 9 10 11],'glm',2);

% for s=[7 8 9 10 11]
%     sss_imana('GLM:design','sn',s,'glm',0);
%     sss_imana('GLM:estimate','sn',s,'glm',0,'fig',0);
%     sss_imana('HRF:fit','sn',s);
% end
    
% for s=[1:3 5:11]
%     sss_imana('HRF:ROI_hrf_get','sn',s,'glm',2,'post',20);
%     figure;
%     for r=1:8
%         subplot(2,4,r);
%         sss_imana('HRF:ROI_hrf_plot','sn',s,'roi',r,'glm',2,'post',20);
%     end
% end

            
            % for s=[8 9]
%     figure;
%     for r=1:10
%         subplot(2,5,r);
%         sss_imana('HRF:ROI_hrf_plot','sn',s,'roi',r,'glm',0,'post',20);
%     end
% end
% % 


% s = 6; g=1;
% sss_imana('GLM:design','sn',s,'glm',g);
% sss_imana('GLM:estimate','sn',s,'glm',g,'fig',0);
% sss_imana('GLM:tcontrast','sn',s,'glm',g);
% for s=[3 5 6]
%     for g=2:4
%         sss_imana('GLM:design','sn',s,'glm',g);
%         sss_imana('GLM:estimate','sn',s,'glm',g,'fig',0);
%         sss_imana('GLM:tcontrast','sn',s,'glm',g);
%      end
% end
%sss_imana('PREP:step0','sn',6);

% sss_imana('GLM:design','sn',2,'glm',1);
% sss_imana('GLM:estimate','sn',2,'glm',1,'fig',0);
% sss_imana('GLM:tcontrast','sn',2,'glm',1);
% sss_imana('GLM:tcontrast','sn',5,'glm',1);

% for s=[1 2 3 5 6]
% %      sss_imana('GLM:design','sn',s,'glm',2);
% %      sss_imana('GLM:estimate','sn',s,'glm',2,'fig',0);
%    sss_imana('GLM:tcontrast','sn',s,'glm',1, 'opt',1);
%   % sss_imana('WB:vol2surf_indiv','sn',s,'glm',1,'map','con');
%    sss_imana('GLM:psc','sn',s,'glm',1);
%    sss_imana('WB:vol2surf_indiv','sn',s,'glm',1,'map','psc');
% end
% sss_imana('WB:vol2surf_group','sn',[1 2 3 5 6],'glm',1);
% sss_imana('WB:vol2surf_stats','sn',[1 2 3 5 6],'glm',1);

% s = 12;
% sss_imana('PREP:step0','sn',s);
% sss_imana('PREP:step1','sn',s);
% sss_imana('PREP:step2','sn',s);
% sss_imana('GLM:design','sn',s,'glm',1);
% sss_imana('GLM:estimate','sn',s,'glm',1,'fig',0);
% 
% sss_imana('GLM:estimate','sn',s,'glm',2,'fig',0);
% sss_imana('GLM:design','sn',s,'glm',3);
% sss_imana('GLM:estimate','sn',s,'glm',3,'fig',0);
% sss_imana('GLM:tcontrast','sn',s,'glm',1, 'opt',0);
% sss_imana('GLM:tcontrast','sn',s,'glm',1, 'opt',1);
% sss_imana('GLM:tcontrast','sn',s,'glm',2, 'opt',0);
% sss_imana('GLM:design','sn',s,'glm',0);
% sss_imana('GLM:estimate','sn',s,'glm',0,'fig',0);