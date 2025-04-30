T = load('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/ROI/psc_glm2_Task_N=12.mat');
nR = 11;  nC = 8; nS = 12;
for i=1:nS 
    psc(i,:,:) =reshape(T.psc(nC*nR*(i-1)+1:nC*nR*i),nC,nR);
end
% 
% sort_idx = [1:2:7 2:2:8 17:2:23 18:2:24 33:2:39 34:2:40 49:2:55 50:2:56 9:2:15 10:2:16 25:2:31 26:2:32 41:2:47 42:2:48 57:2:63 58:2:64]; 
% psc_sort = psc(:,sort_idx,:);
% labels = {'L-SMA','L-PMv','L-PMd','L-M1','L-S1','L-SPLa','L-SPLp'};
% labels = {'L-SMA','L-PMv','L-PMd','L-M1','L-S1','L-SPLa','L-SPLp','L-DSVC','L-MT+','L-VSVC','L-EAC'};
% O = {'psc_BothRep-L.nii','psc_CueRep-L.nii','psc_MotorRep-L.nii','psc_NRep-L.nii',...
%                  'psc_BothRep-S.nii','psc_CueRep-S.nii','psc_MotorRep-S.nii','psc_NRep-S.nii'};
figure;
for r=1:nR
    subplot(3,4,r);
    meandata = squeeze(mean(psc(:,:,r),1));
    semdata = std(squeeze(psc(:,:,r)))/sqrt(nS);
    bar(meandata); hold on;
    errorbar(meandata,semdata,'k', 'linestyle', 'none');
    % Customize plot
%     set(gca, 'XTick', 1:11, 'XTickLabel', labels);
    xlabel('Regions of Interest','Fontsize',14);
    ylabel('RS effect (psc)','Fontsize',14);%     imagesc(reshape(temp(:,r),8,8)');
%     colorbar;
end    
  
    

% subjlist = [1:3 5 6 8:14];
% for i=1:length(subjlist)
%     subjlist(i)