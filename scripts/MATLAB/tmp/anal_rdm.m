close all;
% trialState ((seqID, seqType)=0:(0,0),1:(0,1),2:(1,0),3:(1,1),4:(2,0),5:(2,1),6:(3,0),7:(3,1)
% indices before sorting
samefingIdx = [1 2 3 8 9 14 23 24 25 26 27 28];
difffingIdx = [4:7 10:13 15:22];
samecueIdx = [2 4 6 9 11 13 15 17 20 22 24 27];
diffcueIdx = [1 3 5 7 8 10 12 14 16 18 19 21 23 25 26 28];
sameseqIdx = [1 14 23 28];
diffseqIdx = [2:13 15:22 24:27];

diffseqIdxCue = [2 4 6 9 11 13 15 17 20 22 24 27];
diffcueIdxSeq = [1 14 23 28];
samefingIdxCue = [2 9 24 27];
difffingIdxCue = [4 6 11 13 15 17 20 22];

subjlist = [1:3 5 6 8:14];
group = 'S_all';
% subjlist = subjlist + 14;
% subjlist(subjlist==25) = [];
% subjlist(subjlist==27) = [];
% group = 'R_all';
% for i=1:length(subjlist)
%     sss_imana('ROI_get_prewhitened_beta','sn',subjlist(i));
%     sss_imana('ROI_calc_rdm','sn',subjlist(i));
%     fprintf('Processed S%02d\n',subjlist(i));
% end

pinfo = dload('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/participants.tsv');
dir_rdm = '/srv/diedrichsen/data/SeqSpatialSupp_fMRI/patterns';
cd(dir_rdm);
for i=1:length(subjlist)
    sn = subjlist(i);
    subj_id = char(pinfo.subj_id(pinfo.sn==sn));
    load(sprintf('%s_fRDM.glm3.mat',subj_id));
    for r=1:9
        mDissFing(i,r) = mean(Data.RDM{r}(difffingIdx))-mean(Data.RDM{r}(samefingIdx));
        mDissFing_wCue(i,r) = mean(Data.RDM{r}(difffingIdxCue))-mean(Data.RDM{r}(samefingIdxCue));
        mDissCue(i,r) = mean(Data.RDM{r}(diffcueIdx))-mean(Data.RDM{r}(samecueIdx));
        mDissSeq(i,r) = mean(Data.RDM{r}(diffseqIdx))-mean(Data.RDM{r}(sameseqIdx));
%         mDissSeq_Cue(i,r) = mean(Data.RDM{r}(diffseqIdx))-mean(Data.RDM{r}(sameseqIdx));
        mDissSeq_wCue(i,r) = mean(Data.RDM{r}(diffseqIdxCue));
        mDissCue_wSeq(i,r) = mean(Data.RDM{r}(diffcueIdxSeq));
        mDiss(i,r) = mean(Data.RDM{r});
        mDissSameFing_wCue(i,r) = mean(Data.RDM{r}(samefingIdxCue));
    end
end

% figure;
% plot_rdm_roi(mDissSeq_wCue);

sort_idx = [2 4 6 1 3 5 7 15 17 8 14 16 18 24 10 19 23 25 12 21 26 28 9 11 13 20 22 27];  %% mapping from previous to the current
for i=1:length(subjlist)
    sn = subjlist(i);
    subj_id = char(pinfo.subj_id(pinfo.sn==sn));
%     load(sprintf('S%02d_fRDM.glm3.mat',subjlist(i)));
    load(sprintf('%s_fRDM.glm3.mat',subj_id));
    for r=1:11
       temp = rsa.rdm.squareRDM(Data.RDM{r}(sort_idx));       
       RDM_all(r,i,:,:) = temp;
       temp = rsa.rdm.vectorizeRDM(temp);
       diss_rep_across(r,i) = mean(temp([4 11 17 22]));
       diss_nrep_across(r,i) = mean(temp([5:7 10 12 13 15 16 18:21]));
       diss_rep_within_L(r,i) = mean(temp([1 2 3 8 9 14]));
       diss_rep_within_S(r,i) = mean(temp([23:28]));
       diss_repf_within_L(r,i) = mean(temp([1 14]));
       diss_repf_within_S(r,i) = mean(temp([23 28]));
       diss_nrepf_within_L(r,i) = mean(temp([2 3 8 9]));
       diss_nrepf_within_S(r,i) = mean(temp([24:27]));
       diss_repf_across(r,i) = mean(temp([4 5 10 11 17 18 21 22]));
       diss_nrepf_across(r,i) = mean(temp([6 7 12 13 15 16 19 20]));       
    end
end

figure;subplot(431);
plot_rdm_roi(diss_rep_within_L,'b');
% ylim([0 0.7]);
title('Dissimilarity of different sequences within letter cue','Rotation',5,'Fontsize',10); ylabel('Mean dissimilarity (a.u.)','Fontsize',10);
h = gca; h.XAxis.TickLabelRotation=45; h.XAxis.FontSize=9;

subplot(432);
plot_rdm_roi(diss_rep_within_S,'b');
% ylim([-0.1 0.7]);
title('Dissimilarity of different sequences within spatial cue','Rotation',5,'Fontsize',10); ylabel('Mean dissimilarity (a.u.)','Fontsize',10);
h = gca; h.XAxis.TickLabelRotation=45; h.XAxis.FontSize=9;

subplot(433);
plot_rdm_roi(diss_rep_within_L-diss_rep_within_S,'b');
% ylim([-0.1 0.1]);
title('Dissimilarity of different sequence (Letter-Spatial)','Rotation',5,'Fontsize',10); ylabel('Mean dissimilarity (a.u.)','Fontsize',10);
h = gca; h.XAxis.TickLabelRotation=45; h.XAxis.FontSize=9;

subplot(434);
plot_rdm_roi(diss_rep_across,'b');
% ylim([-0.1 0.7]);
title('Dissimilarity of the same sequence across cues','Rotation',5,'Fontsize',10); ylabel('Mean dissimilarity (a.u.)','Fontsize',10);
h = gca; h.XAxis.TickLabelRotation=45; h.XAxis.FontSize=9;

subplot(435);
plot_rdm_roi(diss_nrep_across,'b');
% ylim([-0.1 0.7]);
title('Dissimilarity of the different sequence across cues','Rotation',5,'Fontsize',10); ylabel('Mean dissimilarity (a.u.)','Fontsize',10);
h = gca; h.XAxis.TickLabelRotation=45; h.XAxis.FontSize=9;

subplot(436);
plot_rdm_roi(diss_nrep_across-diss_rep_across,'b');
% ylim([-0.1 0.1]);
title('Dissimilarity of (Diff-Same) sequence across cues','Rotation',5,'Fontsize',10); ylabel('Mean dissimilarity (a.u.)','Fontsize',10);
h = gca; h.XAxis.TickLabelRotation=45; h.XAxis.FontSize=9;

subplot(437);plot_rdm_roi(diss_repf_within_L,'b');
% ylim([-0.1 0.7]);
title('Mean dissimilarity (same finger) within letter cue','Rotation',5,'Fontsize',10); ylabel('Mean dissimilarity (a.u.)','Fontsize',10);
h = gca; h.XAxis.TickLabelRotation=45; h.XAxis.FontSize=9;

subplot(438);
plot_rdm_roi(diss_nrepf_within_L,'b');
% ylim([-0.1 0.7]);
title('Mean dissimilarity (different finger) within letter cue','Rotation',5,'Fontsize',10); ylabel('Mean dissimilarity (a.u.)','Fontsize',10);
h = gca; h.XAxis.TickLabelRotation=45; h.XAxis.FontSize=9;

subplot(439);
plot_rdm_roi(diss_nrepf_within_L-diss_repf_within_L,'b');
% ylim([-0.1 0.1]);
title('Mean dissimilarity (Diff-Same) within letter cue','Rotation',5,'Fontsize',10); ylabel('Mean dissimilarity (a.u.)','Fontsize',10);
h = gca; h.XAxis.TickLabelRotation=45; h.XAxis.FontSize=9;

subplot(4,3,10);plot_rdm_roi(diss_repf_within_S,'b');
% ylim([-0.1 0.7]);
title('Mean dissimilarity (same finger) within spatial cue','Rotation',5,'Fontsize',10); ylabel('Mean dissimilarity (a.u.)','Fontsize',10);
h = gca; h.XAxis.TickLabelRotation=45; h.XAxis.FontSize=9;

subplot(4,3,11);
plot_rdm_roi(diss_nrepf_within_S,'b');
% ylim([-0.1 0.7]);
title('Mean dissimilarity (different finger) within spatial cue','Rotation',5,'Fontsize',10); ylabel('Mean dissimilarity (a.u.)','Fontsize',10);
h = gca; h.XAxis.TickLabelRotation=45; h.XAxis.FontSize=9;

subplot(4,3,12);
plot_rdm_roi(diss_nrepf_within_S-diss_repf_within_S,'b');
% ylim([-0.1 0.1]);
title('Mean dissimilarity (Diff-Same) within spatial cue','Rotation',5,'Fontsize',10); ylabel('Mean dissimilarity (a.u.)','Fontsize',14);
h = gca; h.XAxis.TickLabelRotation=45; h.XAxis.FontSize=9;

% saveas(gcf, sprintf('%s/plot.anal_rdm.%s.png',dir_rdm,group));
