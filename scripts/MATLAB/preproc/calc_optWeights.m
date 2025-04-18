function optC = calc_optWeights(Ris, X)
% function optC = calc_optCont(SPM, con)
% Calculate weighted design matrix by using a design matrix encoding every trial
% X: design matrix encoding every trial excluding offset regressors
% con: original contrast
[nRun, nTr] = size(Ris);
C = zeros(nTr*nRun,2); 
C(find(Ris'==1)-1,1) = 1;
C(:,2) = ones(nTr*nRun,1)-C(:,1);            
optW = inv(C'*X'*X*C)*C'*X'*X;
optC(1,:) = [optW(1,:)  zeros(1,nRun)];

C = zeros(nTr*nRun,2); 
C(find(Ris'==1),1) = 1;
C(:,2) = ones(nTr*nRun,1)-C(:,1);            
optW = inv(C'*X'*X*C)*C'*X'*X;
optC(2,:) = [optW(1,:)  zeros(1,nRun)];

optC(3,:) = optC(2,:)-optC(1,:);
% load(fullfile('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/glm_1',  sprintf('S%02d', sn), 'SPM.mat')
% nTr = 68; nRun = 8;
% if size(SPM.xX.X,2) ~= nTr*nRun+nRun
%     error('Load trial-wise design matrix');
% end
% X = SPM.xX.X(:,1:end-length(SPM.nscan));
% optW = inv(C'*X'*X*C)*C'*X'*X;
