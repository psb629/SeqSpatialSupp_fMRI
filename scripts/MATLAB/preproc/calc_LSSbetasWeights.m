function optC = calc_LSSbetasWeights(SPM)
% function optC = calc_optCont(SPM, con)
% Calculate weighted design matrix by using a design matrix encoding every trial
% X: design matrix encoding every trial excluding offset regressors
% con: original contrast

nTr = 68; nRun = 8;
if size(SPM.xX.X,2) ~= nTr*nRun+nRun
    error('Load trial-wise design matrix');
end
X = SPM.xX.X(:,1:end-length(SPM.nscan));
CC = zeros(nTr,nRun*nTr,2);  %% did not include offset regressors
optC = zeros(nRun*nTr, nRun*nTr+nRun);
for i=1:nTr*nRun
    CC(i,:,2)=1;
    CC(i,i,1) = 1;
    CC(i,i,2) = 0;
    C = squeeze(CC(i,:,:));
    optW(i,:,:)  = inv(C'*X'*X*C)*C'*X'*X;
    optC(i,1:nRun*nTr) = squeeze(optW(i,1,:));  %% row vector is a contrast
end
fprintf('Process done...\n')



% optC = zeros(2,size(SPM.xX.X,2));
% X = SPM.xX.X(:,1:end-length(SPM.nscan));
% C_motor = zeros(size(X,2),length(find(unique(con)))); % Number of rows of design matrix
% 
% idx{1} = find(con(1,:)==1);
% idx{2} = find(con(1,:)==-1);
% idx{3} = find(con(1,:)==0);
% 
% for i=1:length(idx)
%     C_motor(idx{i},i)=1;
% end
% 
% C = C_motor;
% optW(1,:,:)  = inv(C'*X'*X*C)*C'*X'*X;
% optC(1,1:size(X,2))=squeeze(optW(1,1,:)-optW(1,2,:));
% 
% 
% idx{1} = find(con(2,:)==1);
% idx{2} = find(con(2,:)==-1);
% idx{3} = find(con(2,:)==0);
% 
% for i=1:length(idx)
%     C_cue(idx{i},i)=1;
% end
% 
% 
% 
% C = C_cue;
% optW(2,:,:) = inv(C'*X'*X*C)*C'*X'*X;
% optC(2,1:1:size(X,2))=squeeze(optW(2,1,:)-optW(2,2,:));
% 
