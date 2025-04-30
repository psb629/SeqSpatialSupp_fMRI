function optC = calc_optWeights2(Ris, X)
% function optC = calc_optCont(SPM, con)
% Calculate weighted design matrix by using a design matrix encoding every trial
% X: design matrix encoding every trial excluding offset regressors
% con: original contrast
[nRun, nTr] = size(Ris);
C = zeros(nTr*nRun,2); 
C(find(Ris'==1),1) = 1;
C(:,2) = ones(nTr*nRun,1)-C(:,1);            
optW = inv(C'*X'*X*C)*C'*X'*X;
optC(1,:) = [optW(1,:)  zeros(1,nRun)];
