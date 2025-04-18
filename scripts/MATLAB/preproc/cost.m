function err=cost(p, Y,SPM)
p_hrf = [6 16 1 1 6 0 32];
p_hrf(1:2) = p(1:2);
p_hrf(6) = p(3);
p_hrf(end) = p(4);
% p_hrf = [6 16 1 1 6 0 32];
% p_hrf(1:2) = p(1:2);
% p_hrf(end) = p(3);
% p_hrf = p;
temp = spm_hrf_modified(SPM.xY.RT, p_hrf);
SPM.xBF.bf = temp/sum(temp);
SPM = fMRI_design_changeBF(SPM);
res = spm_sp('r',SPM.xX.xKXs,Y); % get the residual
% constraints
err = sum(sum(res.^2))/numel(res);

        