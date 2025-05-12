function [SPM,Yhat,Yres, p_opt]=fit_spm_hrf(SPM,Yraw)
% Fits the two-gamma hrf fucntion to extracted time 
% function [hrf,P,xBF]=spmj_fit_hrf(SPM,Y_raw,options);
% INPUT: 
%   SPM:    SPM-structure with model to be fitted 
%   Yraw:   raw time series for hrf fitting 
% OUTPUT: 
%   p:      parameters fitted 
%   SPM:    changed SPM
%   Yhat:   predicted time series 
%   Yres:   residual time series 
% p    - parameters of the response function (two Gamma functions)
%                                                           defaults
%                                                          {seconds}
%        p(1) - delay of response (relative to onset)          6
%        p(2) - delay of undershoot (relative to onset)       16
%        p(3) - dispersion of response                         1
%        p(4) - dispersion of undershoot                       1
%        p(5) - ratio of response to undershoot                6
%        p(6) - onset {seconds}                                0
%        p(7) - length of kernel {seconds}                    32

% --------------------------------------------------------------
% fit = [];
% vararginoptions(varargin,{'P0','fit'}); 

% TODO make is more general
% Fill in the default duration of the event, if not given 

SPM.xBF.T        = 16; % time bins per scan
SPM.xBF.T0       = 1; % first time bin
SPM.xBF.UNITS    = 'secs'; % units of the 
SPM.xBF.name     = 'fitted_hrf';
% TODO
SPM.xBF.length   = 50.0625; % support in seconds 
SPM.xBF.order    = 1;
SPM.xBF.Volterra = 1;  % volterra expansion order?
SPM.xBF.dt       = SPM.xY.RT/16;

% Filter and prepare the data 
Y = spm_filter(SPM.xX.K,SPM.xX.W*Yraw);


% optimization options
option.TolX = 1E-10;
option.Display = 'off';
option.MaxFunEvals = 500;


options = optimoptions('fmincon','Algorithm','interior-point'); % run interior-point algorithm
% X = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON) X = fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON)
% 
% p0 = [6 16 1 1 6 0 32]';  %% default parameters from spm_get_defaults('stats.fmri.hrf');
% LB = [2 10 0 0 0 0 16]; UB = [10 20 5 5 10 5 40];
% % LB = p0; UB=p0;
% p_opt = fmincon(@(p) cost(p, Y,SPM), p0, [1 -1 0 0 0 0 0;0 1 0 0 0 0 -1],zeros(2,1), [], [], LB, UB,[],options);
% SPM.xBF.bf = spm_hrf_modified(SPM.xY.RT, p_opt);  % replace the basis function with optimal hrf

p0 = [6 16 0 32]';
p_opt = fmincon(@(p) cost(p, Y,SPM), p0, [1 -1 0 0;0 1 0 -1],zeros(2,1), [], [], [2 5 0 16], [10 20 5 40],[],options);
SPM.xBF.bf = spm_hrf_modified(SPM.xY.RT, [p_opt(1:2)' 1 1 6 p_opt(3:4)']);  % replace the basis function with optimal hrf

% p0 = [6 16 32]';
% LB = [2 5 16]; UB = [10 20 40];
% p_opt = fmincon(@(p) cost(p, Y,SPM), p0, [1 -1 0;0 1 -1],zeros(2,1), [], [], LB, UB,[],options);
% SPM.xBF.bf = spm_hrf_modified(SPM.xY.RT, [p_opt(1:2)' 1 1 6 0 p_opt(3)]);  % replace the basis function with optimal hrf
% p0 = [6 16]';
% LB = [2 5]; UB = [10 20];
% p_opt = fmincon(@(p) cost(p, Y,SPM), p0, [1 -1],0, [], [], LB, UB,[],options);
% SPM.xBF.bf = spm_hrf(SPM.xY.RT/16, [p_opt(1:2)' 1 1 6 0 32],16);  % replace the basis function with optimal hrf


SPM = fMRI_design_changeBF(SPM); 

% return predicted timeseries and residuals 
beta  = SPM.xX.pKX*Y; %-Parameter estimates
Yres  = spm_sp('r',SPM.xX.xKXs,Y); % get the 
reg_interest=[SPM.xX.iH SPM.xX.iC]; 
Yhat   = SPM.xX.xKXs.X(:,reg_interest)*beta(reg_interest,:); %- predicted values 
%Yhat   = SPM.xX.xKXs.X*beta; %- predicted values

function err=cost(p, Y,SPM)
p_hrf = [6 16 1 1 6 0 32];
p_hrf(1:2) = p(1:2);
p_hrf(6) = p(3);
p_hrf(end) = p(4);

% p_hrf(1:2) = p(1:2);
% p_hrf(end) = p(3);

% p_hrf = p;
SPM.xBF.bf = spm_hrf(SPM.xY.RT/16, p_hrf, 16);
SPM = fMRI_design_changeBF(SPM);
res = spm_sp('r',SPM.xX.xKXs,Y); % get the residual
% constraints
err = sum(sum(res.^2))/numel(res);

        