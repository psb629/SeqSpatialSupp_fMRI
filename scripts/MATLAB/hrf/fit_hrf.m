function [SPM,Yhat,Yres]=fit_hrf(SPM,Yraw)
% Fits the two-gamma hrf fucntion to extracted time 
% function [hrf,P,xBF]=spmj_fit_hrf(SPM,Y_raw,options);
% INPUT: 
%   SPM:    SPM-structure with model to be fitted 
%   Yraw:   raw time series for hrf fitting 
% OUTPUT: 
%   P:      parameters fitted 
%   SPM:    changed SPM
%   Yhat:   predicted time series 
%   Yres:   residual time series 
%   P(1) - duration of the event type1                   ? % press
%   P(2) - duration of the event type2                   ? % instruction
% --------------------------------------------------------------
% P0 = [18.75 0.1]';
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

% load the HRF groups - these are the basis function that are going to be
% tested
hrfs = getcanonicalhrflibrary(0,SPM.xBF.dt);

err = zeros(1,size(hrfs,1));
% P_store = zeros(length(P0),size(hrfs,1));

% loop over all hrfs
for h=1:size(hrfs,1)
    
    fprintf('HRF %d\n', h);
    bf = hrfs(h,:)';
    
    % initial value for the optimization
    %th0=P0(fit);
    
%     if isempty(fit)
    err(h) = cost(Y,SPM,bf);
    %  th = [];
%     else
%         [th, err(h)]=fminsearch(@cost,log(th0),option,Y,SPM,bf,fit,P0); 
%     end
    
    %P = P0;
    %P(fit) = exp(th);
    %fprintf('press time = %.2f and intruction time = %.2f\n', P(1), P(2))
    
    % store the parameter
    %P_store(:,h) = P;
    
end

% find the best hrf for this subject
[~,idx] = min(err);

SPM.xBF.bf = hrfs(idx,:)';
%P = P_store(:,idx);
SPM = fMRI_design_changeBF(SPM); 

% return predicted timeseries and residuals 
beta  = SPM.xX.pKX*Y; %-Parameter estimates
Yres  = spm_sp('r',SPM.xX.xKXs,Y); % get the 
reg_interest=[SPM.xX.iH SPM.xX.iC]; 
Yhat   = SPM.xX.xKXs.X(:,reg_interest)*beta(reg_interest,:); %- predicted values 
%Yhat   = SPM.xX.xKXs.X*beta; %- predicted values

function err=cost(Y,SPM,bf)
% function err=cost(th,Y,SPM,bf,fit,P0) % Returns accumulated fitting error 
%P = P0; 
%P(fit) = exp(th);

SPM.xBF.bf = bf;

SPM = fMRI_design_changeBF(SPM);
res = spm_sp('r',SPM.xX.xKXs,Y); % get the residual

% constraints
err = sum(sum(res.^2))/numel(res);


        