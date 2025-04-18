function [P,SPM,Yhat,Yres]=fit_hrfmulti(SPM,Yraw,R,varargin)
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
%   P(1) - duration of the long                          ? % press
%   P(2) - duration of the short                         ? % press
%   P(3) - duration of the instruction                   ? % instruction
% --------------------------------------------------------------
P0 = [1.8 1.25 1 0.5]';
fit = [];
vararginoptions(varargin,{'P0','fit'}); 

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
P_store = zeros(length(P0),size(hrfs,1));

% loop over all hrfs
for h=1:size(hrfs,1)
    
    fprintf('HRF %d\n',h);
    bf = hrfs(h,:)';
    
    % initial value for the optimization
    th0=P0(fit);
    
    if isempty(fit)
        err(h) = cost(log(th0),Y,SPM,bf,fit,P0,R);
        th = [];
    else
        [th,err(h)]=fminsearch(@cost,log(th0),option,Y,SPM,bf,fit,P0,R); 
    end
    
    P = P0;
    P(fit) = exp(th);
    fprintf('long = %.2f and short = %.2f and intruction time = %.2f and onset instruction = %.2f\n',P(1),P(2),P(3),P(4))
    fprintf('err = %.2f\n',err(h))
    
    % store the parameter
    P_store(:,h) = P;
    
end

% find the best hrf for this subject
[~,idx] = min(err);

SPM.xBF.bf = hrfs(idx,:)';
P = P_store(:,idx);
fprintf('long = %.2f and short = %.2f and intruction time = %.2f and onset instruction = %.2f\n',P(1),P(2),P(3),P(4))
fprintf('err = %.2f\n',min(err))
SPM = fMRI_design_changeBF_multi(SPM,'dur1',P(1),'dur2',P(2),'dur3',P(3),'ons3',P(4),'R',R); 

% return predicted timeseries and residuals 
beta  = SPM.xX.pKX*Y; %-Parameter estimates
Yres  = spm_sp('r',SPM.xX.xKXs,Y); % get the 
reg_interest=[SPM.xX.iH SPM.xX.iC]; 
Yhat   = SPM.xX.xKXs.X(:,reg_interest)*beta(reg_interest,:); %- predicted values 
%Yhat   = SPM.xX.xKXs.X*beta; %- predicted values

function err=cost(th,Y,SPM,bf,fit,P0,R) % Returns accumulated fitting error 

P = P0; 
P(fit) = exp(th);

SPM.xBF.bf = bf;

SPM = fMRI_design_changeBF_multi(SPM,'dur1',P(1),'dur2',P(2),'dur3',P(3),'ons3',P(4),'R',R); 
res = spm_sp('r',SPM.xX.xKXs,Y); % get the residual
err = sum(sum(res.^2))/numel(res);
% constraints
% if P(3)>5
%     err = inf;
% else
%     err = sum(sum(res.^2))/numel(res);
% end


        