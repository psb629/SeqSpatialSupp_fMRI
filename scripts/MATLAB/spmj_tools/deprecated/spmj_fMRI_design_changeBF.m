function SPM=spmj_fMRI_design_changeBF(SPM,dur);
% Re-convolves the SPM structure with a new basis function
% If second input argument is given, also changes duration of the event 


Xx    = [];
Xb    = [];
iCs   = [];     % number of the run 
iCc   = [];     % number of the condition
iCb   = [];     % number of the basis function 
iN    = [];     % Index for regressor of non-interest 

% Check if we have to change duration 
if (nargin>1 && ~isempty(dur) && dur~=SPM.Sess(1).U(1).dur(1))
    for s = 1:length(SPM.nscan)
        for u=1:length(SPM.Sess(s).U)
            SPM.Sess(s).U(u).dur=ones(size(SPM.Sess(s).U(u).dur))*dur; 
        end'/
        SPM.Sess(s).U=spm_get_ons(SPM,s);
    end;
end; 


numbasis=size(SPM.xBF.bf,2); 
numscan=length(SPM.nscan); 

for s = 1:numscan
    
    % number of scans for this session
    %----------------------------------------------------------------------
    k = SPM.nscan(s);
    
    % create convolved stimulus functions or inputs
    %==================================================================
    
    % Get inputs, neuronal causes or stimulus functions U
    %------------------------------------------------------------------
    U=SPM.Sess(s).U;
    numcond=length(U); 
   
    % Convolve stimulus functions with basis functions
    %------------------------------------------------------------------
    [X,Xn,Fc] = spm_Volterra(U,SPM.xBF.bf,SPM.xBF.Volterra);
    
    % Resample regressors at acquisition times (32 bin offset)
    %------------------------------------------------------------------
        X = X((0:(k - 1))*SPM.xBF.T + SPM.xBF.T0 + 32,:);
    
    % and orthogonalise (within trial type)
    %------------------------------------------------------------------
    for i = 1:length(Fc)
        X(:,Fc(i).i) = spm_orth(X(:,Fc(i).i));
    end
    
    
    % get user specified regressors
    %==================================================================
    C     = SPM.Sess(s).C.C;
    numreg = size(C,2);
    X      = [X spm_detrend(C)];
    
    
    %-Session structure array
    %----------------------------------------------------------------------
    SPM.Sess(s).row    = size(Xx,1) + (1:k);
    SPM.Sess(s).col    = size(Xx,2) + (1:size(X,2));
       
    % Confounds: Session effects
    %==================================================================
    B      = ones(k,1);
    
    % append into Xx and Xb
    %======================================================================
    Xx    = blkdiag(Xx,X);
    Xb    = blkdiag(Xb,B);
    
    iCs   = [iCs ones(1,size(X,2)+numreg)*s];
    iCc   = [iCc kron([1:numcond],ones(1,numbasis)) zeros(1,numreg)]; 
    iCb   = [iCb kron(ones(1,numcond),[1:numbasis]) zeros(1,numreg)];
    iN    = [iN zeros(1,numcond*numbasis) ones(1,numreg)];
end

% finished
%--------------------------------------------------------------------------
SPM.xX.X      = [Xx Xb];
if isfield(SPM.xX,'W')
    SPM.xX.xKXs   = spm_sp('Set',spm_filter(SPM.xX.K,SPM.xX.W*SPM.xX.X));       % KWX
else
    SPM.xX.xKXs   = spm_sp('Set',spm_filter(SPM.xX.K,SPM.xX.X));       % KWX
end;
SPM.xX.xKXs.X = full(SPM.xX.xKXs.X);
SPM.xX.pKX    = spm_sp('x-',SPM.xX.xKXs);                        % projector

SPM.xX.iC     = 1:size(Xx,2);
SPM.xX.iB     = (1:size(Xb,2)) + size(Xx,2);
% Indices 
SPM.xX.iCs    = [iCs [1:numscan]];   % number of the session
SPM.xX.iCc    = [iCc zeros(1,numscan)];   % number of the condition
SPM.xX.iCb    = [iCb zeros(1,numscan)];  % number of the basis function
SPM.xX.iN     = [iN  ones(1,numscan)*2];   % Flag for the regressors of no interest (+ intercept)
