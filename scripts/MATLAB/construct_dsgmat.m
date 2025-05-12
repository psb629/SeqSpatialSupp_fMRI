function R = construct_dsgmat(behav_data, glm)
%%%%%%%%%%%%%%%%%%%%%
%%%     Input     %%%
%%%%%%%%%%%%%%%%%%%%%
%   behav_data: str
%       The file name.
%       e.g.) F:/SeqSpatialSupp_fMRI/behavDir/sub-S01/ssh__S01.dat
%%%%%%%%%%%%%%%%%%%%%
%%%     output    %%%
%%%%%%%%%%%%%%%%%%%%%
%   R: cell
%       Data field

% condition number:
% 1: seqType (0: Number cue, 1: Spatial visual cue
% 2: seqID (0: 32451, 1:35124, 2:13254, 3:14523)
% 3: trialState ((seqID, seqType)=0:(0,0),1:(0,1),2:(1,0),3:(1,1),4:(2,0),5:(2,1),6:(3,0),7:(3,1)
% 4: transitionState
% 5: isError (0: success, 1: error)

seqID = [32451, 35124, 13254, 14523];
nTr = 68;
% transition ID: 1~68, 65-68: non-interest, 1:(0,0)->(0,0)
temp = [1 3 5 7 10 12 14 16];
repmotorID = [1 2 9 10 19 20 27 28 37 38 45 46 55 56 63 64]; %% transition ID with repated motor condition
rep1stfingID = [1:4 9:12 17:20 25:28 37:40 45:48 53:56 61:64];%% transition ID with repated first finger 
repcueID  = [temp temp+16 temp+32 temp+48]; %% transition ID with repeated cue condition
repbothID = 9*[0 1 2 3 4 5 6 7]+1; %% 0: (0,0)->(0,0), (0,1)->(0,1),...
nintID = [65:68]; %% transition ID with non-interest
nrepID = find(~ismember([1:nTr], [repcueID repmotorID nintID]));

behavDir        = 'behavDir';           % Timing data from the scanner
S = dload(behav_data);
%nRun = length(unique(S.BN)); % 8, 6 for S04
%nTR = length(unique(S.TN)); % 68
nRun = 8;
nTR = 68;
for r=1:nRun
%    dat_format = 'sub-S%.2d/ssh__S%.2d_%.2d.mov';
%    mov = movload(fullfile(workdir, behavDir, sprintf(data_format, s, s, r)));
    idx = find(S.BN==r);
    tS(r,S.cueP(idx) == seqID(1) & S.seqType(idx)==0) = 0;
    tS(r,S.cueP(idx) == seqID(1) & S.seqType(idx)==1) = 1;
    tS(r,S.cueP(idx) == seqID(2) & S.seqType(idx)==0) = 2;
    tS(r,S.cueP(idx) == seqID(2) & S.seqType(idx)==1) = 3;
    tS(r,S.cueP(idx) == seqID(3) & S.seqType(idx)==0) = 4;
    tS(r,S.cueP(idx) == seqID(3) & S.seqType(idx)==1) = 5;
    tS(r,S.cueP(idx) == seqID(4) & S.seqType(idx)==0) = 6;
    tS(r,S.cueP(idx) == seqID(4) & S.seqType(idx)==1) = 7;
    R.isError(r,:) = S.isError(idx);

    for t=1:nTR
        % R.onset(r,t) = (S.startTimeReal(idx(t)) + 1000 + S.RT(idx(t)))*0.001;
        R.onset(r,t) = (S.startTimeReal(idx(t)) + 1000) * 0.001;

        if S.RT(idx(t))==0
            R.dur(r,t)=0;  %% invalid trials 
        else
            R.dur(r,t) = 2; % fix the duration
            % if S.MT(idx(t))==0   
            %     R.dur(r,t) = (3000-S.MT(idx(t)))*0.001;
            % else
            %     R.dur(r,t) = S.MT(idx(t))*0.001;
            % end
        end
    end
    R.isValid(r,:) = (S.RT(idx)~=0);
    R.isValidRep(r,:) = R.isValid(r,:);
    R.isValidRep(r,min(find(R.isValid(r,:)==0)+1,nTR))=0;
    R.isValidRepmotor(r,:) = R.isValidRep(r,:);
    R.isValidRepmotor(r,min(find(R.isError(r,:)==1)+1,nTR))=0;
end

R.tS = tS;

% if glm ==1
trial_idx = [2:17 19:34 36:51 53:68]; % trial numbers for repetition suppression effects + 1
dummy_idx = [1 18 35 52];
for r=1:nRun
    for i=1:length(trial_idx)
        % R.cond: 1~68, transition: (seqID, seqType), (0,0)->(0,0):1, (0,0)->(0,1):2,
        % (0,1)->(0,0):3, (0,1)->(0,1):4,...(3,1)->(3,1):64, 
        % 1st trial:65, 18th: 66, 35th: 67, 52nd: 68 
        R.transID(r,trial_idx(i)) = tS(r,trial_idx(i)-1)*8+tS(r,trial_idx(i))+1; 
    end
    R.transID(r,dummy_idx) = [65:68];
    % R.cond(r,dummy_idx) = [65:68];  % conditions of non-interest
    R.isRepMotor(r,:) = ismember(R.transID(r,:), repmotorID) & R.isValidRepmotor(r,:);
    R.isNRepMotor(r,:) = ~ismember(R.transID(r,:), repmotorID) & R.isValidRepmotor(r,:) & ~ismember([1:nTR],dummy_idx);
    R.isRepCue(r,:) = ismember(R.transID(r,:), repcueID) & R.isValidRep(r,:);
    R.isNRepCue(r,:) = ~ismember(R.transID(r,:), repcueID) & R.isValidRep(r,:) & ~ismember([1:nTR],dummy_idx);
    
    R.isRepBoth(r,:) = ismember(R.transID(r,:), repbothID) & R.isValidRepmotor(r,:);
    R.isNRep(r,:) = ismember(R.transID(r,:), nrepID) & R.isValidRep(r,:);
    % R.isNint(r,:) = ismember(R.transID(r,:), nintID);
    R.isNint(r,:) = ~ismember([1:nTr], find(R.isRepMotor(r,:)+R.isRepCue(r,:)+R.isNRep(r,:)~=0));

    R.MT(r,:) = S.MT(nTR*(r-1)+1:nTR*r);
    R.RT(r,:) = S.RT(nTR*(r-1)+1:nTR*r);
    
end
R.isLetter = zeros(nRun, nTR);
R.isSpatial = zeros(nRun, nTR);
R.isLetter(tS==0 | tS==2 | tS==4 | tS==6) =1;
R.isSpatial(tS==1 | tS==3 | tS==5 | tS==7) =1;
R.isRepBothLetter = R.isRepBoth.*R.isLetter;
R.isRepBothSpatial = R.isRepBoth.*R.isSpatial;
R.isRepCueLetter = R.isRepCue.*R.isLetter;
R.isRepCueSpatial = R.isRepCue.*R.isSpatial;
R.repmotor_transID = repmotorID;
R.repcue_transID = repcueID;
R.repboth_transID = repbothID;
R.nrep_transID = nrepID;
R.nint_transID = nintID;
R.cond = zeros(nRun,nTR);
if glm==0
    R.cond = ones(1,nTR); 
    idx = find(S.BN==9);
%     R.onset = (S.startTimeReal(idx) + 1000 + S.RT(idx))/1000;
    R.onset = S.startTimeReal(idx)/1000;
elseif glm==1
    R.cond = repmat([1:nTR],nRun, 1);
elseif glm==2 % motor or cue repetition
    R.cond(R.isRepMotor==1 & R.isRepBoth==0 & R.isLetter==1)=1;
    R.cond(R.isRepMotor==1 & R.isRepBoth==0 & R.isSpatial==1)=2;

    R.cond(R.isRepCue==1 & R.isRepBoth==0 & R.isLetter==1)=3;
    R.cond(R.isRepCue==1 & R.isRepBoth==0 & R.isSpatial==1)=4;

    R.cond(R.isRepBoth & R.isLetter==1)=5;
    R.cond(R.isRepBoth & R.isSpatial==1)=6;
    R.cond(R.isNRep & R.isLetter==1)=7;
    R.cond(R.isNRep & R.isSpatial==1)=8;
    R.cond(R.isNint)=9;
    
elseif glm==3
    R.cond = tS+1;
    R.cond(R.isError==1 | R.isValid==0)=9;
elseif glm==4
    R.cond(R.isLetter==1 & R.isValid==1)=1; %% to be modifed
    R.cond(R.isSpatial==1 & R.isValid==1)=2;
    R.cond(~R.isValid)=3;
end

% end
% else
% if glm ==2 || glm==3 % glm=2: motor repetition, glm=3: cue repetition
%     con = calc_contrasts(R);
%     for r=1:nRun
%         R.cond(r, find(con(glm-1,68*(r-1)+1:68*r)==1))=1;
%         R.cond(r, find(con(glm-1,68*(r-1)+1:68*r)==-1))=2;
%         R.cond(r, find(con(glm-1,68*(r-1)+1:68*r)==0))=3;
% elseif glm==4  %% letter vs spatial visual cues
%     con = calc_contrasts(R);
%     for r=1:nRun
%         R.cond(r, ismember(R.tS(r,:),[0 2 4 6]))=1;  % Number cue
%         R.cond(r, ismember(R.tS(r,:),[1 3 5 7]))=2;  % Spatial visual cue
%     end

% end



        % R.cond = tS;

