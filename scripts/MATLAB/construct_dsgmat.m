function R = construct_dsgmat(R, glm)
%%%%%%%%%%%%%%%%%%%%%
%%%     Input     %%%
%%%%%%%%%%%%%%%%%%%%%
%   R: struct
%       The data field R created by get_bahav.m
%   glm: GLM number to be perfomred
%
%%%%%%%%%%%%%%%%%%%%%
%%%     output    %%%
%%%%%%%%%%%%%%%%%%%%%
%   R: struct
%       Enhanced data field with GLM regressors
%

% nRuns = R.nRuns;
% nTrials = R.nTrials;

switch glm
    case 1
    %% GLM3 : Trial State
    % 0: (0,0), 1: (0,1), 2: (1,0), 3: (1,1), 4: (2,0), 5: (2,1), 6: (3,0), 7: (3,1)
    
    cond = R.TrialState+1; % The elements of TrialState would be a member of set [1:8].
    cond(R.isError==1|R.isValid==0) = 0; % The 9-th TrialState: Non-interest (Incorrect)

    labelReg = string(R.TrialState);
    labelReg(cond==0) = "N"; % Non-interest

    R.labelreg = labelReg;
    R.cond = cond;
    R.dur = repmat(2, size(cond)); % The durations of all regressors are set to 2 seconds.

    case 2
    %% GLM2 : Repetition Suppression
    % |(s,c)|(0,0)|(0,1)|(1,0)|(1,1)|(2,0)|(2,1)|(3,0)|(3,1)| 
    % |-----|-----|-----|-----|-----|-----|-----|-----|-----|
    % |(0,0)|(B00 | S01 | C02 | N03)| C04 | N05 | C06 | N07 |
    % |(0,1)|(    | B09 | N10 | C11)| N12 | C13 | N14 | C15 |
    % |(1,0)|(    |     | B18 | S19)| C20 | N21 | C22 | N23 |
    % |(1,1)|(    |     |     | B27)| N28 | C29 | N30 | C31 |
    % |(2,0)|     |     |     |     |(B36 | S37 | C38 | N39)|
    % |(2,1)|     |     |     |     |(    | B45 | N46 | C47)|
    % |(3,0)|     |     |     |     |(    |     | B54 | S55)|
    % |(3,1)|     |     |     |     |(    |     |     | B63)|
    % B: Both-Rep, S: Seq-Rep, C: Cue-Rep, N: No-Rep, (First Finger)
    
    idxB = [];
    idxS = [];
    idxC = [];
    idxN = [];
    % regF = []; % the repetition of the 1st finger
    for si=0:3
        for ci=0:1
            TSi = 2*si + ci; % initial Trial State
            for sf=0:3
                for cf=0:1
                    TSf = 2*sf + cf; % final Trial State
                    idx = 8*TSi + TSf; % Transition Index: (si,ci)->(sf,cf)
                    if si==sf % sequences are same
                        if ci==cf % cues are also same
                            idxB = [idxB idx];
                        else % cue are different
                            idxS = [idxS idx];
                        end
                    else % sequences are different
                        if ci==cf % cues are same
                            idxC = [idxC idx];
                        else % cues are different either
                            idxN = [idxN idx];
                        end
                    end
                    % if si<=1 && sf<=1
                    %     regF = [regF idx];
                    % elseif si>=2 && sf>=2
                    %     regF = [regF idx];
                    % end
                end
            end
        end
    end
    
    labelReg = repmat("N",size(R.TransitionState));
    labelReg(ismember(R.TransitionState, idxB)) = "B";
    labelReg(ismember(R.TransitionState, idxS)) = "S";
    labelReg(ismember(R.TransitionState, idxC)) = "C";

    % letter vs. spatial cue? which one? (current or previous)

end

% for r=1:nRuns
%     R.isValidRep(r,:) = R.isValid(r,:);
%     R.isValidRep(r,min(find(R.isValid(r,:)==0)+1,nTrials))=0;
%     R.isValidRepmotor(r,:) = R.isValidRep(r,:);
%     R.isValidRepmotor(r,min(find(R.isError(r,:)==1)+1,nTrials))=0;
% end
% 
% % if glm ==1
% trial_idx = [2:17 19:34 36:51 53:68]; % trial numbers for repetition suppression effects + 1
% dummy_idx = [1 18 35 52];
% for r=1:nRuns
%     for i=1:length(trial_idx)
%         % R.cond: 1~68, transition: (seqID, seqType), (0,0)->(0,0):1, (0,0)->(0,1):2,
%         % (0,1)->(0,0):3, (0,1)->(0,1):4,...(3,1)->(3,1):64, 
%         % 1st trial:65, 18th: 66, 35th: 67, 52nd: 68 
%         R.transID(r,trial_idx(i)) = tS(r,trial_idx(i)-1)*8+tS(r,trial_idx(i))+1; 
%     end
%     R.transID(r,dummy_idx) = [65:68];
%     % R.cond(r,dummy_idx) = [65:68];  % conditions of non-interest
%     R.isRepMotor(r,:) = ismember(R.transID(r,:), repmotorID) & R.isValidRepmotor(r,:);
%     R.isNRepMotor(r,:) = ~ismember(R.transID(r,:), repmotorID) & R.isValidRepmotor(r,:) & ~ismember([1:nTrials],dummy_idx);
%     R.isRepCue(r,:) = ismember(R.transID(r,:), repcueID) & R.isValidRep(r,:);
%     R.isNRepCue(r,:) = ~ismember(R.transID(r,:), repcueID) & R.isValidRep(r,:) & ~ismember([1:nTrials],dummy_idx);
% 
%     R.isRepBoth(r,:) = ismember(R.transID(r,:), repbothID) & R.isValidRepmotor(r,:);
%     R.isNRep(r,:) = ismember(R.transID(r,:), nrepID) & R.isValidRep(r,:);
%     % R.isNint(r,:) = ismember(R.transID(r,:), nintID);
%     R.isNint(r,:) = ~ismember([1:nTr], find(R.isRepMotor(r,:)+R.isRepCue(r,:)+R.isNRep(r,:)~=0));
% 
% end
% 
% R.isRepBothLetter = R.isRepBoth.*R.isLetter;
% R.isRepBothSpatial = R.isRepBoth.*R.isSpatial;
% R.isRepCueLetter = R.isRepCue.*R.isLetter;
% R.isRepCueSpatial = R.isRepCue.*R.isSpatial;
% R.repmotor_transID = repmotorID;
% R.repcue_transID = repcueID;
% R.repboth_transID = repbothID;
% R.nrep_transID = nrepID;
% R.nint_transID = nintID;
% R.cond = zeros(nRun,nTrials);
% if glm==0
%     R.cond = ones(1,nTrials); 
%     idx = find(S.BN==9);
% %     R.onset = (S.startTimeReal(idx) + 1000 + S.RT(idx))/1000;
%     R.onset = S.startTimeReal(idx)/1000;
% elseif glm==1
%     R.cond = repmat([1:nTrials],nRun, 1);
% elseif glm==2 % motor or cue repetition
%     R.cond(R.isRepMotor==1 & R.isRepBoth==0 & R.isLetter==1)=1;
%     R.cond(R.isRepMotor==1 & R.isRepBoth==0 & R.isSpatial==1)=2;
% 
%     R.cond(R.isRepCue==1 & R.isRepBoth==0 & R.isLetter==1)=3;
%     R.cond(R.isRepCue==1 & R.isRepBoth==0 & R.isSpatial==1)=4;
% 
%     R.cond(R.isRepBoth & R.isLetter==1)=5;
%     R.cond(R.isRepBoth & R.isSpatial==1)=6;
%     R.cond(R.isNRep & R.isLetter==1)=7;
%     R.cond(R.isNRep & R.isSpatial==1)=8;
%     R.cond(R.isNint)=9;
% 
% elseif glm==3
%     R.cond = tS+1;
%     R.cond(R.isError==1 | R.isValid==0)=9;
% elseif glm==4
%     R.cond(R.isLetter==1 & R.isValid==1)=1; %% to be modifed
%     R.cond(R.isSpatial==1 & R.isValid==1)=2;
%     R.cond(~R.isValid)=3;
% end
