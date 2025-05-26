function R = get_behav(subj_id)
%%%%%%%%%%%%%%%%%%%%%
%%%     Input     %%%
%%%%%%%%%%%%%%%%%%%%%
%   subj_id: str
%       Subject ID
%       e.g.) "S01", ..., "R14"
%
%%%%%%%%%%%%%%%%%%%%%
%%%     output    %%%
%%%%%%%%%%%%%%%%%%%%%
%   R: cell
%       Data field

%% Note
% Introduction and explanation of the index,
% 1: seqType (0: Letter cue, 1: Spatial cue)
cueID = ["letter", "spatial"]; R.cueID = cueID;
% 2: seqID (0: 32451, 1:35124, 2:13254, 3:14523)
seqID = [32451, 35124, 13254, 14523]; R.seqID = seqID;
% 3: TrialState (seqID, seqType)=0:(0,0),1:(0,1),2:(1,0),3:(1,1),4:(2,0),5:(2,1),6:(3,0),7:(3,1)
% 4: TransitionState 0:(0,0)->(0,0), 1:(0,0)->(0,1), ... , 62:(3,1)->(3,0), 63:(3,1)->(3,1)
% 5: isError (0: success, 1: error)

%% Load the behavioural data
dir_behav = '/Volumes/Diedrichsen_data$/data/SeqSpatialSupp_fMRI/behavDir';
behav_data = fullfile(dir_behav, sprintf('sub-%s/ssh__%s.dat',subj_id,subj_id));
addpath('/Volumes/Diedrichsen_data$/matlab/dataframe')
S = dload(behav_data); % https://github.com/DiedrichsenLab/dataframe.git/util/dload.m

%% The number of trials/runs
%nRun = length(unique(S.BN)); % 8, 6 for S04
%nTR = length(unique(S.TN)); % 68
nRuns = 8; R.nRuns = nRuns;
nTrials = 68; R.nTrials = nTrials;

%% Check trials
for r=1:nRuns
%    dat_format = 'sub-S%.2d/ssh__S%.2d_%.2d.mov';
%    mov = movload(fullfile(workdir, behavDir, sprintf(data_format, s, s, r)));
    idx = find(S.BN==r); % find the row corresponding to the r-th run.
    tS(r,S.cueP(idx) == seqID(1) & S.seqType(idx)==0) = 0; % 0: (seq,cue)=(0,0)
    tS(r,S.cueP(idx) == seqID(1) & S.seqType(idx)==1) = 1; % 1: (seq,cue)=(0,1)
    tS(r,S.cueP(idx) == seqID(2) & S.seqType(idx)==0) = 2; % 2: (seq,cue)=(1,0)
    tS(r,S.cueP(idx) == seqID(2) & S.seqType(idx)==1) = 3; % 3: (seq,cue)=(1,1)
    tS(r,S.cueP(idx) == seqID(3) & S.seqType(idx)==0) = 4; % 4: (seq,cue)=(2,0)
    tS(r,S.cueP(idx) == seqID(3) & S.seqType(idx)==1) = 5; % 5: (seq,cue)=(2,1)
    tS(r,S.cueP(idx) == seqID(4) & S.seqType(idx)==0) = 6; % 6: (seq,cue)=(3,0)
    tS(r,S.cueP(idx) == seqID(4) & S.seqType(idx)==1) = 7; % 7: (seq,cue)=(3,1)

    transS(r,1) = -1; % The 1st trial doesn't have an index for the transition state.

    R.isError(r,:) = S.isError(idx); % The sequence input was incorrect.

    for t=1:nTrials
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

        if t>1
            % The index of transition state for idx_i -> idx_j: 8*i + j
            i = tS(r,t-1);
            j = tS(r,t); 
            transS(r,t) = 8*i + j;
        end
    end
    R.isValid(r,:) = (S.RT(idx)~=0);
    R.MT(r,:) = S.MT(nTrials*(r-1)+1:nTrials*r); % movement time
    R.RT(r,:) = S.RT(nTrials*(r-1)+1:nTrials*r); % reaction time
end

R.TrialState = tS;
R.TransitionState = transS;

R.idxSeq = zeros(nRuns, nTrials);
R.idxSeq(tS==0 | tS==1) = 0;
R.idxSeq(tS==2 | tS==3) = 1;
R.idxSeq(tS==4 | tS==5) = 2;
R.idxSeq(tS==6 | tS==7) = 3;

R.idxCue = zeros(nRuns, nTrials);
R.idxCue(tS==0 | tS==2 | tS==4 | tS==6) = 0;
R.idxCue(tS==1 | tS==3 | tS==5 | tS==7) = 1;