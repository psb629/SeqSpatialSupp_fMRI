function [idx_TS, TS] = get_TS(seq, cue)
    %% sequences = [1, 2, 3, 4]; % 1: 32451, 2:35124, 3:13254, 4:14523    
    idx_s = seq;
    %% cues = ["L", "S"]; % Letter / Spatial cue
    if cue == "L"
        idx_c = 1;
    elseif cue == "S"
        idx_c = 2;
    end
    %% return the Trial State (TS)
    idx_TS = (idx_s-1)*2 + idx_c;
    TS = sprintf("(%d,%s)",seq,cue);
end