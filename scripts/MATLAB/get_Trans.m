function [idx_Trans, Trans] = get_Trans(TS_i, TS_f)
    cues = ["L", "S"]; % Letter / Spatial cue
    sequences = [1, 2, 3, 4]; % 1: 32451, 2:35124, 3:13254, 4:14523

    %% Transition Index: TS_i -> TS_f
    idx_Trans = 8*(TS_i-1) + TS_f;

    %% Category: BothRep, SeqRep, CueRep, NonRep / Letter, Spatial (=Cue_f) 
    list_Trans = [];
    for si = 1:length(sequences)
        for ci = 1:length(cues)
            for sf = 1:length(sequences)
                for cf = 1:length(cues)
                    cue = cues(cf);
                    if si ~= sf
                        if ci ~= cf
                            ts = sprintf("N_%s",cue);
                        else
                            ts = sprintf("C_%s",cue);
                        end
                    else
                        if ci ~= cf
                            ts = sprintf("S_%s",cue);
                        else
                            ts = sprintf("B_%s",cue);
                        end
                    end
                    list_Trans = [list_Trans; ts];
                end
            end
        end
    end
    Trans = list_Trans(idx_Trans);

end