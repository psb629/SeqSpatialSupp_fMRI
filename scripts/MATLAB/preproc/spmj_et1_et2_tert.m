function [et1, et2, tert] = spmj_et1_et2_tert(dataDir, subj_name)
    
    % vararginoptions(varargin,{'sn'});
    
    fmapDir = fullfile(dataDir, "BIDS", ['sub-' subj_name], "fmap");
    funcDir = fullfile(dataDir, "BIDS", ['sub-' subj_name], "func");
    
    phasediff = jsondecode(fileread(fullfile(fmapDir, ['sub-' subj_name '_phasediff.json'])));
    func1st = jsondecode(fileread(fullfile(funcDir, ['sub-' subj_name '_task-task_run-01_bold.json'])));
    
    %% compute tert
    % total EPI readout time = = echo spacing (in ms) * base resolution 
    % (also knows as number of echos). If you use GRAPPA acceleration, 
    % you need to divide the total number of echos by two:
    base_resolution = func1st.BaseResolution;
    echo_spacing = func1st.EffectiveEchoSpacing * 1000;                        % compute echo spacing in milliseconds
    tert = echo_spacing * base_resolution;
    
    %% retrieve et1 et2 (in milliseconds)
    et1 = phasediff.EchoTime1 * 1000;
    et2 = phasediff.EchoTime2 * 1000;
