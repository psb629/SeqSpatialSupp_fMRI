function [SPM] = rename_source(SPM)
    if ispc
        dir_raw = 'F:/SeqSpatialSupp_fMRI/imaging_data';
    elseif ismac
        dir_raw = '/Volumes/Diedrichsen_data$/data/SeqSpatialSupp_fMRI/imaging_data';
    end
    idx = strfind(SPM.swd, '/');
    subj = SPM.swd(idx(end)+1:end);
    for ii = 1:length(SPM.xY.VY)
        fname = SPM.xY.VY(ii).fname;
        idx = strfind(fname, ['/',subj]);
        new_fname = strrep(fname, fname(1:idx(1)-1), dir_raw);
        SPM.xY.VY(ii).fname = new_fname;
    end
end