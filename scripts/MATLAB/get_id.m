function [subj_id, S_id] = get_id(fname, sn, combine)
    if nargin < 3
        combine = false;
    end
    pinfo = dload(fname);
    subj_id = char(pinfo.subj_id(pinfo.sn==sn));
    S_id = strrep(subj_id,'R','S');
    if combine
        subj_id = subj_id(2:end);
    end
end