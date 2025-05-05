function [subj_id, S_id] = get_id(fname, sn)
    pinfo = dload(fname);
    subj_id = char(pinfo.subj_id(pinfo.sn==sn));
    S_id = strrep(subj_id,'R','S');
end