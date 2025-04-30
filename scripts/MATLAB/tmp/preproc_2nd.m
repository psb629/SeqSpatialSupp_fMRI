% s = 11;

for sn=[28]  % S02, S03, S05, S06, // S08, S11, S13 
    sss_imana('BIDS:move_unzip_raw_func','sn',sn);
    sss_imana('BIDS:move_unzip_raw_fmap','sn',sn);
    
    s = sn-14;
    change_filenames(s);
    src_dir = fullfile('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/imaging_data_raw',sprintf('R%02d',s));
    dst_dir = fullfile('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/imaging_data_raw',sprintf('S%02d',s));
    % Form the shell command
    system(sprintf('mv "%s"/* "%s"', src_dir, dst_dir));

    sss_imana('FUNC:make_fmap','sn',s);
    
    
    sss_imana('FUNC:realign_unwarp','sn',s,'rtm',0,'endTR',[410*ones(1,16) 385]);
    sss_imana('FUNC:move_realigned_images','sn',s,'rtm',0);
    sss_imana('FUNC:meanimage_bias_correction','sn',s,'prefix','u','rtm',0);



end
% 
sss_imana('FUNC:coreg','sn',s);
sss_imana('FUNC:make_samealign','sn',s);
sss_imana('FUNC:make_maskImage','sn',s);