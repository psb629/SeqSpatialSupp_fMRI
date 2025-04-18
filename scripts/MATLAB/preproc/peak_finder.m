function [start,peak] = peak_finder(f,npeaks)

prec_speed = 0.05;

x = 1:length(f);
d_f = diff(f);

[pks,loc] = findpeaks(f,x,'SortStr','descend','NPeaks',npeaks,'MinPeakDistance',200);
n_pks = length(pks);

start = zeros(1,n_pks);
peak = zeros(1,n_pks);

loc_sort = sort(loc,'ascend');

for i = 1:n_pks
    
    % find the max speed index
    mask_start_idx = 1;
    mask = 1:loc_sort(1);
    if i==2
        mask_start_idx = loc_sort(1);
        mask = loc_sort(1):loc_sort(2);
    end
    max_speed_idx = find(d_f==max(d_f(mask)),1);
    max_speed = max(d_f(mask));

    % find start press index
    mask = mask_start_idx:max_speed_idx;
    th_speed = max_speed * prec_speed;
    start_idx = (mask_start_idx-1) + find(d_f(mask)<th_speed);

    if i==1
        if isempty(start_idx)
            warning('this should not happen')
            start_idx = 1;
        end
    end
    
    start(i) = start_idx(end);
    peak(i) = loc_sort(i);

end


    


