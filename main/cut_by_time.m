function [mdek_n, mdek_s, sa_n, sa_s, rtk_n, rtk_s] = cut_by_time(mdek_N, mdek_S, SA_N, SA_S, RTK_N, RTK_S, time_interval)
    nums = intersect(find([mdek_N.sec_wk] > time_interval(1)),find([mdek_N.sec_wk] < time_interval(2)));
    mdek_n = mdek_N(nums);
    
    nums = intersect(find([mdek_S.sec_wk] > time_interval(1)),find([mdek_S.sec_wk] < time_interval(2)));
    mdek_s = mdek_S(nums);
    
    nums = intersect(find(SA_N(1,:) > time_interval(1)),find(SA_N(1,:) < time_interval(2)));
    sa_n = SA_N(:,nums);
    
    nums = intersect(find(SA_S(1,:) > time_interval(1)),find(SA_S(1,:) < time_interval(2)));
    sa_s = SA_S(:,nums);
    
    nums = intersect(find(RTK_N(1,:) > time_interval(1)),find(RTK_N(1,:) < time_interval(2)));
    rtk_n = RTK_N(:,nums);
    
    nums = intersect(find(RTK_S(1,:) > time_interval(1)),find(RTK_S(1,:) < time_interval(2)));
    rtk_s = RTK_S(:,nums);
end

