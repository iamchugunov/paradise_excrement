function [mdek_N, mdek_S, RTK_N, RTK_S, SA_N, SA_S, BASE] = open_logs()
    
    load('logs_mat_files\nikita_mdek_full.mat')
    load('logs_mat_files\sereja_mdek_full.mat')

    UTC_epoch_seconds_N = [nikita_mdek_full.PC_time] + 18;
    UTC_offset_N = UTC_epoch_seconds_N/(24*60*60);
    atomTime_N = UTC_offset_N + datenum(1970,1,1);
    [nik_h, nik_min, nik_s] = hms(datetime(atomTime_N,'ConvertFrom','datenum')); 
    [nik_y, nik_m, nik_d] = ymd(datetime(atomTime_N,'ConvertFrom','datenum'));

    for k = 1:length(nik_s)
        [~, nik_sec_wk(k)] = GPSweek(nik_y(k),nik_m(k),nik_d(k),nik_h(k),nik_min(k),nik_s(k));
        nikita_mdek_full(k).sec_wk = nik_sec_wk(k);
    end

    UTC_epoch_seconds_S = [sereja_mdek_full.PC_time] + 18;
    UTC_offset_S = UTC_epoch_seconds_S/(24*60*60);
    atomTime_S = UTC_offset_S + datenum(1970,1,1);
    [ser_h, ser_min, ser_s] = hms(datetime(atomTime_S,'ConvertFrom','datenum')); 
    [ser_y, ser_m, ser_d] = ymd(datetime(atomTime_S,'ConvertFrom','datenum'));

    for k = 1:length(ser_s)
        [~, ser_sec_wk(k)] = GPSweek(ser_y(k),ser_m(k),ser_d(k),ser_h(k),ser_min(k),ser_s(k));
        sereja_mdek_full(k).sec_wk = ser_sec_wk(k);
    end
    
    mdek_N = nikita_mdek_full;
    for i = 1:length(mdek_N)
        mdek_N(i).count = length(find(mdek_N(i).anc_dist));
    end
    mdek_S = sereja_mdek_full;
    for i = 1:length(mdek_S)
        mdek_S(i).count = length(find(mdek_S(i).anc_dist));
    end
    
    load('logs_mat_files\base_rtk.mat')
    BASE = out;
    base_mean_ecef = [mean(BASE(2,:)) mean(BASE(3,:)) mean(BASE(4,:))];
    base_mean_lla = ecef2lla(base_mean_ecef);
    wgs84 = wgs84Ellipsoid;
    
    
    RTK_N = [];
    load('logs_mat_files\nikita_1_rtk.mat')
    RTK_N = [RTK_N out];
    load('logs_mat_files\nikita_2_rtk.mat')
    RTK_N = [RTK_N out];
    load('logs_mat_files\nikita_3_rtk.mat')
    RTK_N = [RTK_N out];
    for i = 1: length(RTK_N)
        [RTK_N(2,i),RTK_N(3,i),RTK_N(4,i)] = ecef2enu(RTK_N(2,i),RTK_N(3,i),RTK_N(4,i),base_mean_lla(1),base_mean_lla(2),base_mean_lla(3), wgs84);
    end
    
    RTK_S = [];
    load('logs_mat_files\sereja_1_rtk.mat')
    RTK_S = [RTK_S out];
    load('logs_mat_files\sereja_2_rtk.mat')
    RTK_S = [RTK_S out];
    load('logs_mat_files\sereja_3_rtk.mat')
    RTK_S = [RTK_S out];
    for i = 1: length(RTK_S)
        [RTK_S(2,i),RTK_S(3,i),RTK_S(4,i)] = ecef2enu(RTK_S(2,i),RTK_S(3,i),RTK_S(4,i),base_mean_lla(1),base_mean_lla(2),base_mean_lla(3), wgs84);
    end
    
    SA_N = [];
    load('logs_mat_files\nikita_1_sa.mat')
    SA_N = [SA_N out];
    load('logs_mat_files\nikita_2_sa.mat')
    SA_N = [SA_N out];
    load('logs_mat_files\nikita_3_sa.mat')
    SA_N = [SA_N out];
    for i = 1: length(SA_N)
        [SA_N(2,i),SA_N(3,i),SA_N(4,i)] = ecef2enu(SA_N(2,i),SA_N(3,i),SA_N(4,i),base_mean_lla(1),base_mean_lla(2),base_mean_lla(3), wgs84);
        [SA_N(5,i),SA_N(6,i),SA_N(7,i)] = ecef2enuv(SA_N(5,i),SA_N(6,i),SA_N(7,i),base_mean_lla(1),base_mean_lla(2));
    end
    SA_N = SA_N([1 2 5 3 6 4 7],:);
    
    SA_S = [];
    load('logs_mat_files\sereja_1_sa.mat')
    SA_S = [SA_S out];
    load('logs_mat_files\sereja_2_sa.mat')
    SA_S = [SA_S out];
    load('logs_mat_files\sereja_3_sa.mat')
    SA_S = [SA_S out];
    for i = 1: length(SA_S)
        [SA_S(2,i),SA_S(3,i),SA_S(4,i)] = ecef2enu(SA_S(2,i),SA_S(3,i),SA_S(4,i),base_mean_lla(1),base_mean_lla(2),base_mean_lla(3), wgs84);
        [SA_S(5,i),SA_S(6,i),SA_S(7,i)] = ecef2enuv(SA_S(5,i),SA_S(6,i),SA_S(7,i),base_mean_lla(1),base_mean_lla(2));
    end
    SA_S = SA_S([1 2 5 3 6 4 7],:);
    
    
    
end

