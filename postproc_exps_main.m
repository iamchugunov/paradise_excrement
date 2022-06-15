clear all
close all
clc

load('logs_mat_files\nikita_2_rtk.mat')

GPS_sec_wk_N = GPS_sec_wk;
GPS_wk_N = GPS_wk;
ECEF_coords_N = out;

clearvars -except GPS_sec_wk_N GPS_wk_N ECEF_coords_N

load('logs_mat_files\sereja_2_rtk.mat')

GPS_sec_wk_S = GPS_sec_wk;
GPS_wk_S = GPS_wk;
ECEF_coords_S = out;

clearvars -except GPS_sec_wk_N GPS_wk_N ECEF_coords_N GPS_sec_wk_S GPS_wk_S ECEF_coords_S

load('logs_mat_files\base_rtk.mat')

GPS_sec_wk_B = GPS_sec_wk;
GPS_wk_B = GPS_wk;
ECEF_coords_B = out;

clearvars -except GPS_sec_wk_N GPS_wk_N ECEF_coords_N GPS_sec_wk_S GPS_wk_S ECEF_coords_S GPS_sec_wk_B GPS_wk_B ECEF_coords_B

load('logs_mat_files\nikita_2_sa.mat')

GPS_sec_wk_N_sa = GPS_sec_wk;
GPS_wk_N_sa = GPS_wk;
ECEF_coords_N_sa = out;

clearvars -except GPS_sec_wk_N GPS_wk_N ECEF_coords_N GPS_sec_wk_S GPS_wk_S ECEF_coords_S GPS_sec_wk_B GPS_wk_B ECEF_coords_B GPS_sec_wk_N_sa GPS_wk_N_sa ECEF_coords_N_sa

load('logs_mat_files\sereja_2_sa.mat')

GPS_sec_wk_S_sa = GPS_sec_wk;
GPS_wk_S_sa = GPS_wk;
ECEF_coords_S_sa = out;

clearvars -except GPS_sec_wk_N GPS_wk_N ECEF_coords_N GPS_sec_wk_S GPS_wk_S ECEF_coords_S GPS_sec_wk_B GPS_wk_B ECEF_coords_B GPS_sec_wk_N_sa GPS_wk_N_sa ECEF_coords_N_sa GPS_sec_wk_S_sa GPS_wk_S_sa ECEF_coords_S_sa

load('logs_mat_files\nikita_mdek_full.mat')
load('logs_mat_files\sereja_mdek_full.mat')

UTC_epoch_seconds_N = [nikita_mdek_full.PC_time] + 18;
UTC_offset_N = UTC_epoch_seconds_N/(24*60*60);
atomTime_N = UTC_offset_N + datenum(1970,1,1);
[nik_h, nik_min, nik_s] = hms(datetime(atomTime_N,'ConvertFrom','datenum')); 
[nik_y, nik_m, nik_d] = ymd(datetime(atomTime_N,'ConvertFrom','datenum'));

for k = 1:length(nik_s)
    [~, nik_sec_wk(k)] = GPSweek(nik_y(k),nik_m(k),nik_d(k),nik_h(k),nik_min(k),nik_s(k));
end

UTC_epoch_seconds_S = [sereja_mdek_full.PC_time] + 18;
UTC_offset_S = UTC_epoch_seconds_S/(24*60*60);
atomTime_S = UTC_offset_S + datenum(1970,1,1);
[ser_h, ser_min, ser_s] = hms(datetime(atomTime_S,'ConvertFrom','datenum')); 
[ser_y, ser_m, ser_d] = ymd(datetime(atomTime_S,'ConvertFrom','datenum'));

for k = 1:length(ser_s)
    [~, ser_sec_wk(k)] = GPSweek(ser_y(k),ser_m(k),ser_d(k),ser_h(k),ser_min(k),ser_s(k));
end

%%%%%%%%%%%%%% Обрезка GPS логов
% GPS_sec_start_dif = ECEF_coords_N_sa(1,1) - ECEF_coords_S_sa(1,1); 
% 
% if GPS_sec_start_dif < 0
%     [~,nik_GPS_start_ind] = min(abs(ECEF_coords_N_sa(1,:) - ECEF_coords_S_sa(1,1)));
%     ECEF_coords_N_sa = ECEF_coords_N_sa(:,nik_GPS_start_ind:end);
% else
%     [~,ser_GPS_start_ind] = min(abs(ECEF_coords_S_sa(1,:) - ECEF_coords_N_sa(1,1)));
%     ECEF_coords_S_sa = ECEF_coords_S_sa(:,ser_GPS_start_ind:end);
% end

% nik_gps_start_ind = find(round(ECEF_coords_N_sa(1,:),1) == 38350);
% nik_gps_end_ind = find(round(ECEF_coords_N_sa(1,:),1) == 38500);
% ser_gps_start_ind = find(round(ECEF_coords_S_sa(1,:),1) == 38350);
% ser_gps_end_ind = find(round(ECEF_coords_S_sa(1,:),1) == 38500);

% ECEF_coords_N_sa = ECEF_coords_N_sa(:,nik_gps_start_ind:nik_gps_end_ind);
% ECEF_coords_S_sa = ECEF_coords_S_sa(:,ser_gps_start_ind:ser_gps_end_ind);
%%%%%%%%%%%%%%%%%%%% лог 1
% nik_gps_start_ind = find(round(ECEF_coords_N_sa(1,:),1) == 39135);
% nik_gps_end_ind = find(round(ECEF_coords_N_sa(1,:),1) == 39425);
% ser_gps_start_ind = find(round(ECEF_coords_S_sa(1,:),1) == 39135);
% ser_gps_end_ind = find(round(ECEF_coords_S_sa(1,:),1) == 39425);
% 
% nik_gps_rtk_start_ind = find(round(ECEF_coords_N(1,:),1) == 39135);
% nik_gps_rtk_end_ind = find(round(ECEF_coords_N(1,:),1) == 39425);
% ser_gps_rtk_start_ind = find(round(ECEF_coords_S(1,:),1) == 39135);
% ser_gps_rtk_end_ind = find(round(ECEF_coords_S(1,:),1) == 39425);
%%%%%%%%%%%%%%%%%%%%% лог 2
nik_gps_start_ind = find(round(ECEF_coords_N_sa(1,:),1) == 39820);
nik_gps_end_ind = find(round(ECEF_coords_N_sa(1,:),1) == 40100);
ser_gps_start_ind = find(round(ECEF_coords_S_sa(1,:),1) == 39820);
ser_gps_end_ind = find(round(ECEF_coords_S_sa(1,:),1) == 40100);

nik_gps_rtk_start_ind = find(round(ECEF_coords_N(1,:),1) == 39820);
nik_gps_rtk_end_ind = find(round(ECEF_coords_N(1,:),1) == 40100);
ser_gps_rtk_start_ind = find(round(ECEF_coords_S(1,:),1) == 39820);
ser_gps_rtk_end_ind = find(round(ECEF_coords_S(1,:),1) == 40100);
%%%%%%%%%%%%%%%%%%%%%
ECEF_coords_N_sa = ECEF_coords_N_sa(:,nik_gps_start_ind:nik_gps_end_ind);
ECEF_coords_S_sa = ECEF_coords_S_sa(:,ser_gps_start_ind:ser_gps_end_ind);

ECEF_coords_N = ECEF_coords_N(:,nik_gps_rtk_start_ind:nik_gps_rtk_end_ind);
ECEF_coords_S = ECEF_coords_S(:,ser_gps_rtk_start_ind:ser_gps_rtk_end_ind);

%%%%%%%%%%%%%%

[~,mdek_start_ind_N] = min(abs(nik_sec_wk - ECEF_coords_N_sa(1,1)));
[~,mdek_end_ind_N] = min(abs(nik_sec_wk - ECEF_coords_N_sa(1,end)));

nik_sec_wk = nik_sec_wk(mdek_start_ind_N:mdek_end_ind_N);
nikita_mdek_full = nikita_mdek_full(mdek_start_ind_N:mdek_end_ind_N);

[~,mdek_start_ind_S] = min(abs(ser_sec_wk - ECEF_coords_S_sa(1,1)));
[~,mdek_end_ind_S] = min(abs(ser_sec_wk - ECEF_coords_S_sa(1,end)));

ser_sec_wk = ser_sec_wk(mdek_start_ind_S:mdek_end_ind_S);
sereja_mdek_full = sereja_mdek_full(mdek_start_ind_S:mdek_end_ind_S);


flag_1meas_N = [];
flag_2meas_N = [];
flag_3meas_N = [];
flag_4meas_N = [];

for k = 1:length(nikita_mdek_full)
    cur_pos = nikita_mdek_full(k).anc_dist; 
    NoM = length(find(cur_pos));
    switch NoM
        case 1
             flag_1meas_N = [flag_1meas_N nik_sec_wk(k)];
        case 2
             flag_2meas_N = [flag_2meas_N nik_sec_wk(k)];
        case 3
             flag_3meas_N = [flag_3meas_N nik_sec_wk(k)];
        case 4
             flag_4meas_N = [flag_4meas_N nik_sec_wk(k)];
    end
end

flag_1meas_S = [];
flag_2meas_S = [];
flag_3meas_S = [];
flag_4meas_S = [];

for k = 1:length(sereja_mdek_full)
    cur_pos = sereja_mdek_full(k).anc_dist; 
    NoM = length(find(cur_pos));
    switch NoM
        case 1
             flag_1meas_S = [flag_1meas_S ser_sec_wk(k)];
        case 2
             flag_2meas_S = [flag_2meas_S ser_sec_wk(k)];
        case 3
             flag_3meas_S = [flag_3meas_S ser_sec_wk(k)];
        case 4
             flag_4meas_S = [flag_4meas_S ser_sec_wk(k)];
    end
end

ECEF_coords_N_sa(1,:) = round(ECEF_coords_N_sa(1,:),1);
ECEF_coords_S_sa(1,:) = round(ECEF_coords_S_sa(1,:),1);

common_time_nik = sort([nik_sec_wk ECEF_coords_N_sa(1,:)]);
common_time_ser = sort([ser_sec_wk ECEF_coords_S_sa(1,:)]);

common_time = sort([nik_sec_wk ECEF_coords_N_sa(1,:)]);

base_mean_ecef = [mean(ECEF_coords_B(2,:)) mean(ECEF_coords_B(3,:)) mean(ECEF_coords_B(4,:))];
base_mean_lla = ecef2lla(base_mean_ecef);
wgs84 = wgs84Ellipsoid;
for i = 1: length(ECEF_coords_N)
    [N_enu(1,i),N_enu(2,i),N_enu(3,i)] = ecef2enu(ECEF_coords_N(2,i),ECEF_coords_N(3,i),ECEF_coords_N(4,i),base_mean_lla(1),base_mean_lla(2),base_mean_lla(3), wgs84);
end
for i = 1: length(ECEF_coords_S)
    [S_enu(1,i),S_enu(2,i),S_enu(3,i)] = ecef2enu(ECEF_coords_S(2,i),ECEF_coords_S(3,i),ECEF_coords_S(4,i),base_mean_lla(1),base_mean_lla(2),base_mean_lla(3), wgs84);
end

for i = 1: length(ECEF_coords_N_sa)
    [N_coord_sa_enu(1,i),N_coord_sa_enu(2,i),N_coord_sa_enu(3,i)] = ecef2enu(ECEF_coords_N_sa(2,i),ECEF_coords_N_sa(3,i),ECEF_coords_N_sa(4,i),base_mean_lla(1),base_mean_lla(2),base_mean_lla(3), wgs84);
    [N_velo_sa_enu(1,i),N_velo_sa_enu(2,i),N_velo_sa_enu(3,i)] = ecef2enuv(ECEF_coords_N_sa(5,i),ECEF_coords_N_sa(6,i),ECEF_coords_N_sa(7,i),base_mean_lla(1),base_mean_lla(2));
end
for i = 1: length(ECEF_coords_S_sa)
    [S_coord_sa_enu(1,i),S_coord_sa_enu(2,i),S_coord_sa_enu(3,i)] = ecef2enu(ECEF_coords_S_sa(2,i),ECEF_coords_S_sa(3,i),ECEF_coords_S_sa(4,i),base_mean_lla(1),base_mean_lla(2),base_mean_lla(3), wgs84);
    [S_velo_sa_enu(1,i),S_velo_sa_enu(2,i),S_velo_sa_enu(3,i)] = ecef2enuv(ECEF_coords_S_sa(5,i),ECEF_coords_S_sa(6,i),ECEF_coords_S_sa(7,i),base_mean_lla(1),base_mean_lla(2));
end

R_V1V2 = sqrt((N_coord_sa_enu(1,:)-S_coord_sa_enu(1,:)).^2 + (N_coord_sa_enu(2,:)-S_coord_sa_enu(2,:)).^2);   
R_V1V2_rtk = sqrt((N_enu(1,:)-S_enu(1,:)).^2 + (N_enu(2,:)-S_enu(2,:)).^2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mdek_postproc

nikita_mdek_kf = [];
x_est = [nikita_mdek_full(1).pos(1); 0; nikita_mdek_full(1).pos(2); 0; nikita_mdek_full(1).pos(3); 0];
D_est = 0.1*eye(6);
nikita_mdek_pos = [];
for k = 1:length(nikita_mdek_full)
%     y_meas_k = nikita_mdek_full(k).anc_dist; 
%     NoM = length(find(y_meas_k));
%     if NoM > 2
%         SatPos = nikita_mdek_full(k).anc_pos';
%         SatPos = SatPos(:,find(y_meas_k));
%         y_meas_k = y_meas_k(find(y_meas_k));
%         [UserPos] = LSM_xyz(SatPos, y_meas_k, InitPos);
%     end
%     InitPos = UserPos;
%     nikita_mdek_LSM = [nikita_mdek_LSM [nik_sec_wk(k); UserPos]];

    y_meas_k = nikita_mdek_full(k).anc_dist; 
    SatPos = nikita_mdek_full(k).anc_pos';
    SatPos = SatPos(:,find(y_meas_k));
    y_meas_k = y_meas_k(find(y_meas_k));
    if k == 1
        dT = 0.1;
    else
        dT = nik_sec_wk(k) - nik_sec_wk(k-1);
    end
    [x_est,D_est] = mdek_range_kalman_filter(x_est,D_est,y_meas_k,dT, SatPos);
    nikita_mdek_kf = [nikita_mdek_kf  [nik_sec_wk(k); x_est]];
    if ~isempty(nikita_mdek_full(k).pos)
        nikita_mdek_pos = [nikita_mdek_pos [nik_sec_wk(k); nikita_mdek_full(k).pos]];
    end
end

sereja_mdek_kf = [];
x_est = [sereja_mdek_full(2).pos(1); 0; sereja_mdek_full(2).pos(2); 0; sereja_mdek_full(2).pos(3); 0];
D_est = 0.1*eye(6);
sereja_mdek_pos = [];
for k = 1:length(sereja_mdek_full)
%     y_meas_k = nikita_mdek_full(k).anc_dist; 
%     NoM = length(find(y_meas_k));
%     if NoM > 2
%         SatPos = nikita_mdek_full(k).anc_pos';
%         SatPos = SatPos(:,find(y_meas_k));
%         y_meas_k = y_meas_k(find(y_meas_k));
%         [UserPos] = LSM_xyz(SatPos, y_meas_k, InitPos);
%     end
%     InitPos = UserPos;
%     nikita_mdek_LSM = [nikita_mdek_LSM [nik_sec_wk(k); UserPos]];

    y_meas_k = sereja_mdek_full(k).anc_dist; 
    SatPos = sereja_mdek_full(k).anc_pos';
    SatPos = SatPos(:,find(y_meas_k));
    y_meas_k = y_meas_k(find(y_meas_k));
    if k == 1
        dT = 0.1;
    else
        dT = ser_sec_wk(k) - ser_sec_wk(k-1);
    end
    [x_est,D_est] = mdek_range_kalman_filter(x_est,D_est,y_meas_k,dT, SatPos);
    sereja_mdek_kf = [sereja_mdek_kf  [ser_sec_wk(k); x_est]];
    if ~isempty(sereja_mdek_full(k).pos)
        sereja_mdek_pos = [sereja_mdek_pos [ser_sec_wk(k); sereja_mdek_full(k).pos]];
    end
end

figure
plot(nikita_mdek_pos(1,:),nikita_mdek_pos(2,:))
hold on
plot(nikita_mdek_kf(1,:),nikita_mdek_kf(2,:))
grid on

figure
plot(nikita_mdek_pos(1,:),nikita_mdek_pos(3,:))
hold on
plot(nikita_mdek_kf(1,:),nikita_mdek_kf(4,:))
grid on

figure
plot(nikita_mdek_pos(1,:),nikita_mdek_pos(4,:))
hold on
plot(nikita_mdek_kf(1,:),nikita_mdek_kf(6,:))
grid on

figure
plot(sereja_mdek_pos(1,:),sereja_mdek_pos(2,:))
hold on
plot(sereja_mdek_kf(1,:),sereja_mdek_kf(2,:))
grid on

figure
plot(sereja_mdek_pos(1,:),sereja_mdek_pos(3,:))
hold on
plot(sereja_mdek_kf(1,:),sereja_mdek_kf(4,:))
grid on

figure
plot(sereja_mdek_pos(1,:),sereja_mdek_pos(4,:))
hold on
plot(sereja_mdek_kf(1,:),sereja_mdek_kf(6,:))
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% static_gnss_ind_N = 1:3692;
% std_n_x_N = 10*std(N_coord_sa_enu(1,static_gnss_ind_N));
% std_n_y_N = 10*std(N_coord_sa_enu(2,static_gnss_ind_N));
% std_n_Vx_N = 10*std(N_velo_sa_enu(1,static_gnss_ind_N));
% std_n_Vy_N = 10*std(N_velo_sa_enu(2,static_gnss_ind_N));

% std_n_x_N = 1;
% std_n_y_N = 1;
% std_n_Vx_N = 0.1;
% std_n_Vy_N = 0.1;
% 
% y_meas_for_kf_N = [N_coord_sa_enu(1,:);N_velo_sa_enu(1,:);N_coord_sa_enu(2,:);N_velo_sa_enu(2,:)];
% t_meas_for_kf_N = ECEF_coords_N_sa(1,:) - ECEF_coords_N_sa(1,1);
% x_est_init_N = [N_coord_sa_enu(1,1);N_velo_sa_enu(1,1); 0; N_coord_sa_enu(2,1);N_velo_sa_enu(2,1); 0];
% D_est_init_N = [1 0 0 0 0 0;
%                 0 1 0 0 0 0;
%                 0 0 1 0 0 0;
%                 0 0 0 1 0 0;
%                 0 0 0 0 1 0;
%                 0 0 0 0 0 1];
% [x_est_stor_N] = gnss_kalman_filter_coord_velo(x_est_init_N,D_est_init_N,y_meas_for_kf_N,t_meas_for_kf_N,std_n_x_N,std_n_y_N,std_n_Vx_N,std_n_Vy_N);
% 
% figure
% plot(t_meas_for_kf_N,x_est_stor_N(1,:))
% hold on
% plot(t_meas_for_kf_N,N_coord_sa_enu(1,:))
% grid on
% 
% figure
% plot(t_meas_for_kf_N,x_est_stor_N(2,:))
% hold on
% plot(t_meas_for_kf_N,N_velo_sa_enu(1,:))
% grid on
% 
% figure
% plot(t_meas_for_kf_N,x_est_stor_N(4,:))
% hold on
% plot(t_meas_for_kf_N,N_coord_sa_enu(2,:))
% grid on
% 
% figure
% plot(t_meas_for_kf_N,x_est_stor_N(5,:))
% hold on
% plot(t_meas_for_kf_N,N_velo_sa_enu(2,:))
% grid on

% A_gnss_coord_N = [max(nikita_mdek_full(1).anc_pos(1,:))/2 max(nikita_mdek_full(1).anc_pos(2,:))/2 0];
% A_gnss_coord_S = [max(sereja_mdek_full(1).anc_pos(1,:))/2 max(sereja_mdek_full(1).anc_pos(2,:))/2 0];

A_gnss_coord_N = [0.57 0.34 0];
A_gnss_coord_S = [0.54 0.67 0];

x_est_kf_top = [N_coord_sa_enu(1,1);N_velo_sa_enu(1,1);0;N_coord_sa_enu(2,1);N_velo_sa_enu(2,1);0;S_coord_sa_enu(1,1);S_velo_sa_enu(1,1);0;S_coord_sa_enu(2,1);S_velo_sa_enu(2,1);0];
D_est_kf_top = eye(12);
x_est_kf_top_stor = [];
for k = 1:length(common_time)
    if k ==493
        fgf = 1;
    end
       if k == 1
           t_prev = common_time(1)-0.1;
           t_cur = common_time(1);
       else
           t_cur = common_time(k);            
       end
    if ~isempty(find(ECEF_coords_N_sa(1,:) == common_time(k)))
        flag_meas = 1;
        gnss_ind = find(ECEF_coords_N_sa(1,:) == common_time(k),1);      
        y_meas = [N_coord_sa_enu(1,gnss_ind);N_velo_sa_enu(1,gnss_ind);N_coord_sa_enu(2,gnss_ind);N_velo_sa_enu(2,gnss_ind);S_coord_sa_enu(1,gnss_ind);S_velo_sa_enu(1,gnss_ind);S_coord_sa_enu(2,gnss_ind);S_velo_sa_enu(2,gnss_ind)];
        [x_est_kf_top,D_est_kf_top] = gnss_mdek_coop_nav_kf(x_est_kf_top,D_est_kf_top,y_meas,t_cur,t_prev,flag_meas);
    else
        flag_meas = 2;
        mdek_ind = find(nik_sec_wk == common_time(k),1);
        y_meas = sqrt((A_gnss_coord_N(1) - nikita_mdek_kf(2,mdek_ind))^2 + (A_gnss_coord_N(2) - nikita_mdek_kf(3,mdek_ind))^2 + (A_gnss_coord_N(3) - nikita_mdek_kf(4,mdek_ind))^2);
        [x_est_kf_top,D_est_kf_top] = gnss_mdek_coop_nav_kf(x_est_kf_top,D_est_kf_top,y_meas,t_cur,t_prev,flag_meas);        
    end
    t_prev = common_time(k);
    x_est_kf_top_stor = [x_est_kf_top_stor x_est_kf_top];
end

figure
plot(common_time,x_est_kf_top_stor(1,:))
hold on
plot(ECEF_coords_N_sa(1,:),N_coord_sa_enu(1,:))
plot(ECEF_coords_N(1,:),N_enu(1,:))
grid on

figure
plot(common_time,x_est_kf_top_stor(4,:))
hold on
plot(ECEF_coords_N_sa(1,:),N_coord_sa_enu(2,:))
plot(ECEF_coords_N(1,:),N_enu(2,:))
grid on

figure
plot(common_time,x_est_kf_top_stor(7,:))
hold on
plot(ECEF_coords_S_sa(1,:),S_coord_sa_enu(1,:))
plot(ECEF_coords_S(1,:),S_enu(1,:))
grid on

figure
plot(common_time,x_est_kf_top_stor(10,:))
hold on
plot(ECEF_coords_S_sa(1,:),S_coord_sa_enu(2,:))
plot(ECEF_coords_S(1,:),S_enu(2,:))
grid on

x_est_R_V1V2 = sqrt((x_est_kf_top_stor(1,:)-x_est_kf_top_stor(7,:)).^2 + (x_est_kf_top_stor(4,:)-x_est_kf_top_stor(10,:)).^2);

figure
plot(common_time,x_est_R_V1V2)
hold on
plot(ECEF_coords_N_sa(1,:),R_V1V2)
plot(ECEF_coords_N(1,:),R_V1V2_rtk)
grid on
%%
figure
plot(N_enu(1,:),N_enu(2,:),'-.')
hold on
plot(S_enu(1,:),S_enu(2,:),'-.')
grid on

figure
plot(N_coord_sa_enu(1,:),N_coord_sa_enu(2,:),'-.')
hold on
plot(S_coord_sa_enu(1,:),S_coord_sa_enu(2,:),'-.')
plot(x_est_kf_top_stor(1,:),x_est_kf_top_stor(4,:))
plot(x_est_kf_top_stor(7,:),x_est_kf_top_stor(10,:))
grid on

figure
plot(N_enu(1,:),N_enu(2,:),'-.')
hold on
plot(N_coord_sa_enu(1,:),N_coord_sa_enu(2,:),'-.')
grid on

figure
plot(S_enu(1,:),S_enu(2,:),'-.')
hold on
plot(S_coord_sa_enu(1,:),S_coord_sa_enu(2,:),'-.')
grid on

figure
plot(ECEF_coords_N_sa(1,:),N_coord_sa_enu(1,:))
hold on
plot(ECEF_coords_N_sa(1,:),N_coord_sa_enu(2,:))
plot(ECEF_coords_N_sa(1,:),N_coord_sa_enu(3,:))
grid on

figure
plot(ECEF_coords_S_sa(1,:),S_coord_sa_enu(1,:))
hold on
plot(ECEF_coords_S_sa(1,:),S_coord_sa_enu(2,:))
plot(ECEF_coords_S_sa(1,:),S_coord_sa_enu(3,:))
grid on

figure
plot(ECEF_coords_N_sa(1,:),N_velo_sa_enu(1,:))
hold on
plot(ECEF_coords_N_sa(1,:),N_velo_sa_enu(2,:))
plot(ECEF_coords_N_sa(1,:),N_velo_sa_enu(3,:))
grid on

figure
plot(ECEF_coords_S_sa(1,:),S_velo_sa_enu(1,:))
hold on
plot(ECEF_coords_S_sa(1,:),S_velo_sa_enu(2,:))
plot(ECEF_coords_S_sa(1,:),S_velo_sa_enu(3,:))
grid on

