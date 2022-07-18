function [filter_coords] = do_kalman_filter_mdek_r(mdek)
    
    anc0 = mean(mdek(1).anc_pos)';
    satpos = mdek(1).anc_pos';
    satpos = satpos -anc0;
    delta_R = [];
    delta_alp = [];
    for i = 1:4
        delta_R(i) = norm(satpos(:,i));
        delta_alp(i) = atan2(satpos(2,i),satpos(1,i));
    end
        
    x_max_anc_pos_nikita = mdek(1).anc_pos(2,1);
    y_max_anc_pos_nikita = mdek(1).anc_pos(2,2);
    delta_alp = zeros(4,1);
    delta_alp(1) = atan(0.5*x_max_anc_pos_nikita/(0.5*y_max_anc_pos_nikita)); atan(0.5*x_max_anc_pos_nikita/(0.5*y_max_anc_pos_nikita));
    delta_alp(2) = pi - delta_alp(1);
    delta_alp(3) = pi + delta_alp(1);
    delta_alp(4) = 2*pi + delta_alp(1);
            
    x_est = [max(mdek(1).anc_dist); 0; 0; 0];
    D_est = eye(4);


    SV = [mdek(1).sec_wk; x_est];

    for i = 2:length(mdek)
        dt = mdek(i).sec_wk - mdek(i - 1).sec_wk;
        N = find(mdek(i).anc_dist);
        y_meas = mdek(i).anc_dist(N);
        [x_est,D_est] = mdek_rel_dist_kf(x_est,D_est,y_meas,delta_R(N),delta_alp(N),dt);        
        SV(:,i) = [mdek(i).sec_wk; x_est];
    end
    filter_coords = SV;
        
end



