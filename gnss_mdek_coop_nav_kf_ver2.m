function [x_est,D_est] = gnss_mdek_coop_nav_kf_ver2(x_est,D_est,y_meas,t_cur,t_prev,flag_meas)
   D_ksi = [2^2  0   0   0;
             0  2^2  0   0;
             0   0  1^2  0;
             0   0   0  1^2];

   I = eye(6); 
   
    dT = t_cur-t_prev;
        
        F = [1 dT   0     0   0    0;
             0  1   0     0   0    0;
             0  0   1    dT   0    0;
             0  0   0     1   0    0;
             0  0   0     0   1    0;
             0  0   0     0   0    1];
        G = [0  0  0  0;
             dT 0  0  0;
             0  0  0  0;
             0 dT  0  0;
             0  0  dT 0;
             0  0  0 dT];
         
        x_ext = F*x_est;
        D_ext = F*D_est*F' + G*D_ksi*G';
        if flag_meas == 1
            sigma_n = [1 0.5 1 0.5];
            Dn = eye(4);
            for q = 1:size(Dn,1)
                Dn(q,q) = sigma_n(q);
            end
            H = [1 0 0 0 0 0;
                 0 1 0 0 0 0;
                 0 0 1 0 0 0;
                 0 0 0 1 0 0];
            f_x_ext = zeros(4,1);
            f_x_ext(1:4,1) = H*x_ext; 
        else
            sigma_n = 0.3;
            Dn = [sigma_n^2 0;
                   0 sigma_n^2];
            H = zeros(2,6);
            H(1,1) = x_ext(5)/sqrt(x_ext(5)^2 + x_ext(6)^2);
            H(1,3) = x_ext(6)/sqrt(x_ext(5)^2 + x_ext(6)^2); 
            H(1,5) = x_ext(5)/sqrt(x_ext(5)^2 + x_ext(6)^2); 
            H(1,6) = x_ext(6)/sqrt(x_ext(5)^2 + x_ext(6)^2);  
            H(2,1) = -x_ext(6)/(x_ext(5)^2 + x_ext(6)^2);
            H(2,3) = x_ext(5)/(x_ext(5)^2 + x_ext(6)^2); 
            H(2,5) = -x_ext(6)/(x_ext(5)^2 + x_ext(6)^2); 
            H(2,6) = x_ext(5)/(x_ext(5)^2 + x_ext(6)^2); 
            f_x_ext = zeros(2,1);
            f_x_ext(1,1) = sqrt(x_ext(5)^2 + x_ext(6)^2);
            f_x_ext(2,1) = atan2(x_ext(6),x_ext(5));
        end
        K = D_ext*H'*(H*D_ext*H'+Dn)^-1;
        D_est = (I - K*H)*D_ext;
        x_est = x_ext + K*(y_meas - f_x_ext);   
end