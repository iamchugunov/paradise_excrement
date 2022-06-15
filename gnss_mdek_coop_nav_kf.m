function [x_est,D_est] = gnss_mdek_coop_nav_kf(x_est,D_est,y_meas,t_cur,t_prev,flag_meas)
   D_ksi = [0.1^2  0   0   0;
             0  0.1^2  0   0;
             0   0  0.1^2  0;
             0   0   0  0.1^2];

   I = eye(12); 
   
    dT = t_cur-t_prev;
        
        F = [1 dT dT^2/2 0   0    0   0  0   0    0   0   0;
             0  1   dT   0   0    0   0  0   0    0   0   0;
             0  0   1    0   0    0   0  0   0    0   0   0;
             0  0   0    1  dT dT^2/2 0  0   0    0   0   0;
             0  0   0    0   1   dT   0  0   0    0   0   0;
             0  0   0    0   0    1   0  0   0    0   0   0;
             0  0   0    0   0    0   1 dT dT^2/2 0   0   0;
             0  0   0    0   0    0   0  1   dT   0   0   0;
             0  0   0    0   0    0   0  0   1    0   0   0;
             0  0   0    0   0    0   0  0   0    1  dT dT^2/2;
             0  0   0    0   0    0   0  0   0    0   1   dT;
             0  0   0    0   0    0   0  0   0    0   0   1];
        G = [0  0  0  0;
             0  0  0  0;
             dT 0  0  0;
             0  0  0  0;
             0  0  0  0;
             0 dT  0  0;
             0  0  0  0;
             0  0  0  0;
             0  0  dT 0;
             0  0  0  0;
             0  0  0  0;
             0  0  0 dT];
         
        x_ext = F*x_est;
        D_ext = F*D_est*F' + G*D_ksi*G';
        if flag_meas == 1
            sigma_n = [1 0.1 1 0.1 1 0.1 1 0.1];
            Dn = eye(8);
            for q = 1:size(Dn,1)
                Dn(q,q) = sigma_n(q);
            end
            H = [1 0 0 0 0 0 0 0 0 0 0 0;
                 0 1 0 0 0 0 0 0 0 0 0 0;
                 0 0 0 1 0 0 0 0 0 0 0 0;
                 0 0 0 0 1 0 0 0 0 0 0 0;
                 0 0 0 0 0 0 1 0 0 0 0 0;
                 0 0 0 0 0 0 0 1 0 0 0 0;
                 0 0 0 0 0 0 0 0 0 1 0 0;
                 0 0 0 0 0 0 0 0 0 0 1 0];
            f_x_ext = zeros(8,1);
            f_x_ext(1:8,1) = H*x_ext; 
        else
            sigma_n = 0.1;
            Dn = sigma_n^2;
            H = zeros(1,12);
            H(1,1) = (x_ext(1)-x_ext(7))/sqrt((x_ext(1)-x_ext(7))^2 + (x_ext(4)-x_ext(10))^2);
            H(1,4) = (x_ext(4)-x_ext(10))/sqrt((x_ext(1)-x_ext(7))^2 + (x_ext(4)-x_ext(10))^2); 
            H(1,7) = -(x_ext(1)-x_ext(7))/sqrt((x_ext(1)-x_ext(7))^2 + (x_ext(4)-x_ext(10))^2); 
            H(1,10) = -(x_ext(4)-x_ext(10))/sqrt((x_ext(1)-x_ext(7))^2 + (x_ext(4)-x_ext(10))^2);  
            f_x_ext = sqrt((x_ext(1)-x_ext(7))^2 + (x_ext(4)-x_ext(10))^2);
        end
        K = D_ext*H'*(H*D_ext*H'+Dn)^-1;
        D_est = (I - K*H)*D_ext;
        x_est = x_ext + K*(y_meas - f_x_ext);   
end