function [x_est_stor] = gnss_kalman_filter_coord_velo(x_est,D_est,y_meas,t,std_n_x_N,std_n_y_N,std_n_Vx_N,std_n_Vy_N)
   D_ksi = [1^2 0;
             0   1^2];
   Dn = [ std_n_x_N^2  0      0      0;
           0   std_n_Vx_N^2   0      0;
           0    0     std_n_y_N^2    0;
           0    0      0    std_n_Vy_N^2];

   I = eye(6); 
   
   x_est_stor = [];
   
    for k = 1:length(y_meas)
        if k == 1
            dT = 0.2;
        else
            dT = t(k)-t(k-1);
        end
        
        F = [1 dT dT^2/2 0   0    0;
             0  1   dT   0   0    0;
             0  0   1    0   0    0;
             0  0   0    1  dT dT^2/2;
             0  0   0    0   1   dT;
             0  0   0    0   0    1];
        G = [0  0;
             0  0;
             dT 0;
             0  0;
             0  0;
             0 dT];
         
        x_ext = F*x_est;
        D_ext = F*D_est*F' + G*D_ksi*G';
        
        H = [1 0 0 0 0 0;
             0 1 0 0 0 0;
             0 0 0 1 0 0;
             0 0 0 0 1 0];
       
        K = D_ext*H'*(H*D_ext*H'+Dn)^-1;
        D_est = (I - K*H)*D_ext;
        x_est = x_ext + K*(y_meas(:,k) - H*x_ext);
        x_est_stor = [x_est_stor x_est];
    end
    
end