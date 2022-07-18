function [x_est,D_est] = mdek_range_kalman_filter(x_est,D_est,y_meas,dT,SatPos, sigma_ksi, sigma_n)
   D_ksi = sigma_ksi^2 * eye(2);
   
   Dn = sigma_n*eye(length(y_meas));

   I = eye(4);   
        
%         F = [1 dT   0    0   0    0;
%              0  1   0    0   0    0;
%              0  0   1    dT  0    0;
%              0  0   0    1   0    0;
%              0  0   0    0   1   dT;
%              0  0   0    0   0    1];
%         G = [0  0  0;
%              dT 0  0;
%              0  0  0;
%              0  dT 0;
%              0  0  0;
%              0  0 dT];
        F = [1 dT   0    0;
             0  1   0    0;
             0  0   1    dT;
             0  0   0    1;];
        G = [0  0;
             dT 0;
             0  0;
             0  dT];
         
        x_ext = F*x_est;
        D_ext = F*D_est*F' + G*D_ksi*G';
        H = zeros(length(y_meas),4);
        f_x = zeros(length(y_meas),1);
        for q = 1:length(y_meas)
            H(q,1) = (x_ext(1) - SatPos(1,q))/sqrt((x_ext(1) - SatPos(1,q))^2 + (x_ext(3) - SatPos(2,q))^2 + (0 - SatPos(3,q))^2);
            H(q,3) = (x_ext(3) - SatPos(2,q))/sqrt((x_ext(1) - SatPos(1,q))^2 + (x_ext(3) - SatPos(2,q))^2 + (0 - SatPos(3,q))^2);
            f_x(q,1) = sqrt((x_ext(1) - SatPos(1,q))^2 + (x_ext(3) - SatPos(2,q))^2 + (0 - SatPos(3,q))^2);
        end
        K = D_ext*H'*(H*D_ext*H'+Dn)^-1;
        D_est = (I - K*H)*D_ext;
        x_est = x_ext + K*(y_meas - f_x);
    
end