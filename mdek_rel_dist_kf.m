function [x_est,D_est] = mdek_rel_dist_kf(x_est,D_est,y_meas,delta_R,delta_alp,dT)
   D_ksi = [2^2  0;
             0  0.5^2;];

   I = eye(4); 
        
        F = [1 dT   0   0;
             0  1   0   0;
             0  0   1  dT;
             0  0   0   1];
        G = [0  0;
             dT 0;
             0  0;
             0 dT];
         
        x_ext = F*x_est;
        D_ext = F*D_est*F' + G*D_ksi*G';
               
        Dn = 0.1*eye(length(y_meas));
        H = zeros(length(y_meas),4);
        f_x_ext = zeros(length(y_meas),1);
        for q = 1:length(y_meas)
            H(q,:) = [(x_ext(1) - delta_R(q)*cos(x_ext(3) + delta_alp(q)))/sqrt(x_ext(1)^2 + delta_R(q)^2 - 2*x_ext(1)*delta_R(q)*cos(x_ext(3) + delta_alp(q))) 0 x_ext(1)*delta_R(q)*sin(x_ext(3) + delta_alp(q))/sqrt(x_ext(1)^2 + delta_R(q)^2 - 2*x_ext(1)*delta_R(q)*cos(x_ext(3) + delta_alp(q))) 0];
            f_x_ext(q,1) = sqrt(x_ext(1)^2 + delta_R(q)^2 - 2*x_ext(1)*delta_R(q)*cos(x_ext(3) + delta_alp(q)));             
        end


%             sigma_n = 0.1;
%             Dn = sigma_n^2;
%             H = zeros(1,12);
%             H(1,1) = (x_ext(1)-x_ext(7))/sqrt((x_ext(1)-x_ext(7))^2 + (x_ext(4)-x_ext(10))^2);
%             H(1,4) = (x_ext(4)-x_ext(10))/sqrt((x_ext(1)-x_ext(7))^2 + (x_ext(4)-x_ext(10))^2); 
%             H(1,7) = -(x_ext(1)-x_ext(7))/sqrt((x_ext(1)-x_ext(7))^2 + (x_ext(4)-x_ext(10))^2); 
%             H(1,10) = -(x_ext(4)-x_ext(10))/sqrt((x_ext(1)-x_ext(7))^2 + (x_ext(4)-x_ext(10))^2);  
%             f_x_ext = sqrt((x_ext(1)-x_ext(7))^2 + (x_ext(4)-x_ext(10))^2);

        K = D_ext*H'*(H*D_ext*H'+Dn)^-1;
        D_est = (I - K*H)*D_ext;
        x_est = x_ext + K*(y_meas - f_x_ext);   
end