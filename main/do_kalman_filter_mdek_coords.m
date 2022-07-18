function [res] = do_kalman_filter_mdek_coords(mdek)
    res = [];
    
    for i = 1:length(mdek)
        if length(find(mdek(i).anc_dist)) == 4
            posts = mdek(i).anc_pos';
            posts = posts - mean(posts')';
            break
        end
    end
    
    i_start = 0;
    for i = 1:length(mdek)
        if length(find(mdek(i).anc_dist)) > 2
            
            t_cur = mdek(i).sec_wk;
            ranges_cur = mdek(i).anc_dist;
            nms = find(ranges_cur);
            ranges_cur1 = ranges_cur(nms);
            posts_cur = posts(:,nms);
            
            [X, dop, nev(i), flag] = NavSolver_D(ranges_cur1, posts_cur, [0;0], 0);
            if flag
                i_start = i+1;
                x_est = [X(1); 0; X(2); 0;];
                D_est = 0.1*eye(4);
                break
            end
        end
    end
    
    sigma_ksi = 1;
    sigma_n = 0.05;
    
    if i_start == 0
        return
    end
    
    k = 1;
    res.t(k) = mdek(i_start - 1).sec_wk;
    res.sv(:,k) = -x_est;
    for i = i_start:length(mdek)
        t_cur = mdek(i).sec_wk;
        ranges_cur = mdek(i).anc_dist;
        nms = find(ranges_cur);
        ranges_cur1 = ranges_cur(nms);
        posts_cur = posts(:,nms);
        dt = mdek(i).sec_wk - mdek(i - 1).sec_wk;
        [x_est,D_est] = mdek_range_kalman_filter(x_est,D_est,ranges_cur1,dt, posts_cur, sigma_ksi, sigma_n);
        k = k + 1;
        res.t(k) = t_cur;
        res.sv(:,k) = -x_est;
    end
    
    for i = 1:k
        res.R(i) = norm(res.sv([1 3],i));
        res.alp(i) = atan2(res.sv(3,i),res.sv(1,i));
        res.V(i) = norm(res.sv([2 4],i));
    end
    
end

