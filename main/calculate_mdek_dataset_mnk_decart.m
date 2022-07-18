function [res] = calculate_mdek_dataset_mnk_decart(mdek)
    t = [];
    ranges = [];
    k = 0;
    nev = [];
    pos = [];
    R = [];
    alp = [];
    
    for i = 1:length(mdek)
        if length(find(mdek(i).anc_dist)) == 4
            posts = mdek(i).anc_pos';
            posts = posts - mean(posts')';
            break
        end
    end
    
    for i = 1:length(mdek)
        if length(find(mdek(i).anc_dist)) > 2
            
            t_cur = mdek(i).sec_wk;
            ranges_cur = mdek(i).anc_dist;
            nms = find(ranges_cur);
            ranges_cur1 = ranges_cur(nms);
            posts_cur = posts(:,nms);
            
            [X, dop, nev(i), flag] = NavSolver_D(ranges_cur1, posts_cur, [0;0], 0);
            if flag
                k = k + 1;
                pos(:,k) = -X;
                t(k) = t_cur;
%                 ranges(:,k) = ranges_cur;
                alp(k) = atan2(pos(2,k),pos(1,k));
                R(k) = norm(pos(1:2,k));
            else
            end
        end
    end
    
    res = [];
    res.t = t;
    res.pos = pos;
    res.R = R;
    res.alp = alp;
    res.k = k;

end




