function [] = show_td(mdek_N, mdek_S, SA_N, SA_S, RTK_N, RTK_S, flag)
    mdek_n = [];
    for i = 1:length(mdek_N)
        mdek_n(:,i) = [mdek_N(i).sec_wk; mdek_N(i).anc_dist];
    end
    
    mdek_s = [];
    for i = 1:length(mdek_S)
        mdek_s(:,i) = [mdek_S(i).sec_wk; mdek_S(i).anc_dist];
    end
    
%     mdek_s(1,:) = mdek_s(1,:) - 4; % пальцем

    if length(mdek_N) > 10
        N1 = length(mdek_N);
        N_1 = (mdek_N(end).sec_wk - mdek_N(1).sec_wk)/0.1;
        [N1; round(N1/N_1*100,1)]
        kol1 = [length(find([mdek_N.count]==1)) length(find([mdek_N.count]==2)) length(find([mdek_N.count]==3)) length(find([mdek_N.count]==4))];
        per1 = round(kol1/N1*100,1);
        C1 = [kol1; per1]
    else
        N1 = length(mdek_N)
    end
    
    if length(mdek_s) > 10
        N2 = length(mdek_S);
        N_2 = (mdek_S(end).sec_wk - mdek_S(1).sec_wk)/0.1;
        [N2; round(N2/N_1*100,1)]
        kol2 = [length(find([mdek_S.count]==1)) length(find([mdek_S.count]==2)) length(find([mdek_S.count]==3)) length(find([mdek_S.count]==4))];
        per2 = round(kol2/N2*100,1);
        C2 = [kol2; per2]
    else
        N2 = length(mdek_S)
    end
    
    R_SA = [];
    X_SA_N = [];
    X_SA_S = [];
    t1 = round(SA_N(1,:),1);
    t2 = round(SA_S(1,:),1);
    t = intersect(t1, t2);
    for i = 1:length(t)
        k1 = find(t1 == t(i));
        k2 = find(t2 == t(i));
        X_SA_N(:,i) = SA_N([2,4],k1);
        X_SA_S(:,i) = SA_S([2,4],k2);
        R_SA(i) = norm(SA_N([2 4 6],k1) - SA_S([2 4 6],k2));
    end
    
    R_RTK = [];
    t1 = round(RTK_N(1,:),1);
    t2 = round(RTK_S(1,:),1);
    t_RTK = intersect(t1, t2);
    for i = 1:length(t_RTK)
        k1 = find(t1 == t_RTK(i));
        k2 = find(t2 == t_RTK(i));
        R_RTK(i) = norm(RTK_N([2 3 4],k1) - RTK_S([2 3 4],k2));
    end
    
%     t0 = min([mdek_n(1,1) mdek_s(1,1) t(1)]);
    t0 = 0;
    t0 = 3.814622000000000e+04;
    
%     if nargin == 7
%         nums = intersect(find(mdek_n(1,:) > time_interval(1)),find(mdek_n(1,:) < time_interval(2)));
%         mdek_n = mdek_n(:,nums);
%         
%         nums = intersect(find(mdek_s(1,:) > time_interval(1)),find(mdek_s(1,:) < time_interval(2)));
%         mdek_s = mdek_s(:,nums);
%         
%         nums = intersect(find(t > time_interval(1)),find(t < time_interval(2)));
%         t = t(:,nums);
%         R_SA = R_SA(:,nums);
%         X_SA_N = X_SA_N(:,nums);
%         X_SA_S = X_SA_S(:,nums);
%         
%         nums = intersect(find(t_RTK > time_interval(1)),find(t_RTK < time_interval(2)));
%         t_RTK = t_RTK(:,nums);
%         R_RTK = R_RTK(:,nums);
%         
%         nums = intersect(find(SA_N(1,:) > time_interval(1)),find(SA_N(1,:) < time_interval(2)));
%         SA_N = SA_N(:,nums);        
%              
%         nums = intersect(find(SA_S(1,:) > time_interval(1)),find(SA_S(1,:) < time_interval(2)));
%         SA_S = SA_S(:,nums);
%    
%     end
    
    if flag 
        figure
        plot(mdek_n(1,:)-t0, mdek_n(2:5,:),'x')
        hold on
        plot(mdek_s(1,:)-t0, mdek_s(2:5,:),'.','linewidth',3)
        p1 = plot(t-t0, R_SA,'.-k');
        p2 = plot(t_RTK-t0, R_RTK,'.-r');
        grid on
        xlabel('t, сек')
        ylabel('R, м')
%         legend([],[],[],[],[],[],[],[],"Stand alone",'RTK')
        legend([p1 p2],{'Stand alone','RTK'})
        return
    end
    
    figure
    subplot(211)
    plot(SA_N(2,:),SA_N(4,:),'r--')
    hold on
    h1 = plot(SA_N(2,1),SA_N(4,1),'rx','linewidth',4);
    plot(SA_S(2,:),SA_S(4,:),'b--')
    h2 = plot(SA_S(2,1),SA_S(4,1),'bx','linewidth',4);
%     daspect([1 1 1])
    grid on
    
    subplot(212)
    ylim([0 70])
    hold on
    grid on
    plot(t-t0, R_SA,'--k')
    h3 = plot(t(1)-t0, R_SA(1),'or');
    for i = 1:length(t)
         set(h1,'XData',X_SA_N(1,i),'YData',X_SA_N(2,i))
         set(h2,'XData',X_SA_S(1,i),'YData',X_SA_S(2,i))
         set(h3,'XData',t(i) - t0,'YData',R_SA(i))
         pause(0.02)
    end
    grid on
end

