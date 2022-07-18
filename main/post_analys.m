function [] = post_analys(mdek_post, rtk_n, rtk_s)
    
    R_RTK = [];
    t1 = round(rtk_n(1,:),1);
    t2 = round(rtk_s(1,:),1);
    t_RTK = intersect(t1, t2);
    for i = 1:length(t_RTK)
        k1 = find(t1 == t_RTK(i));
        k2 = find(t2 == t_RTK(i));
        R_RTK(i) = norm(rtk_n([2 3 4],k1) - rtk_s([2 3 4],k2));
    end
    
    t = mdek_post.res_kf.t;
    RTK_N(1,:) = interp1(t1, rtk_n(2,:), t);
    RTK_N(2,:) = interp1(t1, rtk_n(3,:), t);
    RTK_S(1,:) = interp1(t2, rtk_s(2,:), t);
    RTK_S(2,:) = interp1(t2, rtk_s(3,:), t);
    
    R_RTK_T = interp1(t_RTK, R_RTK, t);
    
    figure
    plot(t-t(1), R_RTK_T,'.-k','linewidth',2)
    hold on
    grid on
    plot(t-t(1), mdek_post.res_kf.R,'r','linewidth',2)
    xlabel('t, сек')
    ylabel('R, м')
    legend('ГНСС-RTK','СШП ЛНС')
    
    delta = mdek_post.res_kf.R - R_RTK_T;
    
    [mean(delta(10:end)) std(delta(10:end))]
    
    figure
    plot(t-t(1), delta, 'k')
    grid on
    xlabel('t, сек')
    ylabel('ΔR, м')
    
    
    
    ref = [56.0812781544862          37.3734438377441          217.698490284383];
    geo_n = [];
    geo_s = [];
    for i = 1:length(t)
        [geo_n(1,i), geo_n(2,i), geo_n(3,i)] = enu2geodetic(RTK_N(1,i), RTK_N(2,i),0,ref(1),ref(2),ref(3),wgs84Ellipsoid);
        [geo_s(1,i), geo_s(2,i), geo_s(3,i)] = enu2geodetic(RTK_S(1,i), RTK_S(2,i),0,ref(1),ref(2),ref(3),wgs84Ellipsoid);
    end
    figure('Position',[681 559 1162 420])
    subplot(121)
    geoplot(geo_n(1,:),geo_n(2,:),'--r')
    hold on
    geoplot(geo_s(1,:),geo_s(2,:),'--b')
    geobasemap 'streets'
    car_n = geoplot(geo_n(1,10),geo_n(2,10),'.r','MarkerSize',12);
    car_s = geoplot(geo_s(1,10),geo_s(2,10),'.b','MarkerSize',12);
    
    subplot(122)
    rectangle('Position',[-1 -1.5 2 3],'FaceColor',[0 0 0])
    axis([-30 30 -30 30])
    axis([-30 30 -30 30]*2)
    grid on
    daspect([1 1 1])
    hold on
    car = plot(mdek_post.res_kf.sv(1,1),mdek_post.res_kf.sv(3,1),'*r','MarkerSize',12);
    xlabel('x, м')
    ylabel('у, м')
    pause(1)
    for i = 10:length(t)
        set(car_n, "LatitudeData", geo_n(1,i), "LongitudeData", geo_n(2,i))
        set(car_s, "LatitudeData", geo_s(1,i), "LongitudeData", geo_s(2,i))
        set(car,"XData",mdek_post.res_kf.sv(1,i),"YData",mdek_post.res_kf.sv(3,i))
        pause(0.01)
    end
    
    
    
    
    
end

