%%
[mdek_N, mdek_S, RTK_N, RTK_S, SA_N, SA_S, BASE] = open_logs();
for i = 1:length(mdek_S)
    mdek_S(i).sec_wk = mdek_S(i).sec_wk - 4;
end
%% построить все
show_td(mdek_N, mdek_S, SA_N, SA_S, RTK_N, RTK_S, 1)

%%
t0 = 38146.22;
% T = [39893 40108]; % первый тестовый, много мдеков
% T = [200 350]; % друг за другом на большом расстоянии
% T = [450 750]; % долгая статика
% T = [450 560]; % урезанная предыдущая статика статика
% T = [828 1000]; % как будто едут рядом
% T = [1000 1200]; % рядом прямо и в поворот | Никита + Серега -
% T = [1300 1700]; %долгая статика

% T = [1700 2135]; % туда сюда друг за другом | Никита + Серега + лучше разделить на два
T = [1755 1900]; % туда сюда друг за другом | Никита + Серега + первый участок
T = [1960 2100]; % туда сюда друг за другом | Никита + Серега + второй участок

% T = [2390 2630]; % туда сюда друг за другом как будто быстрее чем предыдущий | Никита + Серега + лучше разделить на два
% T = [2525 2570]; % параллельная езда по прямой
% T = [2460 2480]; % по прямой друг за другом короткий | Никита + Серега -
% T = [2790 2810]; % проезд мимо стоячего на повороте | Никита + Серега +-
T = [2930 2960]; % проезд мимо стоячего на повороте в обратную сторону | Никита + Серега +
% T = [3170 3190]; % навстречу в повороте | Никита +- Серега +-
% T = [3260 3280]; % навстречу в повороте х2 | Никита - Серега +-
% T = [3355 3375]; % навстречу в повороте х3 | Никита - Серега +-
% T = [2630 3420]; % предыдущие 5 вместе
% T = [3850 3950];%долгая статика
% T = [3942 4450]; % навстречу в повороте c остановкой
% T = [4450 4962]; % 4 обгона в повороте

% T = [1 4962];

[mdek_n, mdek_s, sa_n, sa_s, rtk_n, rtk_s] = cut_by_time(mdek_N, mdek_S, SA_N, SA_S, RTK_N, RTK_S, T + t0); % режем по времени
show_td(mdek_n, mdek_s, sa_n, sa_s, rtk_n, rtk_s, 1)
%%
mdek_s_post = [];
mdek_n_post = [];
mdek_s_post.res_mnk = calculate_mdek_dataset_mnk_decart(mdek_s);
mdek_n_post.res_mnk = calculate_mdek_dataset_mnk_decart(mdek_n);
mdek_s_post.res_kf = do_kalman_filter_mdek_coords(mdek_s);
mdek_n_post.res_kf = do_kalman_filter_mdek_coords(mdek_n);
post_analys(mdek_s_post, rtk_n, rtk_s)