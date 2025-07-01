% ========= 提取原始 CHAZ 的 PTC 轨迹 =========
clear; clc

fpath = '/Users/dianyili/Desktop/SBU/CHAZ/atl_era5/';

% 参数设定
num_ens = 40;  % 集合成员数量
max_x = 700;   % 最大时间步（气旋数量或时间序列长度）
max_y = 125;   % 最大气旋长度（如一天一个点）

% 预分配存储变量
lonex1   = nan(max_x, max_y, num_ens);
latex1   = nan(max_x, max_y, num_ens);
years1   = nan(max_x, max_y, num_ens);
months1  = nan(max_x, max_y, num_ens);
days1    = nan(max_x, max_y, num_ens);
yr1      = nan(max_x, num_ens);  % 每行一个气旋的年份标签

% ===== 循环读取每个集合成员数据 =====
for i = 1:num_ens

    disp(['Processing Ensemble ', num2str(i)])

    % 构建文件名
    if i-1 < 10
        fname = ['atl_2019_2ens00', num2str(i-1), '_ETv1.nc'];
    else
        fname = ['atl_2019_2ens0', num2str(i-1), '_ETv1.nc'];
    end
    ff = [fpath, fname];

    % 读取数据
    lon = ncread(ff, 'longitude');
    lat = ncread(ff, 'latitude');
    yr = ncread(ff, 'year');
    ld = ncread(ff, 'ldmask');
    tm = ncread(ff, 'time');  % days since 1950-01-01
    tm(tm < 0) = nan;  % 移除无效值

    % 时间转换
    startDate = datetime(1950, 1, 1);
    dates = startDate + days(tm);
    [year, month, day] = ymd(dates);

    % 年份筛选
    valid_idx = yr >= 1981 & yr <= 2019;
    lat = lat(valid_idx, :);
    lon = lon(valid_idx, :);
    ld = ld(valid_idx, :);
    year = year(valid_idx, :);
    month = month(valid_idx, :);
    day = day(valid_idx, :);
    yr = yr(valid_idx);

    % 获取尺寸
    [nx, ny] = size(lat);

    % 保存当前成员的数据
    lonex1(1:nx, 1:ny, i) = lon;
    latex1(1:nx, 1:ny, i) = lat;
    years1(1:nx, 1:ny, i) = year;
    months1(1:nx, 1:ny, i) = month;
    days1(1:nx, 1:ny, i) = day;
    yr1(1:nx, i) = yr;
end

% ===== 保存结果 =====
save('PTC_track_ens_1981_2019.mat', 'lonex1', 'latex1', ...
    'years1', 'months1', 'days1', 'yr1');

disp('All data processed and saved!');
