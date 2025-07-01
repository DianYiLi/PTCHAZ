clear; clc

fpath = '/Users/dianyili/Desktop/SBU/CHAZ/atl_era5/';

num_ens = 40; 
max_x = 700;  
max_y = 125;  

lonex1   = nan(max_x, max_y, num_ens);
latex1   = nan(max_x, max_y, num_ens);
years1   = nan(max_x, max_y, num_ens);
months1  = nan(max_x, max_y, num_ens);
days1    = nan(max_x, max_y, num_ens);
yr1      = nan(max_x, num_ens);  

for i = 1:num_ens

    disp(['Processing Ensemble ', num2str(i)])

    if i-1 < 10
        fname = ['atl_2019_2ens00', num2str(i-1), '_ETv1.nc'];
    else
        fname = ['atl_2019_2ens0', num2str(i-1), '_ETv1.nc'];
    end
    ff = [fpath, fname];

    lon = ncread(ff, 'longitude');
    lat = ncread(ff, 'latitude');
    yr = ncread(ff, 'year');
    ld = ncread(ff, 'ldmask');
    tm = ncread(ff, 'time');  % days since 1950-01-01
    tm(tm < 0) = nan; 
    startDate = datetime(1950, 1, 1);
    dates = startDate + days(tm);
    [year, month, day] = ymd(dates);

    valid_idx = yr >= 1981 & yr <= 2019;
    lat = lat(valid_idx, :);
    lon = lon(valid_idx, :);
    ld = ld(valid_idx, :);
    year = year(valid_idx, :);
    month = month(valid_idx, :);
    day = day(valid_idx, :);
    yr = yr(valid_idx);

    [nx, ny] = size(lat);

    lonex1(1:nx, 1:ny, i) = lon;
    latex1(1:nx, 1:ny, i) = lat;
    years1(1:nx, 1:ny, i) = year;
    months1(1:nx, 1:ny, i) = month;
    days1(1:nx, 1:ny, i) = day;
    yr1(1:nx, i) = yr;
end

save('PTC_track_ens_1981_2019.mat', 'lonex1', 'latex1', ...
    'years1', 'months1', 'days1', 'yr1');

disp('All data processed and saved!');
