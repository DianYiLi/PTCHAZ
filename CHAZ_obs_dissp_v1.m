clear;clc

%read obs
fpath='/Users/dianyili/Desktop/SBU/CHAZ/IBTrACS.NA.v04r00.nc';

bseason=ncread(fpath,'season');byr=bseason;
blat=ncread(fpath,'lat');
blon=ncread(fpath,'lon');blon(blon<0)=blon(blon<0)+360;
bwind=ncread(fpath,'wmo_wind');
bpres=ncread(fpath,'wmo_pres');
bstatus=ncread(fpath,'nature');
blandfall=ncread(fpath,'dist2land');
btime=ncread(fpath,'time');

x=find(bseason<1981 | bseason>2019);

blat(:,x)=[];
blon(:,x)=[];
bwind(:,x)=[];
bpres(:,x)=[];
blandfall(:,x)=[];
bstatus(:,:,x)=[];
btime(:,x)=[];
byr(x)=[];
for i=1:length(blon(:,1))
    for j=1:length(blon(1,:))
        if ~isnan(btime(i,j))
            hr=str2double(datestr(datenum('1858-11-17')+btime(i,j),'HH'));
            if hr~=0 && hr~=6 && hr~=12 && hr~=18
                blon(i,j)=nan;
                blat(i,j)=nan;
                bpres(i,j)=nan;
                bwind(i,j)=nan;
                blandfall(i,j)=nan;
                btime(i,j)=nan;

            end
        end
    end
end

blatex=nan(size(blat));blonex=blatex;bwindex=blatex;bpresex=blatex;byrex=nan(1,length(blat(1,:)));
for i=1:length(blat(1,:))
    for j=1:360
        if (bstatus(1,j,i)=='E' && bstatus(2,j,i)=='T') 
            if blandfall(j,i)<111*0.5/2;break;end

            bpresex(j:end,i)=bpres(j:end,i);
            bwindex(j:end,i)=bwind(j:end,i);
            blonex(j:end,i)=blon(j:end,i);
            blatex(j:end,i)=blat(j:end,i); 
            break
        end
    end
end

bpresex1=nan(1,length(blat(1,:)));bwindex1=bpresex1;blonex1=bpresex1;blatex1=bpresex1;
for i=1:length(blat(1,:))
    x=find(~isnan(blonex(:,i)));
    if ~isempty(x)
        bpresex1(i)=bpresex(x(end),i);
        bwindex1(i)=bwindex(x(end),i);
        blonex1(i)=blonex(x(end),i);
        blatex1(i)=blatex(x(end),i);
    end
end

%%



figure; hold on;
plot(bpresex1, blatex1, 'ok', 'DisplayName', 'All Points');
xlabel('Pressure (hPa)');
ylabel('Latitude (°N)');
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 14);
xlim([960 1020])
ylim([0 80])
%%

% ==== 1. 提取有效数据 ====
data_blue = bpresex1(~isnan(bpresex1));
data_blue = data_blue';

% ==== 2. 拟合整体 GEV 分布 ====
pd_gev = fitdist(data_blue, 'GeneralizedExtremeValue');
k_all = pd_gev.k;
sigma_all = pd_gev.sigma;
mu_all = pd_gev.mu;
save('gev_params_obs.mat', 'k_all', 'sigma_all', 'mu_all');

% ==== 3. 按纬度分段拟合局地 GEV ====
lat_list = 33:6:57;
dl = 6;
n_lat = length(lat_list);
p_vals = 900:1:1020;

% 初始化保存局地 GEV 参数
gev_params = nan(n_lat, 4);  % 第1列为纬度中心，2~4为k,sigma,mu

figure; hold on;
cdf_all = cdf(pd_gev, p_vals);
plot(p_vals, cdf_all*100, 'k', 'LineWidth', 3, 'DisplayName', 'NA');
pdata_sorted = sort(data_blue);
n = length(pdata_sorted);
F_emp = (1:n) ./ (n + 1);  % 或者用 Gringorten: ((1:n)-0.44)./(n+0.12)
plot(pdata_sorted, F_emp*100, ':', 'LineWidth', 3, 'color','k', 'HandleVisibility', 'off');

colors = lines(n_lat);
for i = 1:n_lat
    lat0 = lat_list(i);

    if i == 1
        idx = blatex1 < lat0 + dl/2;
        legend_label = sprintf('Lat=%.1f°N', lat0);
    elseif i == n_lat
        idx = blatex1 > lat0 - dl/2;
        legend_label = sprintf('Lat=%.1f°N', lat0);
    else
        idx = blatex1 > lat0 - dl/2 & blatex1 <= lat0 + dl/2;
        legend_label = sprintf('Lat=%.1f°N', lat0);
    end

    pdata = bpresex1(idx);
    pdata = pdata(~isnan(pdata));

    if numel(pdata) >= 20
        disp(numel(pdata))
        col_now = colors(i, :);
        pd_local = fitdist(pdata', 'gev');
        gev_params(i,:) = [lat0, pd_local.k, pd_local.sigma, pd_local.mu];

        % 绘图
        cdf_local = cdf(pd_local, p_vals);
        plot(p_vals, cdf_local*100, 'LineWidth', 3, 'color',col_now,'DisplayName', legend_label);

        pdata_sorted = sort(pdata);
        n = length(pdata_sorted);
        F_emp = (1:n) ./ (n + 1);  % 或者用 Gringorten: ((1:n)-0.44)./(n+0.12)
        plot(pdata_sorted, F_emp*100, ':', 'LineWidth', 3, 'color',col_now, 'HandleVisibility', 'off');
    end
end

xlabel('Pressure (hPa)');
ylabel('CDF(%)');
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 16);
xlim([950 1020]);
exportgraphics(gca, ['Obs_diss_cdf1.jpg'], 'Resolution', 1100);
close

% ==== 4. 纬度插值 ====
lat_interp = 0:1:90;  % 目标纬度
interp_params = nan(length(lat_interp), 4);  % 纬度，k, sigma, mu

valid_rows = ~any(isnan(gev_params(:,2:4)), 2);
lat_valid = gev_params(valid_rows,1);
k_valid = gev_params(valid_rows,2);
sigma_valid = gev_params(valid_rows,3);
mu_valid = gev_params(valid_rows,4);

for j = 1:length(lat_interp)
    lat_j = lat_interp(j);
    interp_params(j,1) = lat_j;

    if lat_j <= min(lat_valid)
        interp_params(j,2:4) = [k_valid(1), sigma_valid(1), mu_valid(1)];
    elseif lat_j >= max(lat_valid)
        interp_params(j,2:4) = [k_valid(end), sigma_valid(end), mu_valid(end)];
    else
        interp_params(j,2) = interp1(lat_valid, k_valid, lat_j);
        interp_params(j,3) = interp1(lat_valid, sigma_valid, lat_j);
        interp_params(j,4) = interp1(lat_valid, mu_valid, lat_j);
    end
end

% ==== 5. 保存结果 ====
% save('gev_params_interp_by_lat.mat', 'interp_params');

% 可选输出确认
% disp('插值后的 GEV 参数格式：');
% disp('   Lat     k       sigma     mu');
% disp(interp_params(1:5,:));  % 前几行预览
