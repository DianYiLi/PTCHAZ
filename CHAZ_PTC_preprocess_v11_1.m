clear; clc

% === 路径参数 ===
epath = '/glade/campaign/collections/rda/data/d633001/e5.moda.an.pl';
load('PTC_track_ens_1981_2019.mat')  % 预存轨迹数据

% === 基本设置 ===
re = 6371;  % 地球半径，单位 km
b = size(lonex1, 1);  % 气旋数
a = size(lonex1, 2);  % 时间步
%ens = 1;              % 指定集合成员
yr_tp = 0;

% === 要提取的变量名和对应 NetCDF 后缀（你需要根据数据定义这些）===
varNames = {'U', 'V', 'PV'};  % 举例：你可以添加更多变量
fileSuffixes = {'128_131_u.ll025uv', '128_132_v.ll025uv', '128_060_pv.ll025sc'};
outputVars = {'u1', 'v1','pv1'};  % 最终输出变量名

% 初始化输出
for vIdx = 1:length(varNames)
    eval([outputVars{vIdx}, ' = nan(a, b, 37);']);  % 假设 37 层
end

disp('prepare variables');

for i = 1:b
    if ~isnan(yr1(i, ens))
        % 只在年变动时重新读取数据
        if yr_tp ~= yr1(i, ens)
            data = struct();
            for vIdx = 1:length(varNames)
                fname = sprintf('%s/%d/e5.moda.an.pl.%s.%d010100_%d120100.nc', ...
                    epath, yr1(i, ens), fileSuffixes{vIdx}, yr1(i, ens), yr1(i, ens));
                data.(varNames{vIdx}) = ncread(fname, varNames{vIdx});
            end
            yr_tp = yr1(i, ens);

            lat = ncread(fname, 'latitude');
            lon = ncread(fname, 'longitude');
            [xx1, yy1] = meshgrid(lon, lat);
        end

        for j = 1:a
            if ~isnan(latex1(i, j, ens))
                dd = days1(i, j, ens);
                mm = months1(i, j, ens);

                interpolatedData = struct();
                for vIdx = 1:length(varNames)
                    varData = data.(varNames{vIdx});

                    if dd < 15 && mm ~= 1
                        interpolatedData.(varNames{vIdx}) = ...
                            varData(:,:,:,mm-1) + (varData(:,:,:,mm) - varData(:,:,:,mm-1)) * (dd+15)/30;
                    elseif dd < 15 && mm == 1
                        fnamePrev = sprintf('%s/%d/e5.moda.an.pl.%s.%d010100_%d120100.nc', ...
                            epath, yr1(i, ens)-1, fileSuffixes{vIdx}, yr1(i, ens)-1, yr1(i, ens)-1);
                        prevData = ncread(fnamePrev, varNames{vIdx});
                        interpolatedData.(varNames{vIdx}) = ...
                            prevData(:,:,:,12) + (varData(:,:,:,mm) - prevData(:,:,:,12)) * (dd+15)/30;
                    elseif dd >= 15 && mm ~= 12
                        interpolatedData.(varNames{vIdx}) = ...
                            varData(:,:,:,mm) + (varData(:,:,:,mm+1) - varData(:,:,:,mm)) * (dd-15)/30;
                    else
                        fnameNext = sprintf('%s/%d/e5.moda.an.pl.%s.%d010100_%d120100.nc', ...
                            epath, yr1(i, ens)+1, fileSuffixes{vIdx}, yr1(i, ens)+1, yr1(i, ens)+1);
                        nextData = ncread(fnameNext, varNames{vIdx});
                        interpolatedData.(varNames{vIdx}) = ...
                            varData(:,:,:,mm) + (nextData(:,:,:,1) - varData(:,:,:,mm)) * (dd-15)/30;
                    end
                end

                % 空间半径内平均
                dx = 2 * pi * re * cosd(latex1(i, j, ens)) * (xx1 - lonex1(i, j, ens)) / 360;
                dy = 2 * pi * re * (latex1(i, j, ens) - yy1) / 360;
                dxy = sqrt(dx.^2 + dy.^2)';

                for k = 1:37
                    for vIdx = 1:length(varNames)
                        varLayer = interpolatedData.(varNames{vIdx})(:,:,k);
                        eval([outputVars{vIdx}, '(j, i, k) = mean(varLayer(dxy <= 600), ''omitnan'');']);
                    end
                end
            end
        end
        disp(['Finished case ', num2str(i)])
    end
end

% 保存结果（可按需修改文件名）
save(['CHAZ_PTC_interpvars1_ens', num2str(ens), '_1981_2019.mat'], ...
    outputVars{:}, 'lonex1', 'latex1', 'yr1');

