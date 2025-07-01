clear;clc
%%
disp('Prepare lat, translation speed, pressure, wind, EX')

fpath='/Users/dianyili/Desktop/SBU/CHAZ/IBTrACS.NA.v04r00.nc';

bseason=ncread(fpath,'season');
blat=ncread(fpath,'lat');
blon=ncread(fpath,'lon');
wind=ncread(fpath,'wmo_wind');
pres=ncread(fpath,'wmo_pres');
status=ncread(fpath,'nature');
landfall=ncread(fpath,'landfall');
time=ncread(fpath,'time');

% yr1=1981;yr2=2010;%used
yr1=1981;yr2=2019;%CHAZ is 1981-2019, extend ibtracs is 1980-2018

x=find(bseason<yr1 | bseason>yr2);
blat(:,x)=[];
blon(:,x)=[];
wind(:,x)=[];
pres(:,x)=[];
time(:,x)=[];
landfall(:,x)=[];
status(:,:,x)=[];
bseason(x)=[];

blon1=nan(size(blon));
blat1=nan(size(blon));
time1=nan(size(blon));
wind1=nan(size(blon));
pres1=nan(size(blon));
landfall1=nan(size(blon));
status1=nan(size(status));

[a,b]=size(blon);
for j=1:b
    n=1;
    for i=1:a
        if ~isnan(time(i,j))
            hr=str2double(datestr(datenum('1858-11-17')+time(i,j),'HH'));
            if hr==0 || hr==6 || hr==12 || hr==18
                blon1(n,j)=blon(i,j);
                blat1(n,j)=blat(i,j);
                wind1(n,j)=wind(i,j);
                pres1(n,j)=pres(i,j);
                time1(n,j)=time(i,j);
                landfall1(n,j)=landfall(i,j);% also remove ET over land?
                status1(:,n,j)=status(:,i,j);
                n=n+1;
            end
        end
    end
end

et=nan(size(blon));
for i=1:b
    for j=1:a
        if (status1(1,j,i)=='E' && status1(2,j,i)=='T')
            et(j,i)=1;
        else
            et(j,i)=0;
        end
    end
end

lat12=nan(a,b);
lon12=nan(a,b);
re=6371;
for j=2:a-1
    lat12(j,:)=(blat1(j+1,:)-blat1(j-1,:))/12;
    lon12(j,:)=(blon1(j+1,:)-blon1(j-1,:))/12;
end
lat12(a,:)=lat12(1,:);
lon12(a,:)=lon12(1,:);
dy=2*pi*re*lat12/360;
dx=2*pi*re*cosd(blat1).*lon12/360;
dxy=sqrt(dx.^2+dy.^2);%km/h
ts=dxy/3.6;

%%
disp('prepare PI')
pipath='/Users/dianyili/Desktop/SBU/CHAZ/PI';
% piname=dir(pipath);
% piname1=piname(34:72);%1981-2019
pi1=nan(a,b);
for i=1:b
    load([pipath,'/VelMax_ERA5_PI_',num2str(bseason(i)),'_2deg_CMIP6.mat'])
    X(X>180)=X(X>180)-360;
    [xx,yy]=meshgrid(X,Y);
    for j=1:a
        if ~isnan(blat1(j,i))
            dd=str2double(datestr(datenum('1858-11-17')+time1(j,i),'dd'));
            mm=str2double(datestr(datenum('1858-11-17')+time1(j,i),'mm'));
            if dd<15 && mm~=1
                piw=VmaxI(:,:,mm-1)+(VmaxI(:,:,mm)-VmaxI(:,:,mm-1))*(dd+15)/30;
            elseif dd<15 && mm==1
                VmaxI1=load([pipath,'/VelMax_ERA5_PI_',num2str(bseason(i)-1),'_2deg_CMIP6.mat'],'VmaxI');
                VmaxI1=VmaxI1.VmaxI;
                piw=VmaxI1(:,:,12)+(VmaxI(:,:,mm)-VmaxI1(:,:,12))*(dd+15)/30;
            elseif dd>=15 && mm~=12
                piw=VmaxI(:,:,mm)+(VmaxI(:,:,mm+1)-VmaxI(:,:,mm))*(dd-15)/30;
            else
                VmaxI1=load([pipath,'/VelMax_ERA5_PI_',num2str(bseason(i)+1),'_2deg_CMIP6.mat'],'VmaxI');
                VmaxI1=VmaxI1.VmaxI;
                piw=VmaxI(:,:,mm)+(VmaxI1(:,:,1)-VmaxI(:,:,mm))*(dd-15)/30;
            end
            ddpi=(blat1(j,i)-yy).^2+(blon1(j,i)-xx).^2;
            [x1,y1]=find(ddpi==min(ddpi,[],'all','omitnan'));
            x1=x1(1);y1=y1(1);
            pi1(j,i)=piw(x1,y1);
            % dy=2*pi*re*(blat1(j,i)-yy)/360;
            % dx=2*pi*re*cosd((blat1(j,i)+yy)/2).*(blon1(j,i)-xx)/360;
            % dxy=sqrt(dx.^2+dy.^2)';
            % pi1(j,i)=mean(piw(dxy<=500),'omitnan');%CHAZ
        end
    end
    % disp(i)
end

%%
disp('Prepare SHR')
shrpath='/Users/dianyili/Desktop/SBU/CHAZ/era5wind';
shr1=nan(a,b);
for i=1:b
    u=ncread([shrpath,'/u_component_of_wind_',num2str(bseason(i)),'.nc'],'u');
    v=ncread([shrpath,'/v_component_of_wind_',num2str(bseason(i)),'.nc'],'v');
    shr=squeeze(sqrt((u(:,:,1,:)-u(:,:,4,:)).^2+(v(:,:,1,:)-v(:,:,4,:)).^2));
    % level=ncread([shrpath,'/v_component_of_wind_',num2str(bseason(i)),'.nc'],'level');
    lat=ncread([shrpath,'/u_component_of_wind_',num2str(bseason(i)),'.nc'],'latitude');
    lon=ncread([shrpath,'/u_component_of_wind_',num2str(bseason(i)),'.nc'],'longitude');
    lon(lon>180)=lon(lon>180)-360;
    [xx1,yy1]=meshgrid(lon,lat);
    for j=1:a
        if ~isnan(blat1(j,i))
            dd=str2double(datestr(datenum('1858-11-17')+time1(j,i),'dd'));
            mm=str2double(datestr(datenum('1858-11-17')+time1(j,i),'mm'));
            if dd<15 && mm~=1
                shrm=shr(:,:,mm-1)+(shr(:,:,mm)-shr(:,:,mm-1))*(dd+15)/30;
            elseif dd<15 && mm==1
                u1=ncread([shrpath,'/u_component_of_wind_',num2str(bseason(i)-1),'.nc'],'u');
                v1=ncread([shrpath,'/v_component_of_wind_',num2str(bseason(i)-1),'.nc'],'v');
                shrm1=squeeze(sqrt((u1(:,:,1,:)-u1(:,:,4,:)).^2+(v1(:,:,1,:)-v1(:,:,4,:)).^2));
                shrm=shrm1(:,:,12)+(shr(:,:,mm)-shrm1(:,:,12))*(dd+15)/30;
            elseif dd>=15 && mm~=12
                shrm=shr(:,:,mm)+(shr(:,:,mm+1)-shr(:,:,mm))*(dd-15)/30;
            else
                u1=ncread([shrpath,'/u_component_of_wind_',num2str(bseason(i)+1),'.nc'],'u');
                v1=ncread([shrpath,'/v_component_of_wind_',num2str(bseason(i)+1),'.nc'],'v');
                shrm1=squeeze(sqrt((u1(:,:,1,:)-u1(:,:,4,:)).^2+(v1(:,:,1,:)-v1(:,:,4,:)).^2));
                shrm=shr(:,:,mm)+(shrm1(:,:,1)-shr(:,:,mm))*(dd-15)/30;
            end
            
            dy=2*pi*re*(blat1(j,i)-yy1)/360;
            dx=2*pi*re*cosd((blat1(j,i)+yy1)/2).*(blon1(j,i)-xx1)/360;
            dxy=sqrt(dx.^2+dy.^2)';
            
            % shr_close=shrm(dxy==min(dxy,[],'all','omitnan'));%Melanie
            % shr1(j,i)=mean(shr_close);
            shr1(j,i)=mean(shrm(dxy>=200 &dxy<=800),'omitnan');%CHAZ
        end
    end
    % disp(i)
end

%%
disp('Calculate delta')

dlat=nan(a,b);
dts=dlat;dpres=dlat;dwind=dlat;dpi=dlat;dshr=dlat;
for j=3:a
    dlat(j,:)=blat1(j,:)-blat1(j-2,:);
    dts(j,:)=ts(j,:)-ts(j-2,:);
    dpres(j,:)=pres1(j,:)-pres1(j-2,:);
    dwind(j,:)=wind1(j,:)-wind1(j-2,:);
    dpi(j,:)=pi1(j,:)-pi1(j-2,:);
    dshr(j,:)=shr1(j,:)-shr1(j-2,:);
end

save('logistic_data.mat')
%%
disp('split training/testing')

total_cyclones = b;  % 假设 b 是气旋总数
testing_ratio = 0.2;
num_testing = floor(total_cyclones * testing_ratio);

rng(23)  % 固定随机种子以保证可重复性
test_indices = randsample(total_cyclones, num_testing);

test_mask = zeros(total_cyclones, 1);
test_mask(test_indices) = 1;

n = 1; m = 1;
for i = 1:total_cyclones
    if test_mask(i) == 0
        et_train(:,n)    = et(:,i);
        lat_train(:,n)   = blat1(:,i);
        p0_train(:,n)    = pres1(:,i);
        w_train(:,n)     = wind1(:,i);
        pi_train(:,n)    = pi1(:,i);
        shr_train(:,n)   = shr1(:,i);
        ts_train(:,n)    = ts(:,i);
        dlat_train(:,n)  = dlat(:,i);
        dp0_train(:,n)   = dpres(:,i);
        dw_train(:,n)    = dwind(:,i);
        dpi_train(:,n)   = dpi(:,i);
        dshr_train(:,n)  = dshr(:,i);
        dts_train(:,n)   = dts(:,i);
        n = n + 1;
    else
        et_test(:,m)     = et(:,i);
        lat_test(:,m)    = blat1(:,i);
        p0_test(:,m)     = pres1(:,i);
        w_test(:,m)      = wind1(:,i);
        pi_test(:,m)     = pi1(:,i);
        shr_test(:,m)    = shr1(:,i);
        ts_test(:,m)     = ts(:,i);
        dlat_test(:,m)   = dlat(:,i);
        dp0_test(:,m)    = dpres(:,i);
        dw_test(:,m)     = dwind(:,i);
        dpi_test(:,m)    = dpi(:,i);
        dshr_test(:,m)   = dshr(:,i);
        dts_test(:,m)    = dts(:,i);
        m = m + 1;
    end
end

        
%%

disp('Remove nan')

dts_train1=dts_train;dp0_train1=dp0_train;dpi_train1=dpi_train;
dts_test1=dts_test;dp0_test1=dp0_test;dpi_test1=dpi_test;

rmnan_test=isnan(dts_test1)|isnan(dp0_test1)|isnan(dpi_test1);
rmnan_train=isnan(dts_train1)|isnan(dp0_train1)|isnan(dpi_train1);

et_train(rmnan_train)=[];
lat_train(rmnan_train)=[];
p0_train(rmnan_train)=[];
w_train(rmnan_train)=[];
pi_train(rmnan_train)=[];
shr_train(rmnan_train)=[];
ts_train(rmnan_train)=[];
dlat_train(rmnan_train)=[];
dp0_train(rmnan_train)=[];
dw_train(rmnan_train)=[];
dpi_train(rmnan_train)=[];
dshr_train(rmnan_train)=[];
dts_train(rmnan_train)=[];

et_test(rmnan_test)=[];
lat_test(rmnan_test)=[];
p0_test(rmnan_test)=[];
w_test(rmnan_test)=[];
pi_test(rmnan_test)=[];
shr_test(rmnan_test)=[];
ts_test(rmnan_test)=[];
dlat_test(rmnan_test)=[];
dp0_test(rmnan_test)=[];
dw_test(rmnan_test)=[];
dpi_test(rmnan_test)=[];
dshr_test(rmnan_test)=[];
dts_test(rmnan_test)=[];



%%
%Melanie's 
% x_train=[lat_train;p0_train;shr_train;pi_train;ts_train;...
%     dlat_train;dp0_train;dshr_train;dpi_train;dts_train]';
% y_train=et_train;
% x_test=[lat_test;p0_test;shr_test;pi_test;ts_test;...
%     dlat_test;dp0_test;dshr_test;dpi_test;dts_test]';
% y_test=et_test;

%new
x_train=[lat_train;w_train;shr_train;pi_train;ts_train;...
    dlat_train;dw_train;dshr_train;dpi_train;dts_train]';
y_train=et_train;
x_test=[lat_test;w_test;shr_test;pi_test;ts_test;...
    dlat_test;dw_test;dshr_test;dpi_test;dts_test]';
y_test=et_test;

xm=mean(x_train);
xs=std(x_train);
save('logistic_nomarlize.mat','xm','xs')
% pre_m=mean(pre);
% pre_s=std(pre);
% for i=1:10
%     pre_n(:,i)=(pre(:,i)-pre_m(i))/pre_s(i);
% end
% 
% [B,FitInfo] = lasso(pre_n,et);

%%
disp('Write to nc')

ff_pre='/Users/dianyili/Desktop/SBU/CHAZ/pre_logistic_v1_1981-2019.nc';

nccreate(ff_pre,'x_train','Dimensions',{"training data length",length(y_train),'predictor no',10})
ncwrite(ff_pre,'x_train',x_train);
ncwriteatt(ff_pre,'x_train','long_name','predictors in training set')

nccreate(ff_pre,'y_train','Dimensions',{"training data length",length(y_train)})
ncwrite(ff_pre,'y_train',y_train);
ncwriteatt(ff_pre,'y_train','long_name','if extratropical in training set')

nccreate(ff_pre,'x_test','Dimensions',{"testing data length",length(y_test),'predictor no',10})
ncwrite(ff_pre,'x_test',x_test);
ncwriteatt(ff_pre,'x_test','long_name','predictors in testing set')

nccreate(ff_pre,'y_test','Dimensions',{"testing data length",length(y_test)})
ncwrite(ff_pre,'y_test',y_test);
ncwriteatt(ff_pre,'y_test','long_name','if extratropical in testing set')
