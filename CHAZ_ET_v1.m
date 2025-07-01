%use  logistic regression with elastic net to get ET timing

clear;clc

fpath='/Users/dianyili/Desktop/SBU/CHAZ/atl_era5/';

%coefficient
% theta=[-3.26172713 2.34943926 -2.00465078  0.93136911 -0.39493953  0.84584208 -0.50226087...
%    0.43716437 -0.3671656   0.19401361 -0.35924508];%v1
% theta=[-2.70944337 2.2511956  -1.4001523   0.82228464 -0.15269258  0.84556124 -0.50513599...
%    0.22858954 -0.20680489 -0.11013684 -0.32325854];%v11
theta=[-2.52806748 2.00320622 -1.37477621  0.76486714 -0.34419957  0.76759218 -0.51293975...
   0.21562304 -0.14949956  0.18873102 -0.37461071];%1981-2019

% std and mean
mu=[28.7480 50.6521 9.7859 54.5131 6.4690 0.9064 0.1145 0.2693 -1.4074 0.1469];
delta=[12.3493 26.0103 5.6671 24.0853 3.7163 1.4236 9.7353 1.3118 4.9011 1.3516];
% mu=[28.7905 50.4720 9.7852 54.4842 6.4689 0.9097 0.1116 0.2681 -1.4136 0.1465];
% delta=[12.3469 25.9428 5.6558 24.0989 3.7128 1.4229 9.7200 1.3077 4.9049 1.3491];%1981-2018

% load('logistic_nomarlize.mat')
% mu=xm;
% delta=xs;

% used for calculate mu and delta
% lati=nan(1,125);
% p0i=nan(1,125,40);
% mwsi=nan(1,125,40);
% shri=nan(1,125);
% piwi=nan(1,125);
% tsi=nan(1,125);
% dlati=nan(1,125);
% dp0i=nan(1,125,40);
% dmwsi=nan(1,125,40);
% dshri=nan(1,125);
% dpiwi=nan(1,125);
% dtsi=nan(1,125);

for i=1:40
    disp(['Ens ',num2str(i)])
    if i-1<10
        fname=['atl_2019_2ens00',num2str(i-1),'_pre.nc'];
    else
        fname=['atl_2019_2ens0',num2str(i-1),'_pre.nc'];
    end
    ff=[fpath,fname];

    lon=ncread(ff,'longitude');
    lat=ncread(ff,'latitude');
    mws=ncread(ff,'Mwspd');
    ushear=ncread(ff,'ushear');
    vshear=ncread(ff,'vshear');
    piw=ncread(ff,'PIWspd')*0.5144;
    yr=ncread(ff,'year');

    
    [a,b]=size(lat);
    lat12=nan(a,b);
    lon12=nan(a,b);
    dlat=nan(a,b);
    dpiw=nan(a,b);
    dshr=nan(a,b);
    dts=nan(a,b);
    dp0=nan(a,b,40);
    dmws=nan(a,b,40);
    hx=nan(a,b);
    ext=nan(1,a);

   
    re=6371;
    for j=2:b-1
        lat12(:,j)=(lat(:,j+1)-lat(:,j-1))/12;
        lon12(:,j)=(lon(:,j+1)-lon(:,j-1))/12;
    end
    lat12(:,b)=lat12(:,1);
    lon12(:,b)=lon12(:,1);
    dy=2*pi*re*lat12/360;
    dx=2*pi*re*cosd(lat).*lon12/360;
    dxy=sqrt(dx.^2+dy.^2);%km/h
    ts=dxy/3.6;
    shr=sqrt(ushear.^2+vshear.^2);

    p0=23.286-0.483*(mws-ts)-((mws-ts)/24.254).^2-12.587*0.49-0.483*lat+1014.25;%Knaff and Zehr (2007)

    for j=3:b
        dlat(:,j)=lat(:,j)-lat(:,j-2);
        dpiw(:,j)=piw(:,j)-piw(:,j-2);
        dshr(:,j)=shr(:,j)-shr(:,j-2);
        dts(:,j)=ts(:,j)-ts(:,j-2);
        dp0(:,j,:)=p0(:,j,:)-p0(:,j-2,:);
        dmws(:,j,:)=mws(:,j,:)-mws(:,j-2,:);
    end
    
    % used for calculate mu and delta
    % idx=find(yr<1981|yr>2018);
    % lat(idx,:)=nan;
    % p0(idx,:,:)=nan;
    % mws(idx,:,:)=nan;
    % shr(idx,:)=nan;
    % piw(idx,:)=nan;
    % ts(idx,:)=nan;
    % dlat(idx,:)=nan;
    % dp0(idx,:,:)=nan;
    % dmws(idx,:,:)=nan;
    % dshr(idx,:)=nan;
    % dpiw(idx,:)=nan;
    % dts(idx,:)=nan;
    % 
    % lati=[lati;lat];
    % p0i=[p0i;p0];
    % mwsi=[mwsi;mws];
    % shri=[shri;shr];
    % piwi=[piwi;piw];
    % tsi=[tsi;ts];
    % dlati=[dlati;dlat];
    % dp0i=[dp0i;dp0];
    % dmwsi=[dmwsi;dmws];
    % dshri=[dshri;dshr];
    % dpiwi=[dpiwi;dpiw];
    % dtsi=[dtsi;dts];

    latn=(lat-mu(1))/delta(1);
    mwsn=(mws-mu(2))/delta(2);

    shrn=(shr-mu(3))/delta(3);
    piwn=(piw-mu(4))/delta(4);
    tsn=(ts-mu(5))/delta(5);
    dlatn=(dlat-mu(6))/delta(6);
    dmwsn=(dmws-mu(7))/delta(7);

    dshrn=(dshr-mu(8))/delta(8);
    dpiwn=(dpiw-mu(9))/delta(9);
    dtsn=(dts-mu(10))/delta(10);


    for j=1:a
        for k=1:b
            for en=1:40
                x=[1 latn(j,k) mwsn(j,k,en) shrn(j,k) piwn(j,k) tsn(j,k)...
                    dlatn(j,k) dmwsn(j,k,en) dshrn(j,k) dpiwn(j,k) dtsn(j,k)];
                hx(j,k,en)=1/(1+exp(-theta*x'));
            end
        end
    end
    ex=hx;
    lev=0.5;
    ex(hx>=lev)=1;
    ex(hx<lev)=0;

    if i-1<10
        fname_et=['atl_2019_2ens00',num2str(i-1),'_ETv1.nc'];
    else
        fname_et=['atl_2019_2ens0',num2str(i-1),'_ETv1.nc'];
    end
    ff_et=[fpath,fname_et];
    copyfile(ff,ff_et)

    nccreate(ff_et,'ex','Dimensions',{"stormID",a,"lifelength",b,"ensembleNum",40})
    ncwrite(ff_et,'ex',ex);
    ncwriteatt(ff_et,'ex','long_name','If extratropic')

    nccreate(ff_et,'p0','Dimensions',{"stormID",a,"lifelength",b,"ensembleNum",40})
    ncwrite(ff_et,'p0',p0);
    ncwriteatt(ff_et,'p0','long_name','pressure (hPa), calculated using W-P relationship')

end

% mu(1)=mean(lat(:),'omitnan');
% mu(2)=mean(mws(:),'omitnan');
% mu(3)=mean(shr(:),'omitnan');
% mu(4)=mean(piw(:),'omitnan');
% mu(5)=mean(ts(:),'omitnan');
% mu(6)=mean(dlat(:),'omitnan');
% mu(7)=mean(dmws(:),'omitnan');
% mu(8)=mean(dshr(:),'omitnan');
% mu(9)=mean(dpiw(:),'omitnan');
% mu(10)=mean(dts(:),'omitnan');
% 
% 
% delta(1)=std(lat(:),'omitnan');
% delta(2)=std(mws(:),'omitnan');
% delta(3)=std(shr(:),'omitnan');
% delta(4)=std(piw(:),'omitnan');
% delta(5)=std(ts(:),'omitnan');
% delta(6)=std(dlat(:),'omitnan');
% delta(7)=std(dmws(:),'omitnan');
% delta(8)=std(dshr(:),'omitnan');
% delta(9)=std(dpiw(:),'omitnan');
% delta(10)=std(dts(:),'omitnan');




