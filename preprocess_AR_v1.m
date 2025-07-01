clear;clc

%read obs
fpath='/glade/u/home/dyli/IBTrACS.NA.v04r00.nc';

bseason=ncread(fpath,'season');byr=bseason;
blat=ncread(fpath,'lat');
blon=ncread(fpath,'lon');%blon(blon<0)=blon(blon<0)+360;
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
            if blandfall(j,i)<111*0.25/2;break;end

            bpresex(j:end,i)=bpres(j:end,i);
            bwindex(j:end,i)=bwind(j:end,i);
            blonex(j:end,i)=blon(j:end,i);
            blatex(j:end,i)=blat(j:end,i); 
            break
        end
    end
end


[a,b]=size(blonex);


re=6371;

%%
bseason=byr;
timeex=btime;

disp('prepare PI')
pipath='/glade/work/dyli/derecho/CHAZ/PI';
% piname=dir(pipath);
% piname1=piname(34:72);%1981-2019
pi1=nan(a,b);
for i=1:b
    load([pipath,'/VelMax_ERA5_PI_',num2str(bseason(i)),'_2deg_CMIP6.mat'])
    X(X>180)=X(X>180)-360;
    [xx,yy]=meshgrid(X,Y);
    for j=1:a
        if ~isnan(blatex(j,i))
            dd=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'dd'));
            mm=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'mm'));
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
            %ddpi=(blatex(j,i)-yy).^2+(blonex(j,i)-xx).^2;
            %[x1,y1]=find(ddpi==min(ddpi,[],'all','omitnan'));
            %x1=x1(1);y1=y1(1);
	    %pi1(j,i)=piw(x1,y1);
	    dx=2*pi*re*cosd(blatex(j,i)).*(blonex(j,i)-xx)/360;
	    dy=2*pi*re*(blatex(j,i)-yy)/360;
	    dxy=sqrt(dx.^2+dy.^2);
	    pi1(j,i)=mean(piw(dxy<=600),'omitnan');
	    
	    %disp(size(piw))
        end
    end
    disp(i)
end


%%
disp('prepare geopotential height')
epath='/glade/campaign/collections/rda/data/d633001/e5.moda.an.pl';
% zpath='/Users/dianyili/Desktop/SBU/';
z1=nan(a,b,37);
gradz_ns1=z1;gradz_ew1=z1;
zzmax1000=z1;zzmin1000=z1;zzmax_x1000=z1;zzmin_x1000=z1;z_avg_below_5_1000=z1;z_avg_above_95_1000=z1;
zzmax1500=z1;zzmin1500=z1;zzmax_x1500=z1;zzmin_x1500=z1;z_avg_below_5_1500=z1;z_avg_above_95_1500=z1;
%z1w5=z1;z1w10=z1;z1e5=z1;z1e10=z1;
for i=1:b
    zname=[epath,'/',num2str(bseason(i)),'/e5.moda.an.pl.128_129_z.ll025sc.',...
        num2str(bseason(i)),'010100_',num2str(bseason(i)),'120100.nc'];
    % zname=[zpath,'e5.moda.an.pl.128_129_z.ll025sc.1981010100_1981120100.nc'];
    z=ncread(zname,'Z')/9.81;
    % level=ncread(zname,'level');
    lat=ncread(zname,'latitude');
    lon=ncread(zname,'longitude');
    lon(lon>180)=lon(lon>180)-360;
    [xx1,yy1]=meshgrid(lon,lat);
    dlat = abs(gradient(lat) * (pi/180) * re);             % Change in latitude (y) in km
    dlon = abs(gradient(lon) * (pi/180) .* cosd(lat') * re); % Change in longitude (x) in km
    for j=1:a
        if ~isnan(blatex(j,i))
            dd=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'dd'));
            mm=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'mm'));
            if dd<15 && mm~=1
                zz=z(:,:,:,mm-1)+(z(:,:,:,mm)-z(:,:,:,mm-1))*(dd+15)/30;
            elseif dd<15 && mm==1
                zname=[epath,'/',num2str(bseason(i)-1),'/e5.moda.an.pl.128_129_z.ll025sc.',...
                    num2str(bseason(i)-1),'010100_',num2str(bseason(i)-1),'120100.nc'];
                zz2=ncread(zname,'Z')/9.81;
                zz=z(:,:,:,12)+(z(:,:,:,mm)-zz2(:,:,:,12))*(dd+15)/30;
            elseif dd>=15 && mm~=12
                zz=z(:,:,:,mm)+(z(:,:,:,mm+1)-z(:,:,:,mm))*(dd-15)/30;
            else
                zname=[epath,'/',num2str(bseason(i)+1),'/e5.moda.an.pl.128_129_z.ll025sc.',...
                    num2str(bseason(i)+1),'010100_',num2str(bseason(i)+1),'120100.nc'];
                zz2=ncread(zname,'Z')/9.81;
                zz=z(:,:,:,mm)+(zz2(:,:,:,1)-z(:,:,:,mm))*(dd-15)/30;
            end
	    
	    dx=2*pi*re*cosd(blatex(j,i)).*(xx1-blonex(j,i))/360;
            dy=2*pi*re*(blatex(j,i)-yy1)/360;
	    dxt=dx';
            dxy=sqrt(dx.^2+dy.^2)';
           
	    gradZ_x_km = zeros(size(zz)); % To store gradient in x direction in km
            gradZ_y_km = zeros(size(zz)); % To store gradient in y direction in km
	    for k = 1:size(zz, 3)
   		 % Calculate gradient for the k-th layer
    		[gradZ_x, gradZ_y] = gradient(zz(:,:,k));
    
    		% Divide by spatial distance to get gradient in km
    		gradZ_x_km(:,:,k) = gradZ_x ./ dlat';
   		gradZ_y_km(:,:,k) = gradZ_y ./ dlon;
		
	    end	


            %ddzz=(blatex(j,i)-yy1).^2+(blonex(j,i)-xx1).^2;
            %[x1,y1]=find(ddzz==min(ddzz,[],'all','omitnan'));
            %x1=x1(1);y1=y1(1);
	    %disp(size(zz))
            %z1(j,i,:)=zz(x1,y1,:);
	    %gradz_x1(j,i,:)=gradZ_x_km(x1,y1,:);
	    %gradz_y1(j,i,:)=gradZ_y_km(x1,y1,:);

	    %dx=2*pi*re*cosd(blatex(j,i))*(blonex(j,i)-lon)/360;
            %dy=2*pi*re*(blatex(j,i)-lat)/360;
            %dxy=sqrt(dx.^2+dy.^2);
	    for k=1:size(zz, 3)
	        zzk=zz(:,:,k);
                z1(j,i,k)=mean(zzk(dxy<=600),'omitnan');
		gradz_x1k=gradZ_x_km(:,:,k);
		gradz_ns1(j,i,k)=mean(gradz_x1k(dxy<=600),'omitnan');
		gradz_y1k=gradZ_y_km(:,:,k);
                gradz_ew1(j,i,k)=mean(gradz_y1k(dxy<=600),'omitnan');

		zzmax1000(j,i,k)=max(zzk(dxy<1000),[],'omitnan');
                zzmin1000(j,i,k)=min(zzk(dxy<1000),[],'omitnan');
                [x,y]=find(zzk==zzmax1000(j,i,k) & dxy<1000);
                zzmax_x1000(j,i,k)=dxt(x(1),y(1));
                [x,y]=find(zzk==zzmin1000(j,i,k) & dxy<1000);
                zzmin_x1000(j,i,k)=dxt(x(1),y(1));

		zzmax1500(j,i,k)=max(zzk(dxy<1500),[],'omitnan');
		zzmin1500(j,i,k)=min(zzk(dxy<1500),[],'omitnan');
		[x,y]=find(zzk==zzmax1500(j,i,k) & dxy<1500);
		zzmax_x1500(j,i,k)=dxt(x(1),y(1));
		[x,y]=find(zzk==zzmin1500(j,i,k) & dxy<1500);
		zzmin_x1500(j,i,k)=dxt(x(1),y(1));
		
		zzk1=zzk(dxy<1000 & ~isnan(zzk));
                percentile_95=prctile(zzk1,95);
                z_above_95=zzk1(zzk1>percentile_95);
                if ~isempty(z_above_95) % Check if there are elements above 95th percentile
                        z_avg_above_95_1000(j,i,k) = mean(z_above_95);
                else
                        z_avg_above_95_1000(j,i,k) = NaN; % Handle case where no values are above 95th percentile
                end
                percentile_5=prctile(zzk1,5);
                z_below_5=zzk1(zzk1<percentile_5);
                if ~isempty(z_below_5) % Check if there are elements below 5th percentile
                        z_avg_below_5_1000(j,i,k) = mean(z_below_5,'omitnan');
                else
                        z_avg_below_5_1000(j,i,k) = NaN; % Handle case where no values are above 95th percentile
                end

		zzk1=zzk(dxy<1500 & ~isnan(zzk));
		percentile_95=prctile(zzk1,95);
		z_above_95=zzk1(zzk1>percentile_95);
		if ~isempty(z_above_95) % Check if there are elements above 95th percentile
   			z_avg_above_95_1500(j,i,k) = mean(z_above_95);
		else
    			z_avg_above_95_1500(j,i,k) = NaN; % Handle case where no values are above 95th percentile
		end
		percentile_5=prctile(zzk1,5);
                z_below_5=zzk1(zzk1<percentile_5);
                if ~isempty(z_below_5) % Check if there are elements below 5th percentile
                        z_avg_below_5_1500(j,i,k) = mean(z_below_5,'omitnan');
                else
                        z_avg_below_5_1500(j,i,k) = NaN; % Handle case where no values are above 95th percentile
                end

	    end

	%    ddzz=(blatex(j,i)-yy1).^2+(blonex(j,i)-(xx1-5)).^2;
        %    [x1,y1]=find(ddzz==min(ddzz,[],'all','omitnan'));
        %    x1=x1(1);y1=y1(1);
        %    z1w5(j,i,:)=zz(x1,y1,:);

	%    ddzz=(blatex(j,i)-yy1).^2+(blonex(j,i)-(xx1-10)).^2;
        %    [x1,y1]=find(ddzz==min(ddzz,[],'all','omitnan'));
        %    x1=x1(1);y1=y1(1);
        %    z1w10(j,i,:)=zz(x1,y1,:);

        %    ddzz=(blatex(j,i)-yy1).^2+(blonex(j,i)-(xx1+5)).^2;
        %    [x1,y1]=find(ddzz==min(ddzz,[],'all','omitnan'));
        %    x1=x1(1);y1=y1(1);
        %    z1e5(j,i,:)=zz(x1,y1,:);

        %    ddzz=(blatex(j,i)-yy1).^2+(blonex(j,i)-(xx1+10)).^2;
        %    [x1,y1]=find(ddzz==min(ddzz,[],'all','omitnan'));
        %    x1=x1(1);y1=y1(1);
        %    z1e10(j,i,:)=zz(x1,y1,:);

        end
    end
    disp(i)
end

%%
disp('prepare temperature')
epath='/glade/campaign/collections/rda/data/d633001/e5.moda.an.pl';
t1=nan(a,b,37);gradt_ns1=t1;gradt_ew1=t1;
for i=1:b
    tname=[epath,'/',num2str(bseason(i)),'/e5.moda.an.pl.128_130_t.ll025sc.',...
        num2str(bseason(i)),'010100_',num2str(bseason(i)),'120100.nc'];
    t=ncread(tname,'T');
    % level=ncread(zname,'level');
    %lat=ncread(tname,'latitude');
    %lon=ncread(tname,'longitude');
    %lon(lon>180)=lon(lon>180)-360;
    %[xx1,yy1]=meshgrid(lon,lat);
    for j=1:a
        if ~isnan(blatex(j,i))
            dd=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'dd'));
            mm=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'mm'));
            if dd<15 && mm~=1
                tt=t(:,:,:,mm-1)+(t(:,:,:,mm)-t(:,:,:,mm-1))*(dd+15)/30;
            elseif dd<15 && mm==1
                tname=[epath,'/',num2str(bseason(i)-1),'/e5.moda.an.pl.128_130_t.ll025sc.',...
                    num2str(bseason(i)-1),'010100_',num2str(bseason(i)-1),'120100.nc'];
                tt2=ncread(tname,'T');
                tt=t(:,:,:,12)+(t(:,:,:,mm)-tt2(:,:,:,12))*(dd+15)/30;
            elseif dd>=15 && mm~=12
                tt=t(:,:,:,mm)+(t(:,:,:,mm+1)-t(:,:,:,mm))*(dd-15)/30;
            else
                tname=[epath,'/',num2str(bseason(i)+1),'/e5.moda.an.pl.128_130_t.ll025sc.',...
                    num2str(bseason(i)+1),'010100_',num2str(bseason(i)+1),'120100.nc'];
                tt2=ncread(tname,'T');
                tt=t(:,:,:,mm)+(tt2(:,:,:,1)-t(:,:,:,mm))*(dd-15)/30;
            end
            %ddtt=(blatex(j,i)-yy1).^2+(blonex(j,i)-xx1).^2;
            %[x1,y1]=find(ddtt==min(ddtt,[],'all','omitnan'));
            %x1=x1(1);y1=y1(1);
            %t1(j,i,:)=tt(x1,y1,:);
	    
	    dx=2*pi*re*cosd(blatex(j,i)).*(blonex(j,i)-xx1)/360;
            dy=2*pi*re*(blatex(j,i)-yy1)/360;
            dxy=sqrt(dx.^2+dy.^2)';

	    gradT_x_km = zeros(size(tt)); % To store gradient in x direction in km
            gradT_y_km = zeros(size(tt)); % To store gradient in y direction in km
            for k = 1:size(tt, 3)
                 % Calculate gradient for the k-th layer
                [gradT_x, gradT_y] = gradient(tt(:,:,k));

                % Divide by spatial distance to get gradient in km
                gradT_x_km(:,:,k) = gradT_x ./ dlat';
                gradT_y_km(:,:,k) = gradT_y ./ dlon;

            end

	    for k=1:size(tt, 3)
                ttk=tt(:,:,k);
                t1(j,i,k)=mean(ttk(dxy<=600),'omitnan');
		gradt_x1k=gradT_x_km(:,:,k);
                gradt_ns1(j,i,k)=mean(gradt_x1k(dxy<=600),'omitnan');
                gradt_y1k=gradT_y_km(:,:,k);
                gradt_ew1(j,i,k)=mean(gradt_y1k(dxy<=600),'omitnan');
            end

        end
    end
    disp(i)
end

disp('prepare SST')
epath='/glade/campaign/collections/rda/data/d633001/e5.moda.an.sfc';
sst1=nan(a,b);
for i=1:b
    sstname=[epath,'/',num2str(bseason(i)),'/e5.moda.an.sfc.128_034_sstk.ll025sc.',...
        num2str(bseason(i)),'010100_',num2str(bseason(i)),'120100.nc'];
    sst=ncread(sstname,'SSTK');
    % level=ncread(zname,'level');
    %lat=ncread(sstname,'latitude');
    %lon=ncread(sstname,'longitude');
    %lon(lon>180)=lon(lon>180)-360;
    %[xx1,yy1]=meshgrid(lon,lat);
    for j=1:a
        if ~isnan(blatex(j,i))
            dd=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'dd'));
            mm=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'mm'));
            if dd<15 && mm~=1
                sstt=sst(:,:,mm-1)+(sst(:,:,mm)-sst(:,:,mm-1))*(dd+15)/30;
            elseif dd<15 && mm==1
                sstname=[epath,'/',num2str(bseason(i)-1),'/e5.moda.an.sfc.128_034_sstk.ll025sc.',...
                    num2str(bseason(i)-1),'010100_',num2str(bseason(i)-1),'120100.nc'];
                sstt2=ncread(sstname,'SSTK');
                sstt=sst(:,:,12)+(sst(:,:,mm)-sstt2(:,:,12))*(dd+15)/30;
            elseif dd>=15 && mm~=12
                sstt=sst(:,:,mm)+(sst(:,:,mm+1)-sst(:,:,mm))*(dd-15)/30;
            else
                sstname=[epath,'/',num2str(bseason(i)+1),'/e5.moda.an.sfc.128_034_sstk.ll025sc.',...
                    num2str(bseason(i)+1),'010100_',num2str(bseason(i)+1),'120100.nc'];
                sstt2=ncread(sstname,'SSTK');
                sstt=sst(:,:,mm)+(sstt2(:,:,1)-sst(:,:,mm))*(dd-15)/30;
            end
            %ddtt=(blatex(j,i)-yy1).^2+(blonex(j,i)-xx1).^2;
            %[x1,y1]=find(ddtt==min(ddtt,[],'all','omitnan'));
            %x1=x1(1);y1=y1(1);
            %sst1(j,i)=sstt(x1,y1);

	    dx=2*pi*re*cosd(blatex(j,i)).*(blonex(j,i)-xx1)/360;
            dy=2*pi*re*(blatex(j,i)-yy1)/360;
            dxy=sqrt(dx.^2+dy.^2)';
            sst1(j,i)=mean(sstt(dxy<=600),'omitnan');
        end
    end
    disp(i)
end

%%
disp('prepare relative humidity')
epath='/glade/campaign/collections/rda/data/d633001/e5.moda.an.pl';
rh1=nan(a,b,37);
for i=1:b
    rhname=[epath,'/',num2str(bseason(i)),'/e5.moda.an.pl.128_157_r.ll025sc.',...
        num2str(bseason(i)),'010100_',num2str(bseason(i)),'120100.nc'];
    rh=ncread(rhname,'R');
    % level=ncread(zname,'level');
    %lat=ncread(rhname,'latitude');
    %lon=ncread(rhname,'longitude');
    %lon(lon>180)=lon(lon>180)-360;
    %[xx1,yy1]=meshgrid(lon,lat);
    for j=1:a
        if ~isnan(blatex(j,i))
            dd=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'dd'));
            mm=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'mm'));
            if dd<15 && mm~=1
                rhm=rh(:,:,:,mm-1)+(rh(:,:,:,mm)-rh(:,:,:,mm-1))*(dd+15)/30;
            elseif dd<15 && mm==1
                rhname=[epath,'/',num2str(bseason(i)-1),'/e5.moda.an.pl.128_157_r.ll025sc.',...
                    num2str(bseason(i)-1),'010100_',num2str(bseason(i)-1),'120100.nc'];
                rh2=ncread(rhname,'R');
                rhm=rh(:,:,:,12)+(rh(:,:,:,mm)-rh2(:,:,:,12))*(dd+15)/30;
            elseif dd>=15 && mm~=12
                rhm=rh(:,:,:,mm)+(rh(:,:,:,mm+1)-rh(:,:,:,mm))*(dd-15)/30;
            else
                rhname=[epath,'/',num2str(bseason(i)+1),'/e5.moda.an.pl.128_157_r.ll025sc.',...
                    num2str(bseason(i)+1),'010100_',num2str(bseason(i)+1),'120100.nc'];
                rh2=ncread(rhname,'R');
                rhm=rh(:,:,:,mm)+(rh2(:,:,:,1)-rh(:,:,:,mm))*(dd-15)/30;
            end
            %ddrh=(blatex(j,i)-yy1).^2+(blonex(j,i)-xx1).^2;
            %[x1,y1]=find(ddrh==min(ddrh,[],'all','omitnan'));
            %x1=x1(1);y1=y1(1);
            %rh1(j,i,:)=rhm(x1,y1,:);
	    dx=2*pi*re*cosd(blatex(j,i)).*(blonex(j,i)-xx1)/360;
            dy=2*pi*re*(blatex(j,i)-yy1)/360;
            dxy=sqrt(dx.^2+dy.^2)';
            for k=1:size(rh, 3)
                rhk=rhm(:,:,k);
                rh1(j,i,k)=mean(rhk(dxy<=600),'omitnan');
            end
        end
    end
    disp(i)
end

%%
disp('prepare u')
epath='/glade/campaign/collections/rda/data/d633001/e5.moda.an.pl';
u1=nan(a,b,37);u1s=u1;
for i=1:b
    uname=[epath,'/',num2str(bseason(i)),'/e5.moda.an.pl.128_131_u.ll025uv.',...
        num2str(bseason(i)),'010100_',num2str(bseason(i)),'120100.nc'];
    u=ncread(uname,'U');
    % level=ncread(zname,'level');
    %lat=ncread(uname,'latitude');
    %lon=ncread(uname,'longitude');
    %lon(lon>180)=lon(lon>180)-360;
    %[xx1,yy1]=meshgrid(lon,lat);
    for j=1:a
        if ~isnan(blatex(j,i))
            dd=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'dd'));
            mm=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'mm'));
            if dd<15 && mm~=1
                uu=u(:,:,:,mm-1)+(u(:,:,:,mm)-u(:,:,:,mm-1))*(dd+15)/30;
            elseif dd<15 && mm==1
                uname=[epath,'/',num2str(bseason(i)-1),'/e5.moda.an.pl.128_131_u.ll025uv.',...
                    num2str(bseason(i)-1),'010100_',num2str(bseason(i)-1),'120100.nc'];
                u2=ncread(uname,'U');
                uu=u(:,:,:,12)+(u(:,:,:,mm)-u2(:,:,:,12))*(dd+15)/30;
            elseif dd>=15 && mm~=12
                uu=u(:,:,:,mm)+(u(:,:,:,mm+1)-u(:,:,:,mm))*(dd-15)/30;
            else
                uname=[epath,'/',num2str(bseason(i)+1),'/e5.moda.an.pl.128_131_u.ll025uv.',...
                    num2str(bseason(i)+1),'010100_',num2str(bseason(i)+1),'120100.nc'];
                u2=ncread(uname,'U');
                uu=u(:,:,:,mm)+(u2(:,:,:,1)-u(:,:,:,mm))*(dd-15)/30;
            end
            %ddu=(blatex(j,i)-yy1).^2+(blonex(j,i)-xx1).^2;
            %[x1,y1]=find(ddu==min(ddu,[],'all','omitnan'));
            %x1=x1(1);y1=y1(1);
            %u1(j,i,:)=uu(x1,y1,:);
	    dx=2*pi*re*cosd(blatex(j,i)).*(blonex(j,i)-xx1)/360;
            dy=2*pi*re*(blatex(j,i)-yy1)/360;
            dxy=sqrt(dx.^2+dy.^2)';
            for k=1:size(uu, 3)
                uuk=uu(:,:,k);
                u1s(j,i,k)=mean(uuk(dxy<=800 & dxy>=300),'omitnan');
		u1(j,i,k)=mean(uuk(dxy<=600),'omitnan');
            end
        end
    end
    disp(i)
end

%%
disp('prepare v')
epath='/glade/campaign/collections/rda/data/d633001/e5.moda.an.pl';
v1=nan(a,b,37);v1s=v1;
for i=1:b
    vname=[epath,'/',num2str(bseason(i)),'/e5.moda.an.pl.128_132_v.ll025uv.',...
        num2str(bseason(i)),'010100_',num2str(bseason(i)),'120100.nc'];
    v=ncread(vname,'V');
    % level=ncread(zname,'level');
    %lat=ncread(vname,'latitude');
    %lon=ncread(vname,'longitude');
    %lon(lon>180)=lon(lon>180)-360;
    %[xx1,yy1]=meshgrid(lon,lat);
    for j=1:a
        if ~isnan(blatex(j,i))
            dd=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'dd'));
            mm=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'mm'));
            if dd<15 && mm~=1
                vv=v(:,:,:,mm-1)+(v(:,:,:,mm)-v(:,:,:,mm-1))*(dd+15)/30;
            elseif dd<15 && mm==1
                vname=[epath,'/',num2str(bseason(i)-1),'/e5.moda.an.pl.128_132_v.ll025uv.',...
                    num2str(bseason(i)-1),'010100_',num2str(bseason(i)-1),'120100.nc'];
                v2=ncread(vname,'V');
                vv=v(:,:,:,12)+(v(:,:,:,mm)-v2(:,:,:,12))*(dd+15)/30;
            elseif dd>=15 && mm~=12
                vv=v(:,:,:,mm)+(v(:,:,:,mm+1)-v(:,:,:,mm))*(dd-15)/30;
            else
                vname=[epath,'/',num2str(bseason(i)+1),'/e5.moda.an.pl.128_132_v.ll025uv.',...
                    num2str(bseason(i)+1),'010100_',num2str(bseason(i)+1),'120100.nc'];
                v2=ncread(vname,'V');
                vv=v(:,:,:,mm)+(v2(:,:,:,1)-v(:,:,:,mm))*(dd-15)/30;
            end
            %ddv=(blatex(j,i)-yy1).^2+(blonex(j,i)-xx1).^2;
            %[x1,y1]=find(ddv==min(ddv,[],'all','omitnan'));
            %x1=x1(1);y1=y1(1);
            %v1(j,i,:)=vv(x1,y1,:);
	    dx=2*pi*re*cosd(blatex(j,i)).*(blonex(j,i)-xx1)/360;
            dy=2*pi*re*(blatex(j,i)-yy1)/360;
            dxy=sqrt(dx.^2+dy.^2)';
            for k=1:size(vv, 3)
                vvk=vv(:,:,k);
                v1s(j,i,k)=mean(vvk(dxy<=800 & dxy>=300),'omitnan');
		v1(j,i,k)=mean(vvk(dxy<=600),'omitnan');
            end
        end
    end
    disp(i)
end

%%
disp('prepare pv')
epath='/glade/campaign/collections/rda/data/d633001/e5.moda.an.pl';
pv1=nan(a,b,37);
for i=1:b
    pvname=[epath,'/',num2str(bseason(i)),'/e5.moda.an.pl.128_060_pv.ll025sc.',...
        num2str(bseason(i)),'010100_',num2str(bseason(i)),'120100.nc'];
    pv=ncread(pvname,'PV');
    % level=ncread(zname,'level');
    %lat=ncread(pvname,'latitude');
    %lon=ncread(pvname,'longitude');
    %lon(lon>180)=lon(lon>180)-360;
    %[xx1,yy1]=meshgrid(lon,lat);
    for j=1:a
        if ~isnan(blatex(j,i))
            dd=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'dd'));
            mm=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'mm'));
            if dd<15 && mm~=1
                pvm=pv(:,:,:,mm-1)+(pv(:,:,:,mm)-pv(:,:,:,mm-1))*(dd+15)/30;
            elseif dd<15 && mm==1
                pvname=[epath,'/',num2str(bseason(i)-1),'/e5.moda.an.pl.128_060_pv.ll025sc.',...
                    num2str(bseason(i)-1),'010100_',num2str(bseason(i)-1),'120100.nc'];
                pv2=ncread(pvname,'PV');
                pvm=pv(:,:,:,12)+(pv(:,:,:,mm)-pv2(:,:,:,12))*(dd+15)/30;
            elseif dd>=15 && mm~=12
                pvm=pv(:,:,:,mm)+(pv(:,:,:,mm+1)-pv(:,:,:,mm))*(dd-15)/30;
            else
                pvname=[epath,'/',num2str(bseason(i)+1),'/e5.moda.an.pl.128_060_pv.ll025sc.',...
                    num2str(bseason(i)+1),'010100_',num2str(bseason(i)+1),'120100.nc'];
                pv2=ncread(pvname,'PV');
                pvm=pv(:,:,:,mm)+(pv2(:,:,:,1)-pv(:,:,:,mm))*(dd-15)/30;
            end
            %ddpv=(blatex(j,i)-yy1).^2+(blonex(j,i)-xx1).^2;
            %[x1,y1]=find(ddpv==min(ddpv,[],'all','omitnan'));
            %x1=x1(1);y1=y1(1);
            %pv1(j,i,:)=pvm(x1,y1,:);
	    dx=2*pi*re*cosd(blatex(j,i)).*(blonex(j,i)-xx1)/360;
            dy=2*pi*re*(blatex(j,i)-yy1)/360;
            dxy=sqrt(dx.^2+dy.^2)';
            for k=1:size(pv, 3)
                pvk=pvm(:,:,k);
                pv1(j,i,k)=mean(pvk(dxy<=600),'omitnan');
            end
        end
    end
    disp(i)
end

%%
disp('prepare vo')
epath='/glade/campaign/collections/rda/data/d633001/e5.moda.an.pl';
vo1=nan(a,b,37);
for i=1:b
    voname=[epath,'/',num2str(bseason(i)),'/e5.moda.an.pl.128_138_vo.ll025sc.',...
        num2str(bseason(i)),'010100_',num2str(bseason(i)),'120100.nc'];
    vo=ncread(voname,'VO');
    % level=ncread(zname,'level');
    %lat=ncread(voname,'latitude');
    %lon=ncread(voname,'longitude');
    %lon(lon>180)=lon(lon>180)-360;
    %[xx1,yy1]=meshgrid(lon,lat);
    for j=1:a
        if ~isnan(blatex(j,i))
            dd=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'dd'));
            mm=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'mm'));
            if dd<15 && mm~=1
                vom=vo(:,:,:,mm-1)+(vo(:,:,:,mm)-vo(:,:,:,mm-1))*(dd+15)/30;
            elseif dd<15 && mm==1
                voname=[epath,'/',num2str(bseason(i)-1),'/e5.moda.an.pl.128_138_vo.ll025sc.',...
                    num2str(bseason(i)-1),'010100_',num2str(bseason(i)-1),'120100.nc'];
                vo2=ncread(voname,'VO');
                vom=vo(:,:,:,12)+(vo(:,:,:,mm)-vo2(:,:,:,12))*(dd+15)/30;
            elseif dd>=15 && mm~=12
                vom=vo(:,:,:,mm)+(vo(:,:,:,mm+1)-vo(:,:,:,mm))*(dd-15)/30;
            else
                voname=[epath,'/',num2str(bseason(i)+1),'/e5.moda.an.pl.128_138_vo.ll025sc.',...
                    num2str(bseason(i)+1),'010100_',num2str(bseason(i)+1),'120100.nc'];
                vo2=ncread(voname,'VO');
                vom=vo(:,:,:,mm)+(vo2(:,:,:,1)-vo(:,:,:,mm))*(dd-15)/30;
            end
            %ddvo=(blatex(j,i)-yy1).^2+(blonex(j,i)-xx1).^2;
            %[x1,y1]=find(ddvo==min(ddvo,[],'all','omitnan'));
            %x1=x1(1);y1=y1(1);
            %vo1(j,i,:)=vom(x1,y1,:);
	    dx=2*pi*re*cosd(blatex(j,i)).*(blonex(j,i)-xx1)/360;
            dy=2*pi*re*(blatex(j,i)-yy1)/360;
            dxy=sqrt(dx.^2+dy.^2)';
            for k=1:size(pv, 3)
                vok=vom(:,:,k);
                vo1(j,i,k)=mean(vok(dxy<=600),'omitnan');
            end
        end
    end
    disp(i)
end

%%
disp('prepare div')
epath='/glade/campaign/collections/rda/data/d633001/e5.moda.an.pl';
div1=nan(a,b,37);
for i=1:b
    divname=[epath,'/',num2str(bseason(i)),'/e5.moda.an.pl.128_155_d.ll025sc.',...
        num2str(bseason(i)),'010100_',num2str(bseason(i)),'120100.nc'];
    div=ncread(divname,'D');
    % level=ncread(zname,'level');
    %lat=ncread(divname,'latitude');
    %lon=ncread(divname,'longitude');
    %lon(lon>180)=lon(lon>180)-360;
    %[xx1,yy1]=meshgrid(lon,lat);
    for j=1:a
        if ~isnan(blatex(j,i))
            dd=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'dd'));
            mm=str2double(datestr(datenum('1858-11-17')+timeex(j,i),'mm'));
            if dd<15 && mm~=1
                divm=div(:,:,:,mm-1)+(div(:,:,:,mm)-div(:,:,:,mm-1))*(dd+15)/30;
            elseif dd<15 && mm==1
                divname=[epath,'/',num2str(bseason(i)-1),'/e5.moda.an.pl.128_155_d.ll025sc.',...
                    num2str(bseason(i)-1),'010100_',num2str(bseason(i)-1),'120100.nc'];
                div2=ncread(divname,'D');
                divm=div(:,:,:,12)+(div(:,:,:,mm)-div2(:,:,:,12))*(dd+15)/30;
            elseif dd>=15 && mm~=12
                divm=div(:,:,:,mm)+(div(:,:,:,mm+1)-div(:,:,:,mm))*(dd-15)/30;
            else
                divname=[epath,'/',num2str(bseason(i)+1),'/e5.moda.an.pl.128_155_d.ll025sc.',...
                    num2str(bseason(i)+1),'010100_',num2str(bseason(i)+1),'120100.nc'];
                div2=ncread(divname,'D');
                divm=div(:,:,:,mm)+(div2(:,:,:,1)-div(:,:,:,mm))*(dd-15)/30;
            end
            %dddiv=(blatex(j,i)-yy1).^2+(blonex(j,i)-xx1).^2;
            %[x1,y1]=find(dddiv==min(dddiv,[],'all','omitnan'));
            %x1=x1(1);y1=y1(1);
            %div1(j,i,:)=divm(x1,y1,:);
	    dx=2*pi*re*cosd(blatex(j,i)).*(blonex(j,i)-xx1)/360;
            dy=2*pi*re*(blatex(j,i)-yy1)/360;
            dxy=sqrt(dx.^2+dy.^2)';
            for k=1:size(div, 3)
                divk=divm(:,:,k);
                div1(j,i,k)=mean(divk(dxy<=600),'omitnan');
            end
        end
    end
    disp(i)
end

%%
save('preprocess_AR_1981-2019_600km_v1.mat','blonex','blatex','blandfall','pi1','u1','v1','u1s','v1s','pv1','vo1','div1','rh1','t1','z1','gradz_ew1','gradz_ns1','gradt_ew1','gradt_ns1','zzmax1000','zzmin1000','zzmax_x1000','zzmin_x1000','zzmax1500','zzmin1500','zzmax_x1500','zzmin_x1500','sst1','z_avg_above_95_1000','z_avg_below_5_1000','z_avg_above_95_1500','z_avg_below_5_1500')









