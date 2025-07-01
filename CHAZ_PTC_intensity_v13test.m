clear;clc

fpath='/Users/dianyili/Desktop/SBU/CHAZ/';
fpathout='/Users/dianyili/Desktop/SBU/CHAZ/atl_era5/';

load 'AR_trainingset_v11.mat'
Xm=Xm_best;Xs=Xs_best;Ym=Xm(2);Ys=Xs(2);

level = [1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225,...
    250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825,...
    850, 875, 900, 925, 950, 975, 1000 ];
pen=nan(125,700,40,40);

dp_fct=nan(125,700,40,40);

rng(42)

load('gev_params_interp_by_lat.mat');
lat_list = interp_params(:,1);
k_list = interp_params(:,2);
sigma_list = interp_params(:,3);
mu_list = interp_params(:,4);

for ens=1:40
    disp(['Run ens ',num2str(ens)])
    if ens-1<10
        fname=['atl_2019_2ens00',num2str(ens-1),'_ETv1.nc'];
    else
        fname=['atl_2019_2ens0',num2str(ens-1),'_ETv1.nc'];
    end
    ff=[fpathout,fname];
    lon=ncread(ff,'longitude')';
    lat=ncread(ff,'latitude')';
    ex=ncread(ff,'ex');
    p0=ncread(ff,'p0');
    
    dx=nan(125,700);dy=nan(125,700);
    for j=1+1:125-1
        dx(j,1:size(lon,2))=(lon(j+1,:)-lon(j-1,:))/360*2*pi*6371.*cosd(lat(j,:))/12;
        dy(j,1:size(lon,2))=(lat(j+1,:)-lat(j-1,:))/360*2*pi*6371/12;
    end
    sp=sqrt(dx.^2+dy.^2);

    load([fpath,'CHAZ_PTC_interpvars1_ens',num2str(ens),'_1981_2019.mat'])
    load([fpath,'CHAZ_PTC_interpvars2_ens',num2str(ens),'_1981_2019.mat'])

    pv850=pv1(:,:,31);
    u700=u1(:,:,26);v700=v1(:,:,26);
    u850=u1(:,:,31);v850=v1(:,:,31);
    uv850=sqrt(u850.^2+v850.^2);
    shr700850=sqrt((u700-u850).^2+(v700-v850).^2);
    
    theta=nan(125,700);
    for i=1:37
        theta(:,:,i)=t1(:,:,i).*(1000/level(i))^(0.286);
    end
    ss=nan(size(theta));
    for i=2:36
        ss(:,:,i)=(theta(:,:,i+1)-theta(:,:,i-1))./(z1(:,:,i+1)-z1(:,:,i-1));
    end
    ss700=ss(:,:,26);
    rh850=rh1(:,:,31);

    duv=nan(125,700);dshr=duv;dss=duv;drh=duv;
    for i=1+2:125
        duv(i,:)=uv850(i,:)-uv850(i-2,:);
        dshr(i,:)=shr700850(i,:)-shr700850(i-2,:);
        dss(i,:)=ss700(i,:)-ss700(i-2,:);
        drh(i,:)=rh850(i,:)-rh850(i-2,:);
    end

    for ensi=1:40
        dpp12=nan(125,700);
        p=nan(125,700);
        p12=nan(125,700);
        idx=nan(1,700);
        
        exensi=ex(:,:,ensi)';
        pensi=p0(:,:,ensi)';
        for j=1:size(exensi,2)
            x=find(exensi(:,j)==1);
            if ~isempty(x)
                if x(1)-2>0 && ~isnan(pensi(x(1)-2,j))
                    p(x(1)-2:end,j)=pensi(x(1)-2:end,j);
                    idx(j)=x(1)-2;
                else
                    p(x(1):end,j)=pensi(x(1):end,j);
                end
            end
        end
        for j=1+2:125
            dpp12(j,:)=p(j,:)-p(j-2,:);
            p12(j,:)=p(j-2,:);
        end

        rd_dis=rand(size(p));

        prep=nan(size(p));
        for i=1:700
            ex_num=find(~isnan(p(:,i)));
            if ~isempty(ex_num)
                if length(ex_num)>=3
                    p(ex_num(1)+1:ex_num(3)-1,i)=nan;
                    prep(ex_num(1:3),i)=p(ex_num(1:3),i);
                    for j=ex_num(3):2:125-2
                        Xp=[prep(j-2,i),prep(j,i)-prep(j-2,i),sp(j,i)];%3
                        % Xp=[prep(j-2,i),prep(j,i)-prep(j-2,i),sp(j,i),pv850(j,i),duv(j,i)];%5
                        % Xp=[prep(j-2,i),prep(j,i)-prep(j-2,i),sp(j,i),pv850(j,i),duv(j,i),dshr(j,i),dss(j,i),drh(j,i)];%8
                        XXp=(Xp-Xm)./Xs;
                        % YYp=(predict(mdl_best,XXp));
                        YYp=(predict(mdl_weighted_best,XXp));
                        r = random(pd_t_best);
                        q_low = icdf(pd_t_best, 0.025);
                        q_high = icdf(pd_t_best, 0.975);
                        r = min(max(r, q_low), q_high);
                        YYp = YYp + res_coef * r;

                        YYp(YYp>5*Ys+Ym)=5*Ys+Ym;%cut at 5 sigma
                        YYp(YYp<-5*Ys+Ym)=-5*Ys+Ym;
                        prep(j+2,i)=YYp+prep(j,i);
                        % dissipation_prob =cdf(pd_gev, prep(j+2,i));
                        latpre=round(lat(j+2,i));
                        latidx=find(abs(latpre-lat_list)<.001);
                        dissipation_prob =gevcdf(prep(j+2,i),k_list(latidx),sigma_list(latidx),mu_list(latidx));

                        if isnan(latpre); break;end
                        if rd_dis(j,i)^(2/5)<dissipation_prob; break;end
                    end
                else
                    prep(:,i)=p(:,i);
                end
            end
            if ~isnan(idx(i)); prep(idx(i),i)=nan;end
        end

        pen(:,:,ensi,ens)=prep;
    end
end

pen1=pen;


%%
for i=1:40
   
    disp(['write intensity to ens ',num2str(i)])


    if i-1<10
        fname_et=['atl_2019_2ens00',num2str(i-1),'_ETv1.nc'];
    else
        fname_et=['atl_2019_2ens0',num2str(i-1),'_ETv1.nc'];
    end
    ff_et=[fpathout,fname_et];
    clon=ncread(ff_et,'longitude');
    clat=ncread(ff_et,'latitude');
    % cex=ncread(ff_et,'ex');
    
    cpres=nan([size(clat),40]);
    for ensi=1:40
        cpres(:,:,ensi)=squeeze(pen1(:,1:length(clat(:,1)),ensi,i))';
    end
    
    cpres_r=nan([size(clat),40]);
    for ensi=1:40
        data=cpres(:,:,ensi);
        [n, nTime] = size(data);
        for Idx = 1:n
            for t = 1:nTime
                if isnan(data(Idx, t))
                    left = find(~isnan(data(Idx, 1:t-1)), 1, 'last');
                    right = find(~isnan(data(Idx, t+1:end)), 1, 'first');
                    if ~isempty(left) && ~isempty(right)
                        leftValue = data(Idx, left);
                        rightValue = data(Idx, t + right);
                        data(Idx, t) = mean([leftValue, rightValue]);
                    else
                        data(Idx, t) = NaN;
                    end
                end
            end
        end
        cpres_r(:,:,ensi)=data;
    end


    [a,b,c]=size(cpres_r);

    if i-1<10
        fname_ptc=['atl_2019_2ens00',num2str(i-1),'_PTCv1_3pred.nc'];
    else
        fname_ptc=['atl_2019_2ens0',num2str(i-1),'_PTCv1_3pred.nc'];
    end
    ff_ptc=[fpathout,fname_ptc];

    copyfile(ff_et,ff_ptc)

    nccreate(ff_ptc,'p0PTC','Dimensions',{"stormID",a,"lifelength",b,"intensitymember",c})
    ncwrite(ff_ptc,'p0PTC',cpres_r);
    ncwriteatt(ff_ptc,'p0PTC','long_name','mslp of post-tropical cyclone')

end




