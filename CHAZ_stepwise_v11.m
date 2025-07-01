%training for PTC intensity
clear;clc

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

blatex1=nan(size(blat));blonex1=blatex1;bwindex=blatex1;presex=blatex1;
byrex=nan(1,length(blat(1,:)));blandfallex=blatex1;
for i=1:length(blat(1,:))
    for j=1:360
        if (bstatus(1,j,i)=='E' && bstatus(2,j,i)=='T') 
            if blandfall(j,i)<111*0.25/2;break;end

            presex(j:end,i)=bpres(j:end,i);
            bwindex(j:end,i)=bwind(j:end,i);
            blonex1(j:end,i)=blon(j:end,i);
            blatex1(j:end,i)=blat(j:end,i); 
            blandfallex(j:end,i)=blandfall(j:end,i); 
            break
        end
    end
end


[a,b]=size(blonex1);
%%

level = [1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225,... 
     250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825,... 
     850, 875, 900, 925, 950, 975, 1000 ];
%cmip6 daily 1000., 850., 700., 500., 250., 100., 50., 10.

load('/Users/dianyili/Desktop/SBU/CHAZ/preprocess_AR_1981-2019_600km_v1.mat')
load('preprocess_EGR_1981-2019_600km_v1.mat')
presex(:,size(presex,2)+1:size(sst1,2))=nan;

uv=sqrt(u1.^2+v1.^2);
uvdir = atan2(u1, v1) * (180/pi);  % Convert radians to degrees
uvdir = mod(uvdir + 360, 360);  % Ensure angle is between 0 and 360
for i=1:37
    theta(:,:,i)=t1(:,:,i).*(1000/level(i))^(0.286);
end
ss=nan(size(theta));
for i=2:36
    ss(:,:,i)=(theta(:,:,i+1)-theta(:,:,i-1))./(z1(:,:,i+1)-z1(:,:,i-1));
end

dp12=nan(size(blonex));dpp12=dp12;presex12=dp12;dx=dp12;dy=dp12;dlon=dp12;dlat=dp12;dpi=dp12;dpi2=dp12;dsp=dp12;dsst=dp12;
drh=nan(size(rh1));duv=drh;dz=drh;dss=drh;dss1=drh;duvdir=drh;dt=drh;dpv=drh;dvo=dpv;
du=dpv;dv=dpv;dgradt=drh;dgradz=drh;ddiv=drh;
degr=nan(size(egr1));
for j=1:a-4
    dp12(j,:)=presex(j+4,:)-presex(j,:);
end
for j=1+2:a-2
    dx(j,:)=(blonex(j+2,:)-blonex(j-2,:))/360*2*pi*6371.*cosd(blatex(j,:))/12;
    dy(j,:)=(blatex(j+2,:)-blatex(j-2,:))/360*2*pi*6371/12;
end
sp=sqrt(dx.^2+dy.^2);

gradt=sqrt(gradt_ns1.^2+gradt_ew1.^2);
gradz=sqrt(gradz_ns1.^2+gradz_ew1.^2);


for j=1+4:a
    dpp12(j,:)=presex(j,:)-presex(j-4,:);
    dlon(j,:)=blonex(j,:)-blonex(j-4,:);
    dlat(j,:)=blatex(j,:)-blatex(j-4,:);
    dpi(j,:)=pi1(j,:)-pi1(j-4,:);
    dpi2(j,:)=pi1(j,:).^2-pi1(j-4,:).^2;
    dsst(j,:)=sst1(j,:)-sst1(j-4,:);
    dsp(j,:)=sp(j,:)-sp(j-4,:);
    presex12(j,:)=presex(j-4,:);
    drh(j,:,:)=rh1(j,:,:)-rh1(j-4,:,:);
    duv(j,:,:)=uv(j,:,:)-uv(j-4,:,:);
    du(j,:,:)=u1(j,:,:)-u1(j-4,:,:);
    dv(j,:,:)=v1(j,:,:)-v1(j-4,:,:);
    dz(j,:,:)=z1(j,:,:)-z1(j-4,:,:);
    dss(j,:,:)=ss(j,:,:)-ss(j-4,:,:);
    dss1(j,:,:)=1./ss(j,:,:)-1./ss(j-4,:,:);
    duvdir(j,:,:)=uvdir(j,:,:)-uvdir(j-4,:,:);
    dt(j,:,:)=t1(j,:,:)-t1(j-4,:,:);
    dpv(j,:,:)=pv1(j,:,:)-pv1(j-4,:,:);
    dvo(j,:,:)=vo1(j,:,:)-vo1(j-4,:,:);
    degr(j,:,:)=egr1(j,:,:)-egr1(j-4,:,:);
    dgradt(j,:,:)=gradt(j,:,:)-gradt(j-4,:,:);
    dgradz(j,:,:)=gradz(j,:,:)-gradz(j-4,:,:);
    ddiv(j,:,:)=div1(j,:,:)-div1(j-4,:,:);
end

%%
%rh
rh100=rh1(:,:,11);
rh250=rh1(:,:,17);
rh500=rh1(:,:,22);
rh700=rh1(:,:,26);
rh850=rh1(:,:,31);
drh100=drh(:,:,11);
drh250=drh(:,:,17);
drh500=drh(:,:,22);
drh700=drh(:,:,26);
drh850=drh(:,:,31);
%uv

uv100=uv(:,:,11);
uv250=uv(:,:,17);
uv500=uv(:,:,22);
uv700=uv(:,:,26);
uv850=uv(:,:,31);
duv100=duv(:,:,11);
duv250=duv(:,:,17);
duv500=duv(:,:,22);
duv700=duv(:,:,26);
duv850=duv(:,:,31);

du100=du(:,:,11);
du250=du(:,:,17);
du500=du(:,:,22);
du700=du(:,:,26);
du850=du(:,:,31);

dv100=dv(:,:,11);
dv250=dv(:,:,17);
dv500=dv(:,:,22);
dv700=dv(:,:,26);
dv850=dv(:,:,31);

uvdir100=uvdir(:,:,11);
uvdir250=uvdir(:,:,17);
uvdir500=uvdir(:,:,22);
uvdir700=uvdir(:,:,26);
uvdir850=uvdir(:,:,31);
 
duvdir100=duvdir(:,:,11);
duvdir250=duvdir(:,:,17);
duvdir500=duvdir(:,:,22);
duvdir700=duvdir(:,:,26);
duvdir850=duvdir(:,:,31);

%shear

u100=u1s(:,:,11);
u250=u1s(:,:,17);
u500=u1s(:,:,22);
u700=u1s(:,:,26);
u850=u1s(:,:,31);
v100=v1s(:,:,11);
v250=v1s(:,:,17);
v500=v1s(:,:,22);
v700=v1s(:,:,26);
v850=v1s(:,:,31);

shr100250=sqrt((u100-u250).^2+(v100-v250).^2);
shr100500=sqrt((u100-u500).^2+(v100-v500).^2);
shr100700=sqrt((u100-u700).^2+(v100-v700).^2);
shr100850=sqrt((u100-u850).^2+(v100-v850).^2);
shr250850=sqrt((u250-u850).^2+(v250-v850).^2);
shr250700=sqrt((u250-u700).^2+(v250-v700).^2);
shr250500=sqrt((u250-u500).^2+(v250-v500).^2);
shr500850=sqrt((u500-u850).^2+(v500-v850).^2);
shr500700=sqrt((u500-u700).^2+(v500-v700).^2);
shr700850=sqrt((u700-u850).^2+(v700-v850).^2);


dshr250850=nan(size(shr250850));
dshr250700=dshr250850;dshr250500=dshr250850;dshr500850=dshr250850;dshr500700=dshr250850;dshr700850=dshr250850;
dshr100250=dshr250850;dshr100500=dshr250850;dshr100700=dshr250850;dshr100850=dshr250850;
for j=1+4:a
    dshr100250(j,:)=shr100250(j,:)-shr100250(j-4,:);
    dshr100500(j,:)=shr100500(j,:)-shr100500(j-4,:);
    dshr100700(j,:)=shr100700(j,:)-shr100700(j-4,:);
    dshr100850(j,:)=shr100850(j,:)-shr100850(j-4,:);
    dshr250850(j,:)=shr250850(j,:)-shr250850(j-4,:);
    dshr250700(j,:)=shr250700(j,:)-shr250700(j-4,:);
    dshr250500(j,:)=shr250500(j,:)-shr250500(j-4,:);
    dshr500850(j,:)=shr500850(j,:)-shr500850(j-4,:);
    dshr500700(j,:)=shr500700(j,:)-shr500700(j-4,:);
    dshr700850(j,:)=shr700850(j,:)-shr700850(j-4,:);  
end

% geopotential height
z100=z1(:,:,11);
z250=z1(:,:,17);
z500=z1(:,:,22);
z700=z1(:,:,26);
z850=z1(:,:,31);
dz100=dz(:,:,11);
dz250=dz(:,:,17);
dz500=dz(:,:,22);
dz700=dz(:,:,26);
dz850=dz(:,:,31);

%statical stability
ss100=ss(:,:,11);
ss250=ss(:,:,17);
ss500=ss(:,:,22);
ss700=ss(:,:,26);
ss850=ss(:,:,31);

dss100=dss(:,:,11);
dss250=dss(:,:,17);
dss500=dss(:,:,22);
dss700=dss(:,:,26);
dss850=dss(:,:,31);

%t
t100=t1(:,:,11);
t250=t1(:,:,17);
t500=t1(:,:,22);
t700=t1(:,:,26);
t850=t1(:,:,31);

dt100=dt(:,:,11);
dt250=dt(:,:,17);
dt500=dt(:,:,22);
dt700=dt(:,:,26);
dt850=dt(:,:,31);

%pv1
pv100=pv1(:,:,11);
pv250=pv1(:,:,17);
pv500=pv1(:,:,22);
pv700=pv1(:,:,26);
pv850=pv1(:,:,31);

dpv100=dpv(:,:,11);
dpv250=dpv(:,:,17);
dpv500=dpv(:,:,22);
dpv700=dpv(:,:,26);
dpv850=dpv(:,:,31);

%vo1
vo100=vo1(:,:,11);
vo250=vo1(:,:,17);
vo500=vo1(:,:,22);
vo700=vo1(:,:,26);
vo850=vo1(:,:,31);

dvo100=dvo(:,:,11);
dvo250=dvo(:,:,17);
dvo500=dvo(:,:,22);
dvo700=dvo(:,:,26);
dvo850=dvo(:,:,31);
%div1
div100=div1(:,:,11);
div250=div1(:,:,17);
div500=div1(:,:,22);
div700=div1(:,:,26);
div850=div1(:,:,31);

ddiv100=ddiv(:,:,11);
ddiv250=ddiv(:,:,17);
ddiv500=ddiv(:,:,22);
ddiv700=ddiv(:,:,26);
ddiv850=ddiv(:,:,31);

%zzmax
% zzmax=zzmax./z1;
zzmax=zzmax1000;
zzmin=zzmin1000;
zzmax_x=zzmax_x1000;
zzmin_x=zzmin_x1000;

zzmax100=zzmax(:,:,11);
zzmax250=zzmax(:,:,17);
zzmax500=zzmax(:,:,22);
zzmax700=zzmax(:,:,26);
zzmax850=zzmax(:,:,31);

% zzmin=zzmin./z1;
zzmin100=zzmin(:,:,11);
zzmin250=zzmin(:,:,17);
zzmin500=zzmin(:,:,22);
zzmin700=zzmin(:,:,26);
zzmin850=zzmin(:,:,31);

%egr1
egr100=egr1(:,:,1);
egr250=egr1(:,:,2);
egr500=egr1(:,:,3);
egr700=egr1(:,:,4);
egr850=egr1(:,:,5);

degr100=degr(:,:,1);
degr250=degr(:,:,2);
degr500=degr(:,:,3);
degr700=degr(:,:,4);
degr850=degr(:,:,5);

%gradient

gradt100=gradt(:,:,11);
gradt250=gradt(:,:,17);
gradt500=gradt(:,:,22);
gradt700=gradt(:,:,26);
gradt850=gradt(:,:,31);

dgradt100=dgradt(:,:,11);
dgradt250=dgradt(:,:,17);
dgradt500=dgradt(:,:,22);
dgradt700=dgradt(:,:,26);
dgradt850=dgradt(:,:,31);

gradz100=gradz(:,:,11);
gradz250=gradz(:,:,17);
gradz500=gradz(:,:,22);
gradz700=gradz(:,:,26);
gradz850=gradz(:,:,31);

dgradz100=dgradz(:,:,11);
dgradz250=dgradz(:,:,17);
dgradz500=dgradz(:,:,22);
dgradz700=dgradz(:,:,26);
dgradz850=dgradz(:,:,31);
%%
%lev=[11,17,22,26,31];
% 初始化与设置
rng(42);  % 保证随机可复现

% 只要某一列在这几个变量中有一个非 NaN 就算有效
valid_cyclone_mask = any(~isnan(dp12) & ~isnan(dpi) & ~isnan(presex12) & ~isnan(sp), 1); % 逻辑是每列有“至少一个”时为 true
valid_cyclones = find(valid_cyclone_mask);

% 划分训练和测试集
nValidCyclones = numel(valid_cyclones);
nTrain = round(0.8 * nValidCyclones);
shuffled = valid_cyclones(randperm(nValidCyclones));
train_idx = shuffled(1:nTrain);
test_idx = shuffled(nTrain+1:end);

nTrainEachIter = round(0.7 * nTrain);
nIter = 100;

% 所有变量名和变量矩阵
varNames = {'rh100','rh250','rh500','rh700','rh850',...
    'drh100','drh250','drh500','drh700','drh850',...
    'uv100','uv250','uv500','uv700','uv850',...
    'duv100','duv250','duv500','duv700','duv850',...
    'u100','u250','u500','u700','u850',...
    'v100','v250','v500','v700','v850',...
    'du100','du250','du500','du700','du850',...
    'dv100','dv250','dv500','dv700','dv850',...
    'shr100250','shr100500','shr100700','shr100850','shr250500',...
    'shr250700','shr250850','shr500700','shr500850','shr700850',...
    'dshr100250','dshr100500','dshr100700','dshr100850','dshr250500',...
    'dshr250700','dshr250850','dshr500700','dshr500850','dshr700850',...
    'ss100','ss250','ss500','ss700','ss850',...
    'dss100','dss250','dss500','dss700','dss850',...
    't100','t250','t500','t700','t850',...
    'dt100','dt250','dt500','dt700','dt850',...
    'pv100','pv250','pv500','pv700','pv850',...
    'dpv100','dpv250','dpv500','dpv700','dpv850',...
    'vo100','vo250','vo500','vo700','vo850',...
    'dvo100','dvo250','dvo500','dvo700','dvo850',...
    'div100','div250','div500','div700','div850',...
    'ddiv100','ddiv250','ddiv500','ddiv700','ddiv850',...
    'egr100','egr250','egr500','egr700',...
    'degr100','degr250','degr500','degr700',...
    'gradt100','gradt250','gradt500','gradt700',...
    'dgradt100','dgradt250','dgradt500','dgradt700',...
    'gradz100','gradz250','gradz500','gradz700',...
    'dgradz100','dgradz250','dgradz500','dgradz700',...
    'zzmax100','zzmax250','zzmax500','zzmax700','zzmax850',...
    'zzmin100','zzmin250','zzmin500','zzmin700','zzmin850',...
    'presex12','dpp12','blatex','blonex','sp','pi1','sst1','dsst','dlat','dlon','dpi'};

X_vars = {...
    rh100,rh250,rh500,rh700,rh850,...
    drh100,drh250,drh500,drh700,drh850,...
    uv100,uv250,uv500,uv700,uv850,...
    duv100,duv250,duv500,duv700,duv850,...
    u100,u250,u500,u700,u850,...
    v100,v250,v500,v700,v850,...
    du100,du250,du500,du700,du850,...
    dv100,dv250,dv500,dv700,dv850,...
    shr100250,shr100500,shr100700,shr100850,shr250500,...
    shr250700,shr250850,shr500700,shr500850,shr700850,...
    dshr100250,dshr100500,dshr100700,dshr100850,dshr250500,...
    dshr250700,dshr250850,dshr500700,dshr500850,dshr700850,...
    ss100,ss250,ss500,ss700,ss850,...
    dss100,dss250,dss500,dss700,dss850,...
    t100,t250,t500,t700,t850,...
    dt100,dt250,dt500,dt700,dt850,...
    pv100,pv250,pv500,pv700,pv850,...
    dpv100,dpv250,dpv500,dpv700,dpv850,...
    vo100,vo250,vo500,vo700,vo850,...
    dvo100,dvo250,dvo500,dvo700,dvo850,...
    div100,div250,div500,div700,div850,...
    ddiv100,ddiv250,ddiv500,ddiv700,ddiv850,...
    egr100,egr250,egr500,egr700,...
    degr100,degr250,degr500,degr700,...
    gradt100,gradt250,gradt500,gradt700,...
    dgradt100,dgradt250,dgradt500,dgradt700,...
    gradz100,gradz250,gradz500,gradz700,...
    dgradz100,dgradz250,dgradz500,dgradz700,...
    zzmax100./z100,zzmax250./z250,zzmax500./z500,zzmax700./z700,zzmax850./z850,...
    zzmin100./z100,zzmin250./z250,zzmin500./z500,zzmin700./z700,zzmin850./z850,...
    presex12,dpp12,blatex,blonex,sp,pi1,sst1,dsst,dlat,dlon,dpi};

% Stepwise regression
select_counts = zeros(1, numel(varNames));
for i = 1:nIter
    fprintf('Running iteration %d/%d...\n', i, nIter);
    sub_idx = randsample(train_idx, nTrainEachIter);
    
    % 将目标变量和输入变量 reshape 成列向量并去除 NaN
    dp12_vec = reshape(dp12(:,sub_idx), [], 1);
    dpi_vec = reshape(dpi(:,sub_idx), [], 1);
    presex12_vec = reshape(presex12(:,sub_idx), [], 1);
    sp_vec = reshape(sp(:,sub_idx), [], 1);
    valid_idx = find(~isnan(dp12_vec) & ~isnan(dpi_vec) & ~isnan(presex12_vec) & ~isnan(sp_vec));
    if numel(valid_idx) < 20, continue; end

    Ysub = dp12_vec(valid_idx);
    Xsub = zeros(numel(valid_idx), numel(X_vars));
    for j = 1:numel(X_vars)
        Xj = reshape(X_vars{j}(:,sub_idx), [], 1);
        Xsub(:,j) = Xj(valid_idx);
    end

    Xm = mean(Xsub); Xs = std(Xsub);
    XX = (Xsub - Xm)./Xs;
    YY = (Ysub - mean(Ysub))./std(Ysub);

    try
        w = 1 + abs(YY - mean(YY));
        mdl = stepwiselm(XX, YY, 'dp12 ~ presex12 + dpp12', ...
            'Upper', 'linear', 'Criterion', 'BIC', ...
            'Weights', w, ...
            'PEnter', 0, 'PRemove', 0.01, ...
            'VarNames', [varNames, {'dp12'}]);

        selected = mdl.PredictorNames;
        for j = 1:length(selected)
            idx = find(strcmp(varNames, selected{j}));
            if ~isempty(idx)
                select_counts(idx) = select_counts(idx) + 1;
            end
        end
    catch
        continue
    end
end

%%
% === 使用 stepwise 排序后的变量顺序，不做筛选 ===
T = table(varNames(:), select_counts(:), 'VariableNames', {'Variable', 'SelectionCount'});
T = sortrows(T, 'SelectionCount', 'descend');
selected_vars = T.Variable;  % 不筛选，保留全部排序结果

% === 构造训练集 Y （固定）===
dp12_train = reshape(dp12(:,train_idx), [], 1);
dpi_train = reshape(dpi(:,train_idx), [], 1);
presex12_train = reshape(presex12(:,train_idx), [], 1);
sp_train = reshape(sp(:,train_idx), [], 1);
valid_idx = find(~isnan(dp12_train) & ~isnan(dpi_train) & ~isnan(presex12_train) & ~isnan(sp_train));
Y_train = dp12_train(valid_idx);

% === 构造测试集 Y （固定）===
dp12_test = reshape(dp12(:,test_idx), [], 1);
dpi_test = reshape(dpi(:,test_idx), [], 1);
presex12_test = reshape(presex12(:,test_idx), [], 1);
sp_test = reshape(sp(:,test_idx), [], 1);
valid_idx_test = find(~isnan(dp12_test) & ~isnan(dpi_test) & ~isnan(presex12_test) & ~isnan(sp_test));
Y_test = dp12_test(valid_idx_test);

% === 初始化 ===
max_var = min(20, length(selected_vars));
results_lin = [];        % 原本已有
results_resid = [];      % 添加这一行！
results_weighted= [];
cv = cvpartition(length(Y_train), 'KFold', 10);
rmse_cv = nan(max_var,1);
for n_var = 1:max_var
    vars_now = selected_vars(1:n_var);

    % === 构造训练数据 ===
    X_train = zeros(numel(valid_idx), n_var);
    for i = 1:n_var
        idx = find(strcmp(varNames, vars_now{i}));
        Xi = reshape(X_vars{idx}(:,train_idx), [], 1);
        X_train(:,i) = Xi(valid_idx);
    end
    Xm = mean(X_train); Xs = std(X_train); Xs(Xs==0) = 1;
    X_train_std = (X_train - Xm) ./ Xs;

    % === 构造测试数据 ===
    X_test = zeros(numel(valid_idx_test), n_var);
    for i = 1:n_var
        idx = find(strcmp(varNames, vars_now{i}));
        Xi = reshape(X_vars{idx}(:,test_idx), [], 1);
        X_test(:,i) = Xi(valid_idx_test);
    end
    X_test_std = (X_test - Xm) ./ Xs;

    % === 权重定义（尾部更高） ===
    Y_train_centered = Y_train - mean(Y_train);
    w_train = 1 + abs(Y_train_centered);  % 可改为 ^2 强调尾部

    % === 线性回归模型 ===
    mdl_lin = fitlm(X_train_std, Y_train, 'VarNames', [vars_now', {'dp12'}]);
    Y_pred_lin = predict(mdl_lin, X_test_std);

    % === 加权线性回归模型 ===
    mdl_weighted = fitlm(X_train_std, Y_train, 'Weights', w_train, 'VarNames', [vars_now', {'dp12'}]);
    Y_pred_weighted = predict(mdl_weighted, X_test_std);
    Y_pred_weighted_train = predict(mdl_weighted, X_train_std);
    rmse_train = sqrt(mean((Y_train - Y_pred_weighted_train).^2));
    rmse_train_all(n_var, 1) = rmse_train;

    rmse_fold = nan(cv.NumTestSets, 1);
    for k = 1:cv.NumTestSets
        idx_train = training(cv, k);
        idx_test = test(cv, k);

        mdl = fitlm(X_train_std(idx_train,:), Y_train(idx_train), ...
            'Weights', w_train(idx_train), 'VarNames', [vars_now', {'dp12'}]);
        Yhat = predict(mdl, X_test_std);
        rmse_fold(k) = sqrt(mean((Y_test - Yhat).^2));
    end

    rmse_cv(n_var) = mean(rmse_fold);

    % === 残差建模：整体拟合 t 分布 ===
    res_train = Y_train - predict(mdl_weighted, X_train_std);
    pd_t = fitdist(res_train, 'tLocationScale');
    res_sample = random(pd_t, length(Y_pred_weighted), 1);

    q_low = icdf(pd_t, 0.025);
    q_high = icdf(pd_t, 0.975);
    res_sample = min(max(res_sample, q_low), q_high);

    Y_pred_res = Y_pred_weighted + res_sample*.1;

    % === 性能评估指标 ===
    calc_stats = @(y, yhat) struct(...
        'n_var', n_var, ...
        'R2', 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2), ...
        'RMSE', sqrt(mean((y - yhat).^2)), ...
        'MAE', mean(abs(y - yhat)), ...
        'SignAcc', mean(sign(y) == sign(yhat)));

    stats_lin = calc_stats(Y_test, Y_pred_lin);
    stats_resid = calc_stats(Y_test, Y_pred_res);
    stats_weighted = calc_stats(Y_test, Y_pred_weighted);

    results_lin = [results_lin; stats_lin];
    results_resid = [results_resid; stats_resid];
    results_weighted = [results_weighted; stats_weighted];

    % === K-S 检验统计 ===
    [~, p_lin, ks_lin] = kstest2(Y_test, Y_pred_lin);
    [~, p_res, ks_res] = kstest2(Y_test, Y_pred_res);
    [~, p_wgt, ks_wgt] = kstest2(Y_test, Y_pred_weighted);

    ks_results(n_var).n_var = n_var;
    ks_results(n_var).KS_lin = ks_lin;
    ks_results(n_var).KS_res = ks_res;
    ks_results(n_var).KS_wgt = ks_wgt;
    ks_results(n_var).p_lin = p_lin;
    ks_results(n_var).p_res = p_res;
    ks_results(n_var).p_wgt = p_wgt;

    if n_var==8
        % 可视化残差 PDF 与 t 分布拟合
        x_range = linspace(min(res_train), max(res_train), 200);
        pdf_empirical = ksdensity(res_train, x_range);
        pdf_t = pdf(pd_t, x_range);

        figure;
        subplot(2,1,1);
        plot(x_range, pdf_empirical, 'k-', 'LineWidth', 2); hold on;
        plot(x_range, pdf_t, 'b--', 'LineWidth', 2);
        legend('Empirical (KDE)', 't Fit');
        xlabel('Residual'); ylabel('Density');
        title('Residual PDF: Empirical vs t Distribution');
        grid on;

        % 可视化观测与三种预测的 PDF
        subplot(2,1,2);
        [f_obs, x_obs] = ksdensity(Y_test);
        [f_lin, x_lin] = ksdensity(Y_pred_lin);
        [f_res, x_res] = ksdensity(Y_pred_res);
        [f_wgt, x_wgt] = ksdensity(Y_pred_weighted);

        plot(x_obs, f_obs, 'k-', 'LineWidth', 2); hold on;
        plot(x_lin, f_lin, 'b--', 'LineWidth', 2);
        plot(x_res, f_res, 'r-.', 'LineWidth', 2);
        plot(x_wgt, f_wgt, 'g:', 'LineWidth', 2);
        legend('Observed', 'Linear', 't Residual', 'Weighted');
        xlabel('dp12'); ylabel('Density');
        title('PDF of Predictions vs Observations');
        grid on;

        figure;hold on
        plot(Y_pred_res(:),Y_test(:),'ok')
        plot([-20 20],[-20 20])
    end
end

% === 转表格 ===
T_lin = struct2table(results_lin);
T_resid = struct2table(results_resid);
T_weighted = struct2table(results_weighted);
T_ks = struct2table(ks_results);


% === 可视化对比 ===
figure;
subplot(4,1,1);
plot(T_lin.n_var, T_lin.R2, '-o', 'LineWidth', 2); hold on;
plot(T_resid.n_var, T_resid.R2, '-s', 'LineWidth', 2);
plot(T_weighted.n_var, T_weighted.R2, '-^', 'LineWidth', 2);
ylabel('R^2'); title('R^2 Comparison'); legend('Linear','t Residual','Weighted'); grid on;

subplot(4,1,2);
plot(T_lin.n_var, T_lin.RMSE, '-o', 'LineWidth', 2); hold on;
plot(T_resid.n_var, T_resid.RMSE, '-s', 'LineWidth', 2);
plot(T_weighted.n_var, T_weighted.RMSE, '-^', 'LineWidth', 2);
ylabel('RMSE'); legend('Linear','t Residual','Weighted'); grid on;

subplot(4,1,3);
plot(T_lin.n_var, T_lin.MAE, '-o', 'LineWidth', 2); hold on;
plot(T_resid.n_var, T_resid.MAE, '-s', 'LineWidth', 2);
plot(T_weighted.n_var, T_weighted.MAE, '-^', 'LineWidth', 2);
ylabel('MAE'); legend('Linear','t Residual','Weighted'); grid on;

subplot(4,1,4);
plot(T_lin.n_var, T_lin.SignAcc, '-o', 'LineWidth', 2); hold on;
plot(T_resid.n_var, T_resid.SignAcc, '-s', 'LineWidth', 2);
plot(T_weighted.n_var, T_weighted.SignAcc, '-^', 'LineWidth', 2);
ylabel('Sign Accuracy'); xlabel('Number of Predictors');
legend('Linear','t Residual','Weighted'); grid on;

sgtitle('Performance Comparison: Linear vs Residual vs Weighted Modeling');

figure
hold on
% plot(1:max_var, rmse_train_all(1:max_var), '-ok', 'LineWidth', 4,'markersize',8, 'markerfacecolor','k');  % 新增部分
plot(T_weighted.n_var, T_weighted.RMSE, '-ok', 'LineWidth', 4,'markersize',8, 'markerfacecolor','k');
ylabel('RMSE (hPa)'); xlabel('Number of Predictors');
% legend('Training', 'Testing');
grid on;
xlim([1 8])
xticks(1:10)
set(gca, 'fontsize', 20);
exportgraphics(gca, ['rmse_predictor.jpg'], 'Resolution', 1100);

figure
plot(1:max_var, rmse_cv, '-ok', 'LineWidth', 4);
ylabel('10-Fold CV RMSE (hPa)');
xlabel('Number of Predictors');
xlim([1 max_var])
xticks(1:max_var)
grid on;
set(gca, 'fontsize', 20);
title('Cross-Validation RMSE vs Number of Predictors');

% === KS 距离随 predictor 数变化 ===
figure;
plot(T_ks.n_var, T_ks.KS_lin, '-o', 'LineWidth', 2); hold on;
plot(T_ks.n_var, T_ks.KS_res, '-s', 'LineWidth', 2);
plot(T_ks.n_var, T_ks.KS_wgt, '-^', 'LineWidth', 2);
xlabel('Number of Predictors'); ylabel('K-S Distance');
legend('Linear', 't Residual', 'Weighted');
title('K-S Distance vs Number of Predictors');
grid on;


%%
best_n = 8;
vars_best = selected_vars(1:best_n);  
vars_best([4;7])=[];
% 构造训练 X（重新构造，以确保干净）
X_train_best = zeros(numel(valid_idx), best_n-2);
for i = 1:best_n-2
    idx = find(strcmp(varNames, vars_best{i}));
    Xi = reshape(X_vars{idx}(:,train_idx), [], 1);
    X_train_best(:,i) = Xi(valid_idx);
end

% 标准化
Xm_best = mean(X_train_best);
Xs_best = std(X_train_best);
Xs_best(Xs_best == 0) = 1;
X_train_best_std = (X_train_best - Xm_best) ./ Xs_best;

% 加权训练
Y_train_centered = Y_train - mean(Y_train);
w_train_best = 1 + abs(Y_train_centered);
mdl_weighted_best = fitlm(X_train_best_std, Y_train, 'Weights', w_train_best, 'VarNames', [vars_best', {'dp12'}]);

% 残差拟合 t 分布
res_train_best = Y_train - predict(mdl_weighted_best, X_train_best_std);
pd_t_best = fitdist(res_train_best, 'tLocationScale');
res_coef = 0.1;

% 保存模型
save('AR_trainingset6sel_v11.mat', ...
    'mdl_weighted_best', 'Xm_best', 'Xs_best', 'vars_best', ...
    'pd_t_best', 'res_coef', 'train_idx','test_idx');










