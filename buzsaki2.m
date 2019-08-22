%% buzsaki - Achilles_11012013 (circular 1D track)
clear
clc

%% set main dir
% main_dir = 'C:\Tamir\work\Courses\MBL\project\figures\buzsaki\opt_UVW_pos_R=1__run2';
% main_dir = 'C:\Tamir\work\Courses\MBL\project\figures\buzsaki\opt_UVW_WX=0_pos_R=1';
% main_dir = 'C:\Tamir\work\Courses\MBL\project\figures\buzsaki\opt_UVW_UV=0_pos_R=0';
main_dir = 'C:\Tamir\work\Courses\MBL\project\figures\buzsaki\opt_UVW_pos_R=2';
mkdir(main_dir)

%% 
set(0,'defaulttextfontsize',15);
set(0,'defaultaxesfontsize',15);

%% load data
load('C:\Tamir\work\Courses\MBL\project\Datasets\buzsaki\Achilles_11012013.tar\Achilles_11012013\Achilles_11012013_sessInfo.mat')
pos = sessInfo.Position;
spikes = sessInfo.Spikes;
epochs = sessInfo.Epochs;

dt = median(diff(pos.TimeStamps));
fs = 1/dt;

%% explore position
figure
% subplot(121)
plot(pos.TwoDLocation(:,1), pos.TwoDLocation(:,2), '.')
axis equal
% subplot(122)
figure
hold on
plot(pos.TimeStamps, pos.OneDLocation, '.');
xline(epochs.MazeEpoch(1));
xline(epochs.MazeEpoch(2));

%% explore spikes
figure
hold on
cell_IDs = spikes.PyrIDs;
% cell_IDs = spikes.IntIDs;
% cell_IDs = [spikes.PyrIDs;spikes.IntIDs];
for ii_cell = 1:length(cell_IDs)
    cellID = cell_IDs(ii_cell);
    cell = struct();
    cell.ID = cellID;
    cell.ts = spikes.SpikeTimes(spikes.SpikeIDs==cellID);
    cell.pos1D = interp1(pos.TimeStamps, pos.OneDLocation, cell.ts, 'linear','extrap');
%     cell.pos1D(isnan(cell.pos1D )) = -1;
    cla
    plot(pos.TimeStamps, pos.OneDLocation, '.k');
    plot(cell.ts, cell.pos1D, '.r');
    title(sprintf('cell %d(%d)',ii_cell,cell.ID));
    pause
end

%% add velocity to position
smooth_p = 1e-2;
pp = csaps(pos.TimeStamps, pos.TwoDLocation', smooth_p);
% pp = csaps(pos.TimeStamps, pos.TwoDLocation');
pp1 = fnder(pp,1);
vel2D_csaps = vecnorm(fnval(pp1, pos.TimeStamps).*fs);
vel_diff = vecnorm(diff(pos.TwoDLocation)'.*fs);

smooth_p = 0.9;
pp = csaps(pos.TimeStamps, pos.OneDLocation', smooth_p);
pp1 = fnder(pp,1);
vel1D_csaps = fnval(pp1,pos.TimeStamps).*fs;

figure
% hold on
subplot(311)
plot(abs(vel1D_csaps))
title('vel1D_csaps')
subplot(312)
plot(vel2D_csaps)
title('vel2D_csaps')
subplot(313)
plot(vel_diff)
title('vel diff')
% legend({'1Dcsaps';'2Dcsaps';'diff'})
% link axes
hax=findall(gcf,'type','axes');
linkaxes(hax,'x')
pos.vel_2D_csaps = vel2D_csaps;
pos.vel_1D_csaps = vel1D_csaps;
pos.vel_2D_diff = [0 vel_diff];

%%
figure
subplot(211)
plot(pos.vel_1D_csaps)
subplot(212)
plot(pos.OneDLocation)

% link axes
hax=findall(gcf,'type','axes');
linkaxes(hax,'x')

%% prepare data (position)
dt_new = .5;
fs_new = 1/dt_new;
pos.t_new = pos.TimeStamps(1) : dt_new : pos.TimeStamps(end);
pos.pos1D_new = interp1(pos.TimeStamps, pos.OneDLocation, pos.t_new, 'linear', 'extrap');
pos.vel_new = interp1(pos.TimeStamps, pos.vel_2D_diff, pos.t_new, 'linear', 'extrap');
nBinPos = 20; % number of spatial positions
T = length(pos.t_new);
X = zeros(nBinPos,T);
[~,EDGES,BIN] = histcounts(pos.pos1D_new,nBinPos);
invalid_IX = find(BIN==0);
BIN(invalid_IX) = 1;
IX = sub2ind(size(X), BIN, 1:T);
X(IX) = 1;
X(:,invalid_IX) = [];
t = pos.t_new;
t(invalid_IX) = [];
vel = interp1(pos.TimeStamps,pos.vel_2D_diff,t, 'linear','extrap');

%% prepare data (spikes)
cell_IDs = spikes.PyrIDs;
% cell_IDs = spikes.IntIDs;
% cell_IDs = [spikes.PyrIDs;spikes.IntIDs];
N = length(cell_IDs);
P = size(X,1);
Y = zeros(N,T);
EDGES = [pos.t_new-dt_new/2 pos.TimeStamps(end)+dt_new/2];
for ii_cell = 1:N
    cellID = cell_IDs(ii_cell);
    cell = struct();
    cell.ts = spikes.SpikeTimes(spikes.SpikeIDs==cellID);
    cell.pos1D = interp1(pos.TimeStamps, pos.OneDLocation, cell.ts, 'linear','extrap');
    cell.invalid_IX = find(isnan(cell.pos1D));
    cell.ts(cell.invalid_IX) = [];
    cell.pos1D(cell.invalid_IX) = [];
    Y(ii_cell,:) = histcounts(cell.ts,EDGES);
end
Y(:,invalid_IX) = [];
N = size(Y,1);
T = size(Y,2);
P = size(X,1);

%% plot X/Y
figure
subplot(211)
imagesc(X)
title('X')
xlabel('time')
ylabel('pos')
colorbar
subplot(212)
imagesc(Y)
title('Y')
xlabel('time')
ylabel('cell')
colorbar
% link axes
hax=findall(gcf,'type','axes');
linkaxes(hax,'x')

%% calculate naive FR maps
time_spent = sum(X');
[r,c] = find(X);
FR_all_naive = [];
for ii_cell = 1:N
    spike_counts = accumarray(r,Y(ii_cell,:),[],@sum);
    FR_all_naive(ii_cell,:) = spike_counts ./ time_spent';
end
FR_all_naive = fs_new .* FR_all_naive; % Hz

% now remove stationary times
X_move = X;
Y_move = Y;
vel_thr = 0.2;
static_IX = find(vel<vel_thr);
X_move(:,static_IX) = [];
Y_move(:,static_IX) = [];

time_spent_move = sum(X_move');
[r,c] = find(X_move);
FR_all_naive_move = [];
for ii_cell = 1:N
    spike_counts = accumarray(r,Y_move(ii_cell,:),[],@sum);
    FR_all_naive_move(ii_cell,:) = spike_counts ./ time_spent_move';
end
FR_all_naive_move = fs_new .* FR_all_naive_move; % Hz

figure('Units','normalized','Position',[0 0 1 1]);
subplot(121)
imagesc(FR_all_naive)
colorbar
xlabel('position')
ylabel('cell')
title('all points');
subplot(122)
imagesc(FR_all_naive_move)
colorbar
xlabel('position')
ylabel('cell')
title('only during movement');

suptitle('Naive FR maps')
fileout = fullfile(main_dir,'naive_FR_maps');
saveas(gcf,fileout,'tif');
saveas(gcf,fileout,'fig');

%% run the model !!!
rng(0);
R = 0;
U = rand(N,R);
V = randn(R,T);
W = rand(N,P);
% W(:,end) = 0; % zero base line
% U = U - mean(U(:));
% V = V - mean(V(:));
nIter = 15;
tic
options = optimoptions('fminunc','SpecifyObjectiveGradient',true);
% U(:)=0;
% V(:)=0;
% W(:)=0;
% X(:)=0;
% V = V2;
% U = U2;

% opttime2=nan(2,nIter);
% fval2=nan(2,nIter);
% for iter=1:nIter
%     disp(iter);
%     UW = UW_join(U,W);
%     [UW fval2(1,iter)] = fminunc(@(UW)optUW(Y,UW,V,X),UW,options); opttime2(1,iter)=toc;
%     [U,W] = UW_split(UW,R,P);
%     [V fval2(2,iter)] = fminunc(@(V)optV(Y,U,V,W,X),V,options); opttime2(2,iter)=toc;
% end

tic
opttime=nan(2,nIter);
fval=nan(2,nIter);
UVW = UVW_join(U,V,W);
for iter=1:nIter
    disp(iter);
%     UW = UW_join(U,W);
%     [UW fval(1,iter)] = fminunc(@(UW)optUW(Y,UW,V,X),UW,options); opttime(1,iter)=toc;
%     [U,W] = UW_split(UW,R,P);
%     [V fval(2,iter)] = fminunc(@(V)optV(Y,U,V,W,X),V,options); opttime(2,iter)=toc;

    [UVW fval(1,iter)] = fminunc(@(UVW)optUVW(Y,UVW,X,R,P,N,T),UVW,options); opttime(1,iter)=toc;
end
[U,V,W] = UVW_split(UVW,R,P,N,T);
toc

Yhat = U*V + W*X;
log_likelihood_final = sum(Yhat.*Y-exp(Yhat),'all');

%% save results!
fileout = fullfile(main_dir,'results');
save(fileout);

%% plot optimization progress
figure
subplot(121)
plot(opttime'./60,'o-');
xlabel('steps')
ylabel('minutes')
ylimits = get(gca,'ylim');
ylimits(1) = 0;
set(gca,'ylim', ylimits);
subplot(122)
plot(fval','o-');
xlabel('steps')
ylabel('loss')
% legend({'UW';'V'})
fileout = fullfile(main_dir,'opt_progress');
saveas(gcf,fileout,'tif');
saveas(gcf,fileout,'fig');


%% plot V / UV / WX / FR / Y
sort_by = 'cell_ID';
% sort_by = 'U';
switch sort_by 
    case 'cell_ID'
        sort_IX = 1:N;
    case 'U'
        [~,sort_IX] = sort(U);
end
figure('Units','normalized','Position',[0 0 1 1]);
subplot(511)
yyaxis left
% plot(V(1,:))
yyaxis right
plot(vel)
hl=legend({'V';'velocity'},'Location','northeastoutside');
hax = gca;
hl.Position = [hax.Position(1)+hax.Position(3)-0.05 hax.Position(2)+hax.Position(4)+0.02 0.02 0.2*hax.Position(4)];
axis tight
title('V vs. velocity')
% xlabel('time')

subplot(512)
sdf=U*V;
imagesc(sdf(sort_IX,:))
title('UV')
% xlabel('time')
ylabel('cells')
hax= gca;
hcl=colorbar('eastoutside');
hcl.Position = [hax.Position(1)+hax.Position(3)+0.02 hax.Position(2) 0.02 hax.Position(4)];

subplot(513)
sdf = W*X;
imagesc(sdf(sort_IX,:))
title('WX')
% xlabel('time')
ylabel('cells')
hax= gca;
hcl=colorbar('eastoutside');
hcl.Position = [hax.Position(1)+hax.Position(3)+0.02 hax.Position(2) 0.02 hax.Position(4)];

subplot(514)
sdf = fs_new*exp(U*V+W*X);
imagesc(sdf(sort_IX,:))
title('estimate FR')
% xlabel('time')
ylabel('cells')
hax= gca;
hcl=colorbar('eastoutside');
hcl.Position = [hax.Position(1)+hax.Position(3)+0.02 hax.Position(2) 0.02 hax.Position(4)];

subplot(515)
sdf = fs_new*Y;
imagesc(sdf(sort_IX,:))
title('Y')
xlabel('time')
ylabel('cells')
hax= gca;
hcl=colorbar('eastoutside');
hcl.Position = [hax.Position(1)+hax.Position(3)+0.02 hax.Position(2) 0.02 hax.Position(4)];

% link axes
hax=findall(gcf,'type','axes');
linkaxes(hax,'x')

fileout = fullfile(main_dir,sprintf('components_breakdown_sort_by_%s',sort_by));
saveas(gcf,fileout,'tif');
saveas(gcf,fileout,'fig');

xlim([2320 2750])
fileout = fullfile(main_dir,sprintf('components_breakdown_sort_by_%s_ZOOM_IN',sort_by));
saveas(gcf,fileout,'tif');
saveas(gcf,fileout,'fig');

%% plot V vs. velocity
figure
plot(vel,V(1,:),'.')
lm = fitlm(vel,V(1,:));
xlabel('velocity')
ylabel('V')
xlim([0 1.5])
text(0.8,0.95,sprintf('R = %.2f',lm.Rsquared.Ordinary),'Units','normalized');
fileout = fullfile(main_dir,'V_vs_speed');
saveas(gcf,fileout,'tif');
saveas(gcf,fileout,'fig');

%% compare W to FR map (single cells)
% cells_list = [2 33 45 47 53 87];
% ex_IX = 6;
% cell_num = cells_list(ex_IX);
mkdir(fullfile(main_dir,'cells_FR'))
figure
hold on
FR_all_W = [];
for cell_num = 1:N
    cell_ID = cell_IDs(cell_num);
    FR_all_W = fs_new*exp(W);
    cla
    plot(FR_all_naive(cell_num,:))
    plot(FR_all_naive_move(cell_num,:))
    plot(FR_all_W(cell_num,:))
    xlabel('position')
    ylabel('Firing rate (Hz)')
    legend({'naive FR (all points)';'naive FR (speed thr)';'W estimation'},...
        'Location','best')
    title(sprintf('cell %d (%d)',cell_num,cell_ID))
    fileout = fullfile(main_dir,'cells_FR',sprintf('cell_%d(%d)',cell_num,cell_ID));
    saveas(gcf,fileout,'tif');
    saveas(gcf,fileout,'fig');
end

%% calculate SI
prob = time_spent./sum(time_spent);
r_i = FR_all_naive;
r_i = r_i - min(r_i,[],2);
meanFR = sum(r_i .* prob,2);
SI = sum((prob.*r_i).*log2((r_i+eps)./meanFR),2);
SI_all_naive = SI ./ meanFR; % convert to bit/spike

prob = time_spent_move./sum(time_spent_move);
r_i = FR_all_naive_move;
r_i = r_i - min(r_i,[],2);
meanFR = sum(r_i .* prob,2);
SI = sum((prob.*r_i).*log2((r_i+eps)./meanFR),2);
SI_all_naive_move = SI ./ meanFR; % convert to bit/spike

prob = time_spent./sum(time_spent);
r_i = FR_all_W;
r_i = r_i - min(r_i,[],2);
meanFR = sum(r_i .* prob,2);
SI = sum((prob.*r_i).*log2((r_i+eps)./meanFR),2);
SI_all_W = SI ./ meanFR; % convert to bit/spike

%% compare W to FR map (population)
rho1 = corr(FR_all_W',FR_all_naive');
rho2 = corr(FR_all_W',FR_all_naive_move');
edges = linspace(-1,1,50);
figure('Units','normalized','Position',[0.1 0.2 0.7 0.5])
subplot(131)
hold on
histogram(diag(rho1),edges)
histogram(diag(rho2),edges)
legend({'all points';'speed thr'},'Location','northoutside')
xlabel('correlation')
ylabel('counts')
title('corr W vs. naive')
subplot(132)
hold on
ksdensity(SI_all_naive)
ksdensity(SI_all_naive_move)
ksdensity(SI_all_W)
xlabel('SI')
ylabel('density')
title('spatial information')
legend({'naive (all points)';'naive (speed thr)';'W'},'Location','northoutside')
subplot(133)
hold on
ksdensity(max(FR_all_naive))
ksdensity(max(FR_all_naive_move))
ksdensity(max(FR_all_W))
xlabel('max FR (Hz)')
ylabel('density')
title('Max FR')
legend({'naive (all points)';'naive (speed thr)';'W'},'Location','northoutside')
fileout = fullfile(main_dir,'compare_FR_W_vs_naive_population');
saveas(gcf,fileout,'tif');
saveas(gcf,fileout,'fig');

%% plot W
WW = Y/X;
figure
imagesc(exp(W))
title('W');
colorbar
figure
imagesc(WW)
title('WW');
colorbar

%% UV SVD
[UU,S,VV] = svd(U*V);
figure
imagesc(UU)
colorbar
figure
imagesc(VV)
colorbar

%%
figure
yyaxis left
hold on
plot(VV(:,1),'b')
plot(VV(:,2),'r')
yyaxis right
plot(vel,'k')

%% plot xcorr between V1/V2/velocity
figure
subplot(221)
[c,lags] = xcorr(VV(:,1),vel);
plot(lags.*dt_new,c)
title('VV1 vs. speed')

subplot(222)
[c,lags] = xcorr(VV(:,2),vel);
plot(lags.*dt_new,c)
title('VV2 vs. speed')

subplot(223)
[c,lags] = xcorr(VV(:,1),VV(:,2));
plot(lags.*dt_new,c)
title('VV1 vs. VV2')

% link axes
hax=findall(gcf,'type','axes');
linkaxes(hax,'x')
xlim([-20 20])

%% optimization functions
function [f,g] = optW(Y,U,V,W,X)
    Yhat = U*V + W*X;
    f = -sum(Yhat.*Y-exp(Yhat),'all');
    g = -(Y-exp(Yhat))*X';
    % add L2 regularization for W
    lambda = 1;
    f = f + lambda*norm(W,'fro')^2;
    g = g + 2*lambda*W;
end
function [f,g] = optU(Y,U,V,W,X)
    Yhat = U*V + W*X;
    f = -sum(Yhat.*Y-exp(Yhat),'all');
    g = -(Y-exp(Yhat))*V';
    % add L2 regularization for U
    lambda = 1;
    f = f + lambda*norm(U,'fro')^2;
    g = g + 2*lambda*U;
end
function [f,g] = optV(Y,U,V,W,X)
    Yhat = U*V + W*X;
    f = -sum(Yhat.*Y-exp(Yhat),'all');
    g = -U'*(Y-exp(Yhat));
    % add smoothness regularization on V
    Vpad = padarray(V,[0 1],'replicate','both');
    lambda = 1e3;
    f = f + lambda*sum(diff(Vpad,1,2).^2,'all');
    g = g + 2*lambda.*( 2*V -Vpad(:,1:end-2) -Vpad(:,3:end) );
    % add L2 regularization for V
    lambda = 1;
    f = f + lambda*norm(V,'fro')^2;
    g = g + 2*lambda*V;
end
function [f,g] = optUW(Y,UW,V,X)
    R = size(V,1);
    P = size(X,1);
    [U,W] = UW_split(UW,R,P);
    Yhat = U*V + W*X;
    f = -sum(Yhat.*Y-exp(Yhat),'all');
    gW = -(Y-exp(Yhat))*X';
    gU = -(Y-exp(Yhat))*V';
    % add L2 regularization for U
    lambda = 1;
    f = f + lambda*norm(U,'fro')^2;
    gU = gU + 2*lambda*U;
    % add L2 regularization for W
    lambda = 1;
    f = f + lambda*norm(U,'fro')^2;
    gW = gW + 2*lambda*W;
    g = UW_join(gU,gW);
end
function UW = UW_join(U,W)
    UW = [U W];
end
function [U,W] = UW_split(UW,R,P)
    U = UW(:,1:R);
    W = UW(:,(R+1):end);
end
function [f,g] = optUVW(Y,UVW,X,R,P,N,T)
    [U,V,W] = UVW_split(UVW,R,P,N,T);
    Yhat = U*V + W*X;
    f = -sum(Yhat.*Y-exp(Yhat),'all');
    gW = -(Y-exp(Yhat))*X';
    gU = -(Y-exp(Yhat))*V';
    gV = -U'*(Y-exp(Yhat));
    % add L2 regularization for U
    lambda = 1;
    f = f + lambda*norm(U,'fro')^2;
    gU = gU + 2*lambda*U;
    % add L2 regularization for W
    lambda = 1;
    f = f + lambda*norm(U,'fro')^2;
    gW = gW + 2*lambda*W;
    % add smoothness regularization on V
    Vpad = padarray(V,[0 1],'replicate','both');
    lambda = 1e3;
    f = f + lambda*sum(diff(Vpad,1,2).^2,'all');
    gV = gV + 2*lambda.*( 2*V -Vpad(:,1:end-2) -Vpad(:,3:end) );
    % add L2 regularization for V
    lambda = 1;
    f = f + lambda*norm(V,'fro')^2;
    gV = gV + 2*lambda*V;
    
    g = UVW_join(gU,gV,gW);
end
function UVW = UVW_join(U,V,W)
    UVW = [U(:);V(:);W(:)];
end
function [U,V,W] = UVW_split(UVW,R,P,N,T)
    U = UVW(              1:(N*R)  );
    V = UVW( (N*R)     + (1:(R*T)) );
    W = UVW( (N*R+R*T) + (1:(N*P)) );
    U = reshape(U,[N R]);
    V = reshape(V,[R T]);
    W = reshape(W,[N P]);
end



%%
